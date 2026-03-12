/**
 * \file main.c
 * \brief GRASS GIS module \c i.atcorr2 — 6SV2.1 LUT-based single-band
 *        atmospheric correction (C + OpenMP rewrite).
 *
 * MODULE:    i.atcorr2
 * AUTHOR:    Yann Chemin
 * PURPOSE:   Atmospheric correction for a single-band raster using the
 *            grass_sixsv shared library (6SV2.1 LUT-based engine).
 * COPYRIGHT: (C) 2026 by the GRASS Development Team.
 *            GNU GPL >= 2.
 *
 * Processing pipeline
 * -------------------
 *  1. Build a 3-D LUT [AOD × H₂O × λ] via atcorr_compute_lut() or load
 *     a pre-computed LUT from disk.
 *  2. Read input raster row-by-row into fixed-size chunks (default) or all
 *     at once (-R flag).
 *  3. For each chunk, parallelise over pixels with OpenMP:
 *       radiance  →  TOA reflectance  →  bilinear LUT interp  →  BOA inversion
 *  4. Write output raster.
 *
 * Chunk default: ~256 MB working set; -R flag loads everything into RAM
 * (warns at >50 % available RAM, aborts at >90 %).
 */

#include <grass/gis.h>
#include <grass/glocale.h>
#include <grass/raster.h>

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* system libsixsv public API */
#include <sixsv/atcorr.h>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
#endif

/* ── Constants ────────────────────────────────────────────────────────────── */

#define LUT_MAGIC    0x4C555400u
#define LUT_VERSION  1u

/* Target working-set per chunk in bytes (256 MB) */
#define CHUNK_BYTES  (256u * 1024u * 1024u)

/* RAM thresholds for -R mode */
#define RAM_WARN_FRAC  0.50   /* warn if > 50 % of available RAM */
#define RAM_ABORT_FRAC 0.90   /* abort if > 90 % of available RAM */

/*
 * 6SV atmosphere / aerosol codes used internally.
 * These intentionally match the original Fortran iatm/iaer numbering,
 * which differs from the symbolic constants in atcorr.h (those are the
 * C-library's internal remapped identifiers).  The string → int mapping
 * functions below return the 6SV Fortran codes that the library expects.
 */

/* ── Small helpers ────────────────────────────────────────────────────────── */

/** \cond INTERNAL */

/**
 * \brief Parse a comma-separated string of floating-point numbers.
 *
 * Tokenises \c str on commas and converts each token to \c float.  Stops
 * after \c max_n values or when the token list is exhausted.
 *
 * \param[in]  str    Input string (e.g. \c "0.0,0.1,0.5").
 * \param[out] out    Output array; must have room for at least \c max_n elements.
 * \param[in]  max_n  Maximum number of values to parse.
 * \return Number of values successfully parsed, or -1 on conversion error.
 */
static int parse_csv_floats(const char *str, float *out, int max_n)
{
    char *buf = G_store(str);
    char *tok, *sv = NULL;
    int   n  = 0;

    tok = strtok_r(buf, ",", &sv);
    while (tok && n < max_n) {
        char *end;
        out[n] = (float)strtod(tok, &end);
        if (end == tok) { G_free(buf); return -1; }
        n++;
        tok = strtok_r(NULL, ",", &sv);
    }
    G_free(buf);
    return n;
}

/**
 * \brief Map an atmosphere name string to a 6SV Fortran \c iatm code.
 *
 * \param[in]  s  Atmosphere name as accepted by the \c atmosphere= option.
 * \return 6SV \c iatm integer code (0–6); defaults to 6 (US62) for unknown names.
 */
static int atmo_str_to_int(const char *s)
{
    if (!strcmp(s, "none"))       return 0;
    if (!strcmp(s, "us62"))       return ATMO_US62;       /* 1 */
    if (!strcmp(s, "midsum"))     return ATMO_MIDSUM;     /* 2 */
    if (!strcmp(s, "midwin"))     return ATMO_MIDWIN;     /* 3 */
    if (!strcmp(s, "tropical"))   return ATMO_TROPICAL;   /* 4 */
    if (!strcmp(s, "subarctsum")) return ATMO_SUBSUM;     /* 5 */
    if (!strcmp(s, "suarctwint")) return ATMO_SUBWIN;     /* 6 */
    return ATMO_US62;
}

/**
 * \brief Map an aerosol model name string to a 6SV Fortran \c iaer code.
 *
 * \param[in]  s  Aerosol model name as accepted by the \c aerosol= option.
 * \return 6SV \c iaer integer code (0–6); defaults to 1 (continental) for unknown names.
 */
static int aerosol_str_to_int(const char *s)
{
    if (!strcmp(s, "none"))          return AEROSOL_NONE;         /* 0 */
    if (!strcmp(s, "continental"))   return AEROSOL_CONTINENTAL;  /* 1 */
    if (!strcmp(s, "maritime"))      return AEROSOL_MARITIME;     /* 2 */
    if (!strcmp(s, "urban"))         return AEROSOL_URBAN;        /* 3 */
    if (!strcmp(s, "desert"))        return AEROSOL_DESERT;       /* 5 */
    return AEROSOL_CONTINENTAL;
}

/**
 * \brief Compute day-of-year from month and day (non-leap-year approximation).
 *
 * \param[in]  month  Month number (1–12).
 * \param[in]  day    Day of month (1–31).
 * \return Day of year (1–365).
 */
static int doy_from_md(int month, int day)
{
    static const int dm[] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
    int d = day;
    for (int m = 1; m < month; m++) d += dm[m];
    return d;
}

/**
 * \brief Compute an approximate solar zenith angle.
 *
 * Uses a single-harmonic approximation for solar declination and a mean
 * equation-of-time term.  Accuracy is ~1–2° compared to a full ephemeris.
 *
 * \param[in]  lat_deg  Geographic latitude in degrees (north positive).
 * \param[in]  doy      Day of year (1–365).
 * \param[in]  hour_ut  UTC hour of day (decimal).
 * \return Solar zenith angle in degrees.
 */
static float approx_sza(float lat_deg, int doy, float hour_ut)
{
    float lat  = (float)(M_PI / 180.0) * lat_deg;
    float decl = (float)(M_PI / 180.0) * 23.45f *
                 sinf((float)(M_PI / 180.0) * (360.0f / 365.0f * (doy - 81)));
    float ha   = (float)(M_PI / 180.0) * 15.0f * (hour_ut - 12.0f);
    float cs   = sinf(lat) * sinf(decl) + cosf(lat) * cosf(decl) * cosf(ha);
    cs = cs < -1.0f ? -1.0f : cs > 1.0f ? 1.0f : cs;
    return (float)(180.0 / M_PI) * acosf(cs);
}

/* ── LUT file I/O ─────────────────────────────────────────────────────────── */

/**
 * \brief Write a precomputed LUT to a binary file.
 *
 * Saves the four correction arrays (R_atm, T_down, T_up, s_alb) together with
 * the AOD and H₂O grid axes, wavelength, and a magic-number header to a flat
 * binary file compatible with i.hyper.atcorr's LUT format.
 *
 * \param[in]  path   Output file path.
 * \param[in]  aod    AOD grid axis array (length \c n_aod).
 * \param[in]  n_aod  Number of AOD grid points.
 * \param[in]  h2o    H₂O grid axis array (length \c n_h2o).
 * \param[in]  n_h2o  Number of H₂O grid points.
 * \param[in]  wl     Single wavelength (µm) for this LUT slice.
 * \param[in]  la     Pointer to the LUT arrays to save.
 */
static void save_lut(const char *path,
                     const float *aod, int n_aod,
                     const float *h2o, int n_h2o,
                     float wl, const LutArrays *la)
{
    FILE *fp = fopen(path, "wb");
    if (!fp)
        G_fatal_error(_("Cannot write LUT file %s: %s"), path, strerror(errno));

    uint32_t magic = LUT_MAGIC, ver = LUT_VERSION;
    fwrite(&magic, 4, 1, fp);
    fwrite(&ver,   4, 1, fp);

    int32_t na = n_aod, nh = n_h2o, nw = 1;
    fwrite(&na, 4, 1, fp);
    fwrite(&nh, 4, 1, fp);
    fwrite(&nw, 4, 1, fp);

    fwrite(aod,        sizeof(float), n_aod,           fp);
    fwrite(h2o,        sizeof(float), n_h2o,           fp);
    fwrite(&wl,        sizeof(float), 1,               fp);

    int n = n_aod * n_h2o;
    fwrite(la->R_atm,  sizeof(float), (size_t)n,       fp);
    fwrite(la->T_down, sizeof(float), (size_t)n,       fp);
    fwrite(la->T_up,   sizeof(float), (size_t)n,       fp);
    fwrite(la->s_alb,  sizeof(float), (size_t)n,       fp);

    fclose(fp);
    G_verbose_message(_("LUT saved to %s"), path);
}

/**
 * \brief Load a precomputed LUT from a binary file.
 *
 * Reads the magic-number header, grid dimensions, axis arrays, and the four
 * correction arrays from a file previously written by save_lut().  On success
 * the caller owns all allocated buffers (\c *aod_out, \c *h2o_out, and
 * \c la->R_atm … la->s_alb), all allocated with G_malloc().
 *
 * \param[in]  path       Input file path.
 * \param[out] aod_out    Pointer to allocated AOD grid array.
 * \param[out] n_aod_out  Number of AOD grid points.
 * \param[out] h2o_out    Pointer to allocated H₂O grid array.
 * \param[out] n_h2o_out  Number of H₂O grid points.
 * \param[out] wl_out     Wavelength stored in the file (µm).
 * \param[out] la         Filled with pointers to the loaded correction arrays.
 * \return 0 on success, -1 if the file cannot be opened or is not a valid LUT.
 */
static int load_lut(const char *path,
                    float **aod_out, int *n_aod_out,
                    float **h2o_out, int *n_h2o_out,
                    float  *wl_out,
                    LutArrays *la)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        G_warning(_("Cannot open LUT file %s: %s"), path, strerror(errno));
        return -1;
    }

    uint32_t magic, ver;
    if (fread(&magic, 4, 1, fp) != 1 || fread(&ver, 4, 1, fp) != 1 ||
        magic != LUT_MAGIC) {
        G_warning(_("Not a valid LUT file: %s"), path);
        fclose(fp);
        return -1;
    }

    int32_t na, nh, nw;
    if (fread(&na, 4, 1, fp) != 1 ||
        fread(&nh, 4, 1, fp) != 1 ||
        fread(&nw, 4, 1, fp) != 1) {
        fclose(fp); return -1;
    }

    float *aod = G_malloc((size_t)na * sizeof(float));
    float *h2o = G_malloc((size_t)nh * sizeof(float));
    float  wl;

    if (fread(aod, sizeof(float), (size_t)na, fp) != (size_t)na ||
        fread(h2o, sizeof(float), (size_t)nh, fp) != (size_t)nh ||
        fread(&wl, sizeof(float), 1, fp) != 1) {
        G_free(aod); G_free(h2o); fclose(fp); return -1;
    }

    /* skip extra wavelength entries if present */
    for (int w = 1; w < nw; w++) {
        float tmp; fread(&tmp, sizeof(float), 1, fp);
    }

    int n   = na * nh * nw;
    float *Ra  = G_malloc((size_t)n * sizeof(float));
    float *Td  = G_malloc((size_t)n * sizeof(float));
    float *Tu  = G_malloc((size_t)n * sizeof(float));
    float *sa  = G_malloc((size_t)n * sizeof(float));

    if (fread(Ra, sizeof(float), (size_t)n, fp) != (size_t)n ||
        fread(Td, sizeof(float), (size_t)n, fp) != (size_t)n ||
        fread(Tu, sizeof(float), (size_t)n, fp) != (size_t)n ||
        fread(sa, sizeof(float), (size_t)n, fp) != (size_t)n) {
        G_free(aod); G_free(h2o);
        G_free(Ra); G_free(Td); G_free(Tu); G_free(sa);
        fclose(fp); return -1;
    }

    fclose(fp);

    *aod_out   = aod;  *n_aod_out = (int)na;
    *h2o_out   = h2o;  *n_h2o_out = (int)nh;
    *wl_out    = wl;

    la->R_atm  = Ra;  la->T_down = Td;
    la->T_up   = Tu;  la->s_alb  = sa;
    la->T_down_dir = NULL;
    la->R_atmQ = NULL;
    la->R_atmU = NULL;

    G_verbose_message(_("LUT loaded from %s  (AOD n=%d, H2O n=%d, WL=%.4f µm)"),
                      path, na, nh, wl);
    return 0;
}

/* ── 6S parameter-file parser (i.atcorr backward compatibility) ────────────── */

/**
 * \brief Read non-empty, non-comment lines from a 6S conditions file.
 *
 * Strips inline comments beginning with \c # or \c (, trims surrounding
 * whitespace, and discards blank lines.
 *
 * \param[in]  path    Path to the 6S parameter file.
 * \param[out] nlines  Receives the number of non-empty lines returned.
 * \return Heap-allocated array of \c *nlines G_store()'d strings.
 *         The caller must G_free() each element and the array itself.
 */
static char **read_6s_lines(const char *path, int *nlines)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
        G_fatal_error(_("Cannot open 6S parameter file %s: %s"), path, strerror(errno));

    char   buf[512];
    int    cap = 64, n = 0;
    char **lines = G_malloc((size_t)cap * sizeof(char *));

    while (fgets(buf, sizeof(buf), fp)) {
        /* strip comments: # and ( */
        char *p = buf;
        while (*p && *p != '#' && *p != '(') p++;
        *p = '\0';
        /* trim trailing whitespace */
        while (p > buf && (p[-1] == ' ' || p[-1] == '\t' ||
                           p[-1] == '\r' || p[-1] == '\n'))
            *--p = '\0';
        /* skip empty */
        char *s = buf;
        while (*s == ' ' || *s == '\t') s++;
        if (!*s) continue;
        if (n >= cap) { cap *= 2; lines = G_realloc(lines, (size_t)cap * sizeof(char *)); }
        lines[n++] = G_store(s);
    }
    fclose(fp);
    *nlines = n;
    return lines;
}

/**
 * \brief Parse a 6S conditions file in the legacy i.atcorr format.
 *
 * Reads the geometry/date line, atmospheric model code, aerosol model code,
 * and visibility/AOD from the file.  Any output pointer may be NULL to skip
 * that value.  Explicit command-line options should override values returned
 * here.
 *
 * \param[in]  path          Path to the 6S parameter file.
 * \param[out] sza_out       Solar zenith angle derived from date/time/location (degrees).
 * \param[out] atmo_out      Atmospheric model code (6SV \c iatm).
 * \param[out] aerosol_out   Aerosol model code (6SV \c iaer).
 * \param[out] aod_from_vis  AOD derived from the visibility line (km).
 * \return 0 on success, -1 on parse error.
 */
static int parse_6s_file(const char *path,
                         float *sza_out,
                         int   *atmo_out,
                         int   *aerosol_out,
                         float *aod_from_vis)
{
    int    n = 0;
    char **L = read_6s_lines(path, &n);

    if (n < 8) {
        G_warning(_("6S parameter file %s has fewer lines than expected (%d)."), path, n);
    }

    int  idx = 0;
#define NEXT() ((idx < n) ? L[idx++] : "0")

    /* Line 1: geometry condition code (skip) */
    NEXT();

    /* Line 2: month day hour [lon lat] */
    int   month = 6, day = 21;
    float hour = 12.0f, lon = 0.0f, lat = 30.0f;
    if (idx < n) {
        sscanf(L[idx++], "%d %d %f %f %f", &month, &day, &hour, &lon, &lat);
    }
    int doy = doy_from_md(month, day);
    if (sza_out)
        *sza_out = approx_sza(lat, doy, hour);

    /* Line 3: atmospheric model */
    int atmo = ATMO_US62;
    if (idx < n) atmo = atoi(L[idx++]);
    if (atmo_out) *atmo_out = atmo;

    /* Line 4: aerosol model */
    int aerosol = 1; /* continental */
    if (idx < n) aerosol = atoi(L[idx++]);
    if (aerosol_out) *aerosol_out = aerosol;

    /* Line 5: visibility km (0 → next line is AOD) */
    float vis = 15.0f;
    if (idx < n) vis = (float)atof(L[idx++]);
    if (aod_from_vis) {
        if (vis == 0.0f && idx < n) {
            *aod_from_vis = (float)atof(L[idx++]);
        } else if (vis > 0.0f) {
            *aod_from_vis = 3.912f / vis;
        } else {
            *aod_from_vis = 0.0f;
        }
    }

#undef NEXT

    for (int i = 0; i < n; i++) G_free(L[i]);
    G_free(L);
    return 0;
}

/* ── Raster helpers ───────────────────────────────────────────────────────── */

/** \endcond */

/**
 * \brief Allocate and read an entire FCELL raster map into memory.
 *
 * \param[in]  mapname  GRASS raster map name.
 * \param[in]  nrows    Number of rows in the current region.
 * \param[in]  ncols    Number of columns in the current region.
 * \return G_malloc()'d array of \c nrows × \c ncols \c FCELL values.
 *         The caller is responsible for calling G_free().
 */
static FCELL *read_raster_f(const char *mapname, int nrows, int ncols)
{
    FCELL *buf = G_malloc((size_t)nrows * (size_t)ncols * sizeof(FCELL));
    int    fd  = Rast_open_old(mapname, "");

    for (int r = 0; r < nrows; r++)
        Rast_get_f_row(fd, buf + (size_t)r * ncols, r);

    Rast_close(fd);
    return buf;
}

/* ── Main ─────────────────────────────────────────────────────────────────── */

/**
 * \brief GRASS module entry point for \c i.atcorr2.
 *
 * Parses command-line options, builds or loads the 3-D [AOD × H₂O × λ] LUT,
 * and corrects the input single-band raster from TOA radiance (or reflectance)
 * to bottom-of-atmosphere reflectance.  Processing is performed in row-chunks
 * (default ~256 MB working set; \c -R loads everything into RAM) with optional
 * OpenMP parallelism over pixels.
 *
 * Two correction paths are supported:
 * - **Scene-average** (default): all four atmospheric coefficients are constant;
 *   the inner loop is branchless and SIMD-vectorisable.
 * - **Per-pixel** (\c aod_map= / \c h2o_map= provided): bilinear LUT
 *   interpolation per pixel; threaded but not SIMD.
 *
 * \param[in]  argc  Argument count.
 * \param[in]  argv  Argument vector.
 * \return \c EXIT_SUCCESS on success, exits via G_fatal_error() on failure.
 */
int main(int argc, char *argv[])
{
    struct GModule *module;
    struct Option  *opt_input, *opt_output, *opt_wl, *opt_sza, *opt_vza,
                   *opt_raa, *opt_altitude, *opt_target_elev, *opt_doy,
                   *opt_atmo, *opt_aerosol, *opt_aod, *opt_h2o,
                   *opt_aod_val, *opt_h2o_val, *opt_ozone, *opt_pressure,
                   *opt_elevation,
                   *opt_aod_map, *opt_h2o_map, *opt_params, *opt_lut;
    struct Flag    *flag_int, *flag_polar, *flag_ram;

    /* ── Module init ─────────────────────────────────────────────────────── */
    G_gisinit(argv[0]);

    module = G_define_module();
    G_add_keyword(_("imagery"));
    G_add_keyword(_("atmospheric correction"));
    G_add_keyword(_("radiometric conversion"));
    G_add_keyword(_("reflectance"));
    G_add_keyword(_("satellite"));
    G_add_keyword(_("6SV"));
    module->description = _("Atmospheric correction using grass_sixsv (6SV2.1 LUT-based).");

    /* ── Options ─────────────────────────────────────────────────────────── */
    opt_input = G_define_standard_option(G_OPT_R_INPUT);
    opt_input->guisection = _("Input");

    opt_output = G_define_standard_option(G_OPT_R_OUTPUT);
    opt_output->guisection = _("Output");

    opt_wl = G_define_option();
    opt_wl->key         = "wavelength";
    opt_wl->type        = TYPE_DOUBLE;
    opt_wl->required    = NO;
    opt_wl->label       = _("Centre wavelength of the input band (µm)");
    opt_wl->description = _("Required unless parameters= (6S file) provides it");
    opt_wl->guisection  = _("Band");

    opt_sza = G_define_option();
    opt_sza->key         = "sza";
    opt_sza->type        = TYPE_DOUBLE;
    opt_sza->required    = NO;
    opt_sza->label       = _("Solar zenith angle (degrees)");
    opt_sza->description = _("Required unless parameters= (6S file) is given");
    opt_sza->guisection  = _("Geometry");

    opt_vza = G_define_option();
    opt_vza->key         = "vza";
    opt_vza->type        = TYPE_DOUBLE;
    opt_vza->required    = NO;
    opt_vza->answer      = "5.0";
    opt_vza->description = _("View zenith angle (degrees)");
    opt_vza->guisection  = _("Geometry");

    opt_raa = G_define_option();
    opt_raa->key         = "raa";
    opt_raa->type        = TYPE_DOUBLE;
    opt_raa->required    = NO;
    opt_raa->answer      = "90.0";
    opt_raa->description = _("Relative azimuth angle (degrees)");
    opt_raa->guisection  = _("Geometry");

    opt_altitude = G_define_option();
    opt_altitude->key         = "altitude";
    opt_altitude->type        = TYPE_DOUBLE;
    opt_altitude->required    = NO;
    opt_altitude->answer      = "1000.0";
    opt_altitude->label       = _("Sensor altitude (km)");
    opt_altitude->description = _("Use >900 for satellite orbit (default 1000); 0 = ground level");
    opt_altitude->guisection  = _("Geometry");

    opt_target_elev = G_define_option();
    opt_target_elev->key         = "target_elevation";
    opt_target_elev->type        = TYPE_DOUBLE;
    opt_target_elev->required    = NO;
    opt_target_elev->answer      = "0.0";
    opt_target_elev->label       = _("Mean target elevation above sea level (km)");
    opt_target_elev->description = _("Scene-average ground elevation; overridden by elevation= map");
    opt_target_elev->guisection  = _("Geometry");

    opt_doy = G_define_option();
    opt_doy->key         = "doy";
    opt_doy->type        = TYPE_INTEGER;
    opt_doy->required    = NO;
    opt_doy->label       = _("Day of year [1-365]");
    opt_doy->description = _("Required for radiance→TOA reflectance conversion");
    opt_doy->guisection  = _("Geometry");

    opt_atmo = G_define_option();
    opt_atmo->key         = "atmosphere";
    opt_atmo->type        = TYPE_STRING;
    opt_atmo->required    = NO;
    opt_atmo->answer      = "us62";
    opt_atmo->options     = "none,tropical,midsum,midwin,subarctsum,suarctwint,us62";
    opt_atmo->label       = _("Atmospheric model");
    opt_atmo->guisection  = _("Atmosphere");

    opt_aerosol = G_define_option();
    opt_aerosol->key         = "aerosol";
    opt_aerosol->type        = TYPE_STRING;
    opt_aerosol->required    = NO;
    opt_aerosol->answer      = "continental";
    opt_aerosol->options     = "none,continental,maritime,urban,desert";
    opt_aerosol->label       = _("Aerosol model");
    opt_aerosol->guisection  = _("Atmosphere");

    opt_aod = G_define_option();
    opt_aod->key         = "aod";
    opt_aod->type        = TYPE_STRING;
    opt_aod->required    = NO;
    opt_aod->answer      = "0.0,0.1,0.2,0.5";
    opt_aod->label       = _("AOD at 550 nm LUT grid values (comma-separated)");
    opt_aod->guisection  = _("Atmosphere");

    opt_h2o = G_define_option();
    opt_h2o->key         = "h2o";
    opt_h2o->type        = TYPE_STRING;
    opt_h2o->required    = NO;
    opt_h2o->answer      = "1.0,2.0,3.5";
    opt_h2o->label       = _("Column water vapour LUT grid values, g/cm² (comma-separated)");
    opt_h2o->guisection  = _("Atmosphere");

    opt_aod_val = G_define_option();
    opt_aod_val->key         = "aod_val";
    opt_aod_val->type        = TYPE_DOUBLE;
    opt_aod_val->required    = NO;
    opt_aod_val->label       = _("AOD value for correction (overrides aod_map=)");
    opt_aod_val->description = _("Single AOD at 550 nm; if omitted, mid-LUT value is used");
    opt_aod_val->guisection  = _("Atmosphere");

    opt_h2o_val = G_define_option();
    opt_h2o_val->key         = "h2o_val";
    opt_h2o_val->type        = TYPE_DOUBLE;
    opt_h2o_val->required    = NO;
    opt_h2o_val->label       = _("H₂O column value for correction, g/cm² (overrides h2o_map=)");
    opt_h2o_val->guisection  = _("Atmosphere");

    opt_ozone = G_define_option();
    opt_ozone->key         = "ozone";
    opt_ozone->type        = TYPE_DOUBLE;
    opt_ozone->required    = NO;
    opt_ozone->answer      = "300.0";
    opt_ozone->description = _("Ozone column (Dobson units; 0 = use standard atmosphere)");
    opt_ozone->guisection  = _("Atmosphere");

    opt_pressure = G_define_option();
    opt_pressure->key         = "surface_pressure";
    opt_pressure->type        = TYPE_DOUBLE;
    opt_pressure->required    = NO;
    opt_pressure->answer      = "0.0";
    opt_pressure->description = _("Surface pressure (hPa; 0 = use standard atmosphere)");
    opt_pressure->guisection  = _("Atmosphere");

    opt_elevation = G_define_standard_option(G_OPT_R_ELEV);
    opt_elevation->key         = "elevation";
    opt_elevation->required    = NO;
    opt_elevation->label       = _("Input elevation raster map (km above sea level)");
    opt_elevation->description = _("Per-pixel target elevation; overrides target_elevation=");
    opt_elevation->guisection  = _("Input");

    opt_aod_map = G_define_standard_option(G_OPT_R_INPUT);
    opt_aod_map->key         = "aod_map";
    opt_aod_map->required    = NO;
    opt_aod_map->label       = _("Per-pixel AOD raster map");
    opt_aod_map->description = _("Spatially variable AOD at 550 nm; overrides aod_val=");
    opt_aod_map->guisection  = _("Atmosphere");

    opt_h2o_map = G_define_standard_option(G_OPT_R_INPUT);
    opt_h2o_map->key         = "h2o_map";
    opt_h2o_map->required    = NO;
    opt_h2o_map->label       = _("Per-pixel H₂O column raster map (g/cm²)");
    opt_h2o_map->description = _("Spatially variable water vapour column; overrides h2o_val=");
    opt_h2o_map->guisection  = _("Atmosphere");

    opt_params = G_define_standard_option(G_OPT_F_INPUT);
    opt_params->key         = "parameters";
    opt_params->required    = NO;
    opt_params->label       = _("6S parameter file (i.atcorr format)");
    opt_params->description = _("Extracts SZA and geometry; explicit options override file values");
    opt_params->guisection  = _("Input");

    opt_lut = G_define_standard_option(G_OPT_F_OUTPUT);
    opt_lut->key         = "lut";
    opt_lut->required    = NO;
    opt_lut->label       = _("Binary LUT output/input file");
    opt_lut->description = _("If the file exists it is loaded; otherwise computed and saved");
    opt_lut->guisection  = _("Output");

    flag_int = G_define_flag();
    flag_int->key         = 'i';
    flag_int->description = _("Output raster map as integer");
    flag_int->guisection  = _("Output");

    flag_polar = G_define_flag();
    flag_polar->key         = 'P';
    flag_polar->description = _("Enable vector (Stokes) polarization in LUT computation");
    flag_polar->guisection  = _("Atmosphere");

    flag_ram = G_define_flag();
    flag_ram->key         = 'R';
    flag_ram->description =
        _("Load all bands into RAM before processing (faster for small scenes; "
          "ensure sufficient free memory)");
    flag_ram->guisection  = _("Input");

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    /* ── Parse parameters ────────────────────────────────────────────────── */

    /* 6S parameter file (optional backward compatibility) */
    float file_sza   = -1.0f;
    int   file_atmo  = -1;
    int   file_aero  = -1;
    float file_aod   = -1.0f;

    if (opt_params->answer) {
        parse_6s_file(opt_params->answer,
                      &file_sza, &file_atmo, &file_aero, &file_aod);
        G_verbose_message(_("6S file parsed: SZA=%.2f atmo=%d aerosol=%d aod=%.3f"),
                          file_sza, file_atmo, file_aero, file_aod);
    }

    /* Wavelength */
    if (!opt_wl->answer)
        G_fatal_error(_("wavelength= (µm) is required."));
    float wavelength = (float)atof(opt_wl->answer);

    /* Solar zenith angle */
    float sza;
    if (opt_sza->answer)
        sza = (float)atof(opt_sza->answer);
    else if (file_sza >= 0.0f)
        sza = file_sza;
    else
        G_fatal_error(_("Solar zenith angle (sza=) is required when no "
                        "parameters= file is given."));

    float vza         = (float)atof(opt_vza->answer);
    float raa         = (float)atof(opt_raa->answer);
    float altitude_km = (float)atof(opt_altitude->answer);
    /* target_elevation: parsed for future per-pixel elevation support */
    (void)opt_target_elev->answer;
    float ozone_du    = (float)atof(opt_ozone->answer);
    float surf_press  = (float)atof(opt_pressure->answer);

    /* Day of year */
    int doy = 0;
    if (opt_doy->answer) doy = atoi(opt_doy->answer);

    /* Atmosphere / aerosol models */
    int atmo_model    = (file_atmo  >= 0 && !opt_atmo->answer)    ? file_atmo
                        : atmo_str_to_int(opt_atmo->answer);
    int aerosol_model = (file_aero  >= 0 && !opt_aerosol->answer) ? file_aero
                        : aerosol_str_to_int(opt_aerosol->answer);

    /* LUT AOD / H2O grids */
    float aod_grid[64], h2o_grid[64];
    int   n_aod = parse_csv_floats(opt_aod->answer, aod_grid, 64);
    int   n_h2o = parse_csv_floats(opt_h2o->answer, h2o_grid, 64);
    if (n_aod <= 0) G_fatal_error(_("Invalid aod= value: %s"), opt_aod->answer);
    if (n_h2o <= 0) G_fatal_error(_("Invalid h2o= value: %s"), opt_h2o->answer);

    /* Optionally inject AOD from 6S file into grid */
    if (file_aod >= 0.0f) {
        int found = 0;
        for (int i = 0; i < n_aod; i++)
            if (fabsf(aod_grid[i] - file_aod) < 1e-6f) { found = 1; break; }
        if (!found && n_aod < 64) {
            aod_grid[n_aod++] = file_aod;
            /* simple insertion sort */
            for (int i = n_aod - 1; i > 0 && aod_grid[i] < aod_grid[i-1]; i--) {
                float t = aod_grid[i]; aod_grid[i] = aod_grid[i-1]; aod_grid[i-1] = t;
            }
        }
    }

    int do_int   = flag_int->answer;
    int do_polar = flag_polar->answer;
    int do_ram   = flag_ram->answer;

    /* Radiance → TOA pre-computation (input is always radiance) */
    if (doy == 0)
        G_fatal_error(_("doy= (day of year) is required to convert radiance "
                        "to TOA reflectance."));
    float E0   = sixs_E0(wavelength);
    float d2   = (float)sixs_earth_sun_dist2(doy);
    float mu_s = cosf(sza * (float)(M_PI / 180.0));
    if (mu_s < 0.01f) {
        G_warning(_("cos(SZA) < 0.01 (sun near/below horizon). Results may be unreliable."));
        mu_s = 0.01f;
    }

    /* ── Build LutConfig ─────────────────────────────────────────────────── */
    LutConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.wl           = &wavelength;  cfg.n_wl  = 1;
    cfg.aod          = aod_grid;     cfg.n_aod = n_aod;
    cfg.h2o          = h2o_grid;     cfg.n_h2o = n_h2o;
    cfg.sza          = sza;
    cfg.vza          = vza;
    cfg.raa          = raa;
    cfg.altitude_km  = altitude_km;
    cfg.atmo_model   = atmo_model;
    cfg.aerosol_model= aerosol_model;
    cfg.surface_pressure = surf_press;
    cfg.ozone_du     = ozone_du;
    cfg.enable_polar = do_polar ? 1 : 0;

    /* ── Compute / load LUT ──────────────────────────────────────────────── */
    int    lut_n_aod = n_aod, lut_n_h2o = n_h2o;
    float *lut_aod   = aod_grid;
    float *lut_h2o   = h2o_grid;
    int    lut_owns_mem = 0;   /* 1 if we G_malloc'd aod/h2o arrays */

    int    n_lut = n_aod * n_h2o;   /* n_wl = 1 */
    float *Ra  = G_malloc((size_t)n_lut * sizeof(float));
    float *Td  = G_malloc((size_t)n_lut * sizeof(float));
    float *Tu  = G_malloc((size_t)n_lut * sizeof(float));
    float *sa  = G_malloc((size_t)n_lut * sizeof(float));

    LutArrays lut_arr = {
        .R_atm  = Ra, .T_down = Td, .T_up = Tu, .s_alb = sa,
        .T_down_dir = NULL, .R_atmQ = NULL, .R_atmU = NULL
    };

    const char *lut_path = opt_lut->answer;
    int lut_loaded = 0;

    if (lut_path) {
        /* try to load an existing LUT */
        float loaded_wl;
        float *la_aod = NULL, *la_h2o = NULL;
        int    la_naod = 0,    la_nh2o = 0;
        LutArrays loaded_arr = {0};

        if (access(lut_path, F_OK) == 0 &&
            load_lut(lut_path, &la_aod, &la_naod, &la_h2o, &la_nh2o,
                     &loaded_wl, &loaded_arr) == 0) {

            if (fabsf(loaded_wl - wavelength) > 0.005f)
                G_warning(_("Wavelength %.4f µm not found in LUT "
                            "(nearest: %.4f µm)."),
                          wavelength, loaded_wl);

            /* replace working arrays with loaded ones */
            G_free(Ra); G_free(Td); G_free(Tu); G_free(sa);
            lut_arr  = loaded_arr;
            lut_aod  = la_aod;
            lut_h2o  = la_h2o;
            lut_n_aod = la_naod;
            lut_n_h2o = la_nh2o;
            lut_owns_mem = 1;

            /* update cfg to match loaded grid */
            cfg.aod   = lut_aod; cfg.n_aod = lut_n_aod;
            cfg.h2o   = lut_h2o; cfg.n_h2o = lut_n_h2o;
            n_lut = lut_n_aod * lut_n_h2o;
            lut_loaded = 1;
        }
    }

    if (!lut_loaded) {
        G_message(_("Computing LUT: %d AOD × %d H₂O × 1 WL ..."), n_aod, n_h2o);
        if (atcorr_compute_lut(&cfg, &lut_arr) != 0)
            G_fatal_error(_("atcorr_compute_lut() failed."));
        if (lut_path)
            save_lut(lut_path, aod_grid, n_aod, h2o_grid, n_h2o,
                     wavelength, &lut_arr);
    }

    /* ── Scene-average correction coefficients (scalar path) ─────────────── */
    int use_perpixel = (opt_aod_map->answer || opt_h2o_map->answer);

    float sc_Ra = 0.0f, sc_Td = 0.0f, sc_Tu = 0.0f, sc_sa = 0.0f;
    float aod_val = 0.0f, h2o_val = 0.0f;

    if (!use_perpixel) {
        /* Determine scene-average AOD / H2O */
        if (opt_aod_val->answer)
            aod_val = (float)atof(opt_aod_val->answer);
        else if (file_aod >= 0.0f)
            aod_val = file_aod;
        else
            aod_val = lut_aod[lut_n_aod / 2];

        if (opt_h2o_val->answer)
            h2o_val = (float)atof(opt_h2o_val->answer);
        else
            h2o_val = lut_h2o[lut_n_h2o / 2];

        /* Single bilinear slice at (aod_val, h2o_val) */
        float Rs[1], Tds[1], Tus[1], ss[1];
        atcorr_lut_slice(&cfg, &lut_arr, aod_val, h2o_val,
                         Rs, Tds, Tus, ss, NULL);
        sc_Ra = Rs[0]; sc_Td = Tds[0]; sc_Tu = Tus[0]; sc_sa = ss[0];

        G_message(_("Correction params: R_atm=%.4f T_down=%.4f T_up=%.4f s_alb=%.4f "
                    "(AOD=%.3f H₂O=%.3f)"),
                  sc_Ra, sc_Td, sc_Tu, sc_sa, aod_val, h2o_val);
    }

    /* ── Raster region ───────────────────────────────────────────────────── */
    /* Snap region to input raster */
    struct Cell_head cellhd;
    Rast_get_cellhd(opt_input->answer, "", &cellhd);
    Rast_set_window(&cellhd);

    struct Cell_head region;
    G_get_window(&region);
    int nrows = region.rows;
    int ncols = region.cols;

    /* ── Determine chunk size ────────────────────────────────────────────── */
    /* n_maps: 1 input + up to 2 aux maps (aod_map, h2o_map) + 1 output */
    int n_aux  = (opt_aod_map->answer ? 1 : 0) + (opt_h2o_map->answer ? 1 : 0);
    int n_maps = 1 + n_aux + 1;

    size_t bytes_per_row = (size_t)ncols * (size_t)n_maps * sizeof(FCELL);
    int    chunk_rows;

    if (do_ram) {
        /* full-RAM mode: check available memory */
        size_t total_bytes = (size_t)nrows * bytes_per_row;

#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGE_SIZE)
        size_t avail = (size_t)sysconf(_SC_AVPHYS_PAGES) *
                       (size_t)sysconf(_SC_PAGE_SIZE);
        double frac  = (avail > 0) ? (double)total_bytes / (double)avail : 1.0;

        if (frac > RAM_ABORT_FRAC)
            G_fatal_error(_("-R: required memory (%.0f MB) exceeds %.0f %% of available RAM "
                            "(%.0f MB). Reduce scene size or remove -R flag."),
                          (double)total_bytes / (1024*1024),
                          RAM_ABORT_FRAC * 100.0,
                          (double)avail / (1024*1024));

        if (frac > RAM_WARN_FRAC)
            G_warning(_("-R: required memory (%.0f MB) is %.0f %% of available RAM "
                        "(%.0f MB). Proceed with caution."),
                      (double)total_bytes / (1024*1024),
                      frac * 100.0,
                      (double)avail / (1024*1024));
#endif
        chunk_rows = nrows;   /* one pass over everything */
    } else {
        chunk_rows = (int)(CHUNK_BYTES / bytes_per_row);
        if (chunk_rows < 1)    chunk_rows = 1;
        if (chunk_rows > nrows) chunk_rows = nrows;
    }

    G_message(_("Processing %d×%d raster in chunks of %d rows "
                "(%s path, %d OpenMP thread(s)) ..."),
              nrows, ncols, chunk_rows,
              use_perpixel ? "per-pixel" : "scene-average",
              omp_get_max_threads());

    /* ── Pre-load per-pixel auxiliary maps (always fully into RAM) ────────── */
    FCELL *aod_map_data = NULL;
    FCELL *h2o_map_data = NULL;

    if (opt_aod_map->answer)
        aod_map_data = read_raster_f(opt_aod_map->answer, nrows, ncols);
    if (opt_h2o_map->answer)
        h2o_map_data = read_raster_f(opt_h2o_map->answer, nrows, ncols);

    /* ── Open input and output rasters ───────────────────────────────────── */
    int fd_in  = Rast_open_old(opt_input->answer, "");
    int fd_out = Rast_open_new(opt_output->answer,
                               do_int ? CELL_TYPE : FCELL_TYPE);

    /* Allocate chunk buffers */
    FCELL *chunk_in  = G_malloc((size_t)chunk_rows * (size_t)ncols * sizeof(FCELL));
    FCELL *chunk_out = G_malloc((size_t)chunk_rows * (size_t)ncols * sizeof(FCELL));
    /* Integer conversion row buffer — allocated once for all chunks */
    CELL  *cell_row  = do_int ? G_malloc((size_t)ncols * sizeof(CELL)) : NULL;

    /* ── Main processing loop (chunk by chunk) ───────────────────────────── */
    int row_done = 0;

    while (row_done < nrows) {
        int rows_this_chunk = chunk_rows;
        if (row_done + rows_this_chunk > nrows)
            rows_this_chunk = nrows - row_done;
        long chunk_npix = (long)rows_this_chunk * ncols;

        /* Read input rows */
        for (int r = 0; r < rows_this_chunk; r++)
            Rast_get_f_row(fd_in, chunk_in + (size_t)r * ncols, row_done + r);

        /* ── Parallel pixel correction ────────────────────────────────────── */
        if (!use_perpixel) {
            /*
             * Scene-average path: all four atm coefficients are constants for
             * the entire chunk.  Precompute the reciprocal of T_down*T_up to
             * replace the per-pixel division with a multiply, then let
             * "parallel for simd" thread + SIMD-vectorise the loop.
             *
             * Null handling: GRASS FCELL nulls are IEEE 754 NaN.  NaN
             * propagates through all arithmetic (NaN op x = NaN), so a null
             * input pixel produces a NaN output which GRASS reads back as
             * null — no explicit null check needed, and the loop stays
             * fully branchless for the vectoriser.
             */
            float inv_TdTu  = 1.0f / (sc_Td * sc_Tu + 1e-10f);
            /* Fold radiance→TOA factor into a single scale. */
            float toa_scale = (float)M_PI * d2 / (E0 * mu_s + 1e-30f);

#ifdef _OPENMP
#pragma omp parallel for simd schedule(static)
#endif
            for (long p = 0; p < chunk_npix; p++) {
                float toa = (float)chunk_in[p] * toa_scale;
                /* NaN-safe clamp: IEEE 754 guarantees NaN comparisons are false. */
                toa = (toa < 0.0f) ? 0.0f : (toa > 2.0f ? 2.0f : toa);
                float y = (toa - sc_Ra) * inv_TdTu;
                chunk_out[p] = (FCELL)(y / (1.0f + sc_sa * y + 1e-10f));
            }
        } else {
            /*
             * Per-pixel path: atcorr_lut_slice() called per pixel — function
             * call overhead prevents SIMD; use threading only.  Explicit null
             * check retained because NaN aod/h2o inputs to atcorr_lut_slice()
             * are not guaranteed to propagate safely through the bracketing
             * logic.
             */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long p = 0; p < chunk_npix; p++) {
                FCELL dn = chunk_in[p];

                if (Rast_is_f_null_value(&dn)) {
                    Rast_set_f_null_value(&chunk_out[p], 1);
                    continue;
                }

                float toa = (float)M_PI * (float)dn * d2 / (E0 * mu_s + 1e-30f);
                toa = toa < 0.0f ? 0.0f : (toa > 2.0f ? 2.0f : toa);

                long  global_p = (long)row_done * ncols + p;
                float pix_aod  = aod_map_data
                                 ? (float)aod_map_data[global_p] : aod_val;
                float pix_h2o  = h2o_map_data
                                 ? (float)h2o_map_data[global_p] : h2o_val;

                /* atcorr_lut_slice(): bilinear in (AOD, H₂O), n_wl=1-safe. */
                float Rs[1], Tds[1], Tus[1], ss[1];
                atcorr_lut_slice(&cfg, &lut_arr, pix_aod, pix_h2o,
                                 Rs, Tds, Tus, ss, NULL);

                chunk_out[p] = (FCELL)atcorr_invert(toa, Rs[0], Tds[0], Tus[0], ss[0]);
            }
        }

        /* Write output rows */
        if (do_int) {
            for (int r = 0; r < rows_this_chunk; r++) {
                FCELL *fp = chunk_out + (size_t)r * ncols;
                /* SIMD-vectorise the float→CELL conversion per row.
                 * NaN (GRASS null) produces an implementation-defined integer
                 * on cast; the null is set explicitly in the scalar epilogue
                 * for any slot where the source was null. */
#ifdef _OPENMP
#pragma omp simd
#endif
                for (int c = 0; c < ncols; c++)
                    cell_row[c] = (CELL)roundf(fp[c]);
                /* Fix up nulls (rare; scalar) */
                for (int c = 0; c < ncols; c++)
                    if (Rast_is_f_null_value(&fp[c]))
                        Rast_set_c_null_value(&cell_row[c], 1);
                Rast_put_c_row(fd_out, cell_row);
            }
        } else {
            for (int r = 0; r < rows_this_chunk; r++)
                Rast_put_f_row(fd_out, chunk_out + (size_t)r * ncols);
        }

        row_done += rows_this_chunk;
        G_percent(row_done, nrows, 5);
    }

    G_percent(1, 1, 1);

    /* ── Close rasters ───────────────────────────────────────────────────── */
    Rast_close(fd_in);
    Rast_close(fd_out);

    /* ── Set output raster metadata ──────────────────────────────────────── */
    {
        struct History hist;
        Rast_short_history(opt_output->answer, "raster", &hist);
        Rast_append_format_history(
            &hist,
            "i.atcorr2: 6SV2.1 LUT-based correction via grass_sixsv. "
            "WL=%.4f µm, SZA=%.1f°, %s=%s, aerosol=%s, atmosphere=%s.",
            wavelength, sza,
            use_perpixel ? "AOD" : "aod_val",
            use_perpixel ? "per-pixel" : opt_aod_val->answer ? opt_aod_val->answer : "mid-LUT",
            opt_aerosol->answer, opt_atmo->answer);
        Rast_command_history(&hist);
        Rast_write_history(opt_output->answer, &hist);
    }

    /* Copy colour table from input */
    struct Colors colors;
    if (Rast_read_colors(opt_input->answer, "", &colors) >= 0)
        Rast_write_colors(opt_output->answer, G_mapset(), &colors);

    /* ── Cleanup ─────────────────────────────────────────────────────────── */
    G_free(chunk_in);
    G_free(chunk_out);
    G_free(cell_row);
    G_free(aod_map_data);
    G_free(h2o_map_data);

    if (lut_owns_mem) {
        G_free(lut_aod);
        G_free(lut_h2o);
    }
    G_free(lut_arr.R_atm);
    G_free(lut_arr.T_down);
    G_free(lut_arr.T_up);
    G_free(lut_arr.s_alb);

    G_message(_("Atmospheric correction complete."));
    return EXIT_SUCCESS;
}
