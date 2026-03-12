#!/usr/bin/env python3
"""
MODULE:    i.atcorr2
AUTHOR:    Yann Chemin
PURPOSE:   Atmospheric correction of a single-band raster using
           the libsixsv Python bindings (6SV2.1 LUT-based engine).
COPYRIGHT: (C) 2026 by the GRASS Development Team
           This program is free software under the GNU General Public
           License (>=v2). Read the file COPYING that comes with GRASS
           for details.

NOTES:
  This module delegates all radiative-transfer computations to
  libgrass_sixsv.so (the shared library from libsixsv) via the
  Python ctypes bindings in libsixsv/python/atcorr.py.

  Key differences from i.atcorr (the classic C++ 6S module):
    - Input must be radiance (W m⁻² sr⁻¹ µm⁻¹); converted to TOA reflectance
      internally using the Thuillier solar spectrum and Earth-Sun distance.
    - Parameters are explicit GRASS options, no 6S parameter text file.
    - LUT-based approach: a 3-D grid [AOD × H2O × wavelength] is computed
      once, then bilinearly interpolated per pixel for spatially varying
      atmospheric corrections.
    - LUT can be saved to disk and reused across multiple bands.
    - Per-pixel AOD / H2O raster maps are supported.
    - The 6S parameter file from i.atcorr can still be supplied via
      parameters= for backward compatibility (geometry is extracted from it).
"""

# %module
# % description: Atmospheric correction using libsixsv Python bindings (6SV2.1 LUT-based)
# % keyword: imagery
# % keyword: atmospheric correction
# % keyword: radiometric conversion
# % keyword: radiance
# % keyword: reflectance
# % keyword: satellite
# % keyword: 6SV
# %end

# %option G_OPT_R_INPUT
# % required: no
# % guisection: Input
# %end

# %option G_OPT_R_OUTPUT
# % required: no
# % guisection: Output
# %end

# %option
# % key: sensor
# % type: string
# % required: no
# % label: Satellite/sensor name (optional, auto-detected from band=)
# % description: If omitted, the sensor is inferred from the band name; only needed when the band name is ambiguous across sensors
# % options: amazonia1,aster,avnir,cbers4a_mux,eo1_ali,geoeye1,ikonos,landsat1_mss,landsat2_mss,landsat3_mss,landsat4_mss,landsat4_tm,landsat5_mss,landsat5_tm,landsat7_etm,landsat8,landsat9_oli2,modis_terra,planetscope_0c_0d,planetscope_0e,planetscope_0f_10,planetscope_dove_r,planetscope_superdove,pleiades1a,pleiades1b,pleiades_neo,prism_b,prism_f,prism_n,quickbird2,rapideye,sentinel2a,sentinel2b,sentinel2c,skysat_1_2,skysat_3_13,skysat_14_19,spot6,spot7,vgt1_spot4,vgt2_spot5,worldview2,worldview3,worldview4
# % guisection: Band
# %end

# %option
# % key: band
# % type: string
# % required: no
# % label: Band name within the selected sensor
# % description: Exact band column name from the sensor SRF table; run with -l flag (and sensor=) to list available names
# % guisection: Band
# %end

# %option
# % key: wavelength
# % type: double
# % required: no
# % label: Centre wavelength of the input band (µm)
# % description: Overrides the automatic lookup from band=; required when band= is not given and parameters= is absent
# % guisection: Optional
# %end

# %option
# % key: sza
# % type: double
# % required: no
# % label: Solar zenith angle (degrees)
# % description: Required unless parameters= (6S file) is given
# % guisection: Geometry
# %end

# %option
# % key: vza
# % type: double
# % required: no
# % answer: 5.0
# % description: View zenith angle (degrees)
# % guisection: Geometry
# %end

# %option
# % key: raa
# % type: double
# % required: no
# % answer: 90.0
# % description: Relative azimuth angle (degrees)
# % guisection: Geometry
# %end

# %option
# % key: altitude
# % type: double
# % required: no
# % answer: 1000.0
# % label: Sensor altitude (km)
# % description: Use >900 for satellite orbit (default 1000); 0 = ground level
# % guisection: Geometry
# %end

# %option
# % key: target_elevation
# % type: double
# % required: no
# % answer: 0.0
# % label: Mean target elevation above sea level (km)
# % description: Scene-average ground elevation in km; matches the 6S parameter file convention (e.g. -0.600 for 600 m asl); overridden by elevation= map
# % guisection: Geometry
# %end

# %option
# % key: doy
# % type: integer
# % required: no
# % label: Day of year [1-365]
# % description: Used to compute Earth-Sun distance for radiance to TOA reflectance conversion
# % guisection: Geometry
# %end

# %option
# % key: atmosphere
# % type: string
# % required: no
# % answer: us62
# % label: Atmospheric model
# % options: none,tropical,midsum,midwin,subarctsum,suarctwint,us62
# % descriptions: none;No gaseous absorption;tropical;Tropical;midsum;Mid-latitude summer;midwin;Mid-latitude winter;subarctsum;Sub-arctic summer;suarctwint;Sub-arctic winter;us62;US Standard 1962 (default)
# % guisection: Atmosphere
# %end

# %option
# % key: aerosol
# % type: string
# % required: no
# % answer: continental
# % label: Aerosol model
# % options: none,continental,maritime,urban,desert
# % descriptions: none;No aerosols;continental;Continental (default);maritime;Maritime;urban;Urban;desert;Desert dust
# % guisection: Atmosphere
# %end

# %option
# % key: aod
# % type: string
# % required: no
# % answer: 0.0,0.1,0.2,0.5
# % label: AOD at 550 nm LUT grid values (comma-separated)
# % description: Aerosol optical depth values for the correction LUT
# % guisection: Atmosphere
# %end

# %option
# % key: h2o
# % type: string
# % required: no
# % answer: 1.0,2.0,3.5
# % label: Column water vapour LUT grid values, g/cm² (comma-separated)
# % guisection: Atmosphere
# %end

# %option
# % key: aod_val
# % type: double
# % required: no
# % label: AOD value for correction (overrides aod_map=)
# % description: Single AOD at 550 nm for the scene; if omitted, mid-LUT value is used
# % guisection: Atmosphere
# %end

# %option
# % key: h2o_val
# % type: double
# % required: no
# % label: H2O column value for correction, g/cm² (overrides h2o_map=)
# % description: If omitted, mid-LUT value is used
# % guisection: Atmosphere
# %end

# %option
# % key: ozone
# % type: double
# % required: no
# % answer: 300.0
# % description: Ozone column (Dobson units, 0 = use standard atmosphere)
# % guisection: Atmosphere
# %end

# %option
# % key: surface_pressure
# % type: double
# % required: no
# % answer: 0.0
# % description: Surface pressure (hPa, 0 = use standard atmosphere)
# % guisection: Atmosphere
# %end

# %option G_OPT_R_ELEV
# % key: elevation
# % required: no
# % label: Input elevation raster map (m above sea level)
# % description: Per-pixel target elevation in metres; overrides target_elevation=
# % guisection: Input
# %end

# %option G_OPT_R_INPUT
# % key: visibility
# % required: no
# % label: Input visibility raster map (km)
# % description: Per-pixel meteorological visibility; overrides visibility_val=
# % guisection: Atmosphere
# %end

# %option
# % key: visibility_val
# % type: double
# % required: no
# % label: Scene visibility (km)
# % description: Meteorological visibility used to estimate AOD (Koschmieder approximation); ignored when aod_val= or aod_map= is given
# % guisection: Atmosphere
# %end

# %option G_OPT_R_INPUT
# % key: aod_map
# % required: no
# % label: Per-pixel AOD raster map
# % description: Spatially variable AOD at 550 nm; overrides aod_val=
# % guisection: Atmosphere
# %end

# %option G_OPT_R_INPUT
# % key: h2o_map
# % required: no
# % label: Per-pixel H2O column raster map (g/cm²)
# % description: Spatially variable water vapour column; overrides h2o_val=
# % guisection: Atmosphere
# %end

# %option G_OPT_F_INPUT
# % key: parameters
# % required: no
# % label: 6S parameter file (i.atcorr format)
# % description: When provided, extracts SZA and geometry from the file; explicit options override file values
# % guisection: Optional
# %end

# %option G_OPT_F_OUTPUT
# % key: lut
# % required: no
# % label: Binary LUT output/input file
# % description: If the file exists it is loaded; otherwise the LUT is computed and saved
# % guisection: Output
# %end

# %flag
# % key: l
# % description: List available band names for the selected sensor and exit
# % guisection: Band
# %end

# %flag
# % key: P
# % description: Enable vector (Stokes) polarization in LUT computation
# % guisection: Atmosphere
# %end

import os
import sys
import struct
import math

import grass.script as gs
import numpy as np

import sensors as _sensors

# ── grass_sixsv Python API ────────────────────────────────────────────────────

def _find_atcorr_api():
    """Locate atcorr.py installed by the grass_sixsv library (libsixsv).

    Search order:
      1. $GISBASE/scripts/          — system or g.extension install
      2. $GRASS_ADDON_BASE/scripts/ — per-user g.extension install
      3. ../libsixsv/python/        — sibling build-tree checkout
      4. $LIBSIXSV_PYTHON           — explicit env-var override
    """
    _gisbase   = os.environ.get("GISBASE", "")
    _addonbase = os.environ.get("GRASS_ADDON_BASE", "")
    _here      = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        # 1. Installed GRASS location (preferred at runtime)
        os.path.join(_gisbase, "scripts"),
        # 2. Per-user g.extension location
        os.path.join(_addonbase, "scripts"),
        # 3. Sibling build-tree (developer checkout: dev/libsixsv/python/)
        os.path.join(_here, "..", "libsixsv", "python"),
        # 4. Explicit env-var override
        os.environ.get("LIBSIXSV_PYTHON", ""),
    ]
    for d in candidates:
        if d and os.path.isfile(os.path.join(d, "atcorr.py")):
            return d
    return None


# ── 6S parameter file parser (i.atcorr backward compatibility) ────────────────

_ATMO_NAME_TO_INT = {
    "none": 0, "us62": 1, "midsum": 2, "midwin": 3,
    "tropical": 4, "subarctsum": 5, "suarctwint": 6,
}
_AEROSOL_NAME_TO_INT = {
    "none": 0, "continental": 1, "maritime": 2, "urban": 3, "desert": 5,
}

# Koschmieder approximation: AOD_550 ≈ 3.912/V - 0.01162  (V in km)
# Clamped to [0.001, 5.0]; matches 6S internal visibility→extinction logic.
def _vis_km_to_aod(vis_km):
    """Convert meteorological visibility (km) to approximate AOD at 550 nm.

    Uses the Koschmieder law: σ_ext = 3.912 / V, with a small correction
    term for background molecular scattering.

    Parameters
    ----------
    vis_km : float or numpy.ndarray
        Visibility in kilometres (scalar or array).

    Returns
    -------
    float or numpy.ndarray
        Estimated AOD at 550 nm (same shape as input), clipped to [0.001, 5].
    """
    aod = 3.912 / np.maximum(vis_km, 0.1) - 0.01162
    return np.clip(aod, 0.001, 5.0)


def _parse_6s_param_file(path):
    """Parse the 6S conditions file used by i.atcorr.

    Returns a dict with keys: sza, month, day, lon, lat, atmo_model,
    aerosol_model, visibility, target_elevation_km, sensor_altitude_km,
    iwave, aod (if given explicitly).

    Only the geometry / atmospheric-state parameters are extracted;
    the spectral band code (iwave) is read but wavelength lookup is
    left to the caller.
    """
    result = {}
    try:
        with open(path) as fh:
            lines = [l.split("#")[0].split("(")[0].strip()
                     for l in fh if l.strip() and not l.strip().startswith("#")]
    except OSError as exc:
        gs.fatal(f"Cannot read 6S parameter file: {exc}")

    idx = 0
    def _next():
        nonlocal idx
        while idx < len(lines) and not lines[idx]:
            idx += 1
        if idx >= len(lines):
            return None
        val = lines[idx]
        idx += 1
        return val

    # Line 1: geometrical condition code
    geom_code = int(_next() or 0)
    result["geom_code"] = geom_code

    # Line 2: month, day, hh.ddd, lon, lat  (for most sensor types)
    geom_line = (_next() or "").split()
    if len(geom_line) >= 5:
        result["month"] = int(geom_line[0])
        result["day"]   = int(geom_line[1])
        result["hour"]  = float(geom_line[2])
        result["lon"]   = float(geom_line[3])
        result["lat"]   = float(geom_line[4])
        # Approximate SZA from latitude + DOY (simple solar geometry)
        # Only used as a fallback if the user does not supply sza= explicitly.
        doy = _doy_from_month_day(result["month"], result["day"])
        result["doy"] = doy
        result["sza"] = _approx_sza(result["lat"], doy, result["hour"])
    elif len(geom_line) >= 3:
        result["month"] = int(geom_line[0])
        result["day"]   = int(geom_line[1])
        result["hour"]  = float(geom_line[2])
        doy = _doy_from_month_day(result["month"], result["day"])
        result["doy"] = doy

    # Line 3: atmospheric model
    result["atmo_model"] = int(_next() or 6)

    # Line 4: aerosol model
    result["aerosol_model"] = int(_next() or 1)

    # Line 5: visibility km  (or 0 if AOD given on next line)
    vis_str = _next() or "15"
    vis = float(vis_str)
    result["visibility_km"] = vis
    if vis == 0.0:
        # Next line is AOD at 550 nm
        aod_str = _next() or "0.2"
        result["aod"] = float(aod_str)
    elif vis < 0:
        # vis < 0 means no aerosols
        result["aod"] = 0.0
    else:
        # Convert visibility to rough AOD (Koschmieder approximation)
        result["aod"] = 3.912 / vis if vis > 0 else 0.2

    # Line 6: target elevation (negative km)
    xps_str = _next() or "0"
    xps = float(xps_str)
    result["target_elevation_km"] = abs(xps) if xps <= 0 else 0.0

    # Line 7: sensor altitude (-1000 = satellite)
    xpp_str = _next() or "-1000"
    xpp = float(xpp_str)
    result["sensor_altitude_km"] = 1000.0 if xpp == -1000 else abs(xpp)

    # Line 8: iwave (spectral band code)
    iwave_str = _next() or "-1"
    result["iwave"] = int(iwave_str)

    return result


def _doy_from_month_day(month, day):
    """Convert a calendar date to day-of-year (no leap-year correction).

    Parameters
    ----------
    month : int
        Month (1–12).
    day : int
        Day of month (1–31).

    Returns
    -------
    int
        Day of year (1–365), assuming a non-leap year.
    """
    days_in_month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return sum(days_in_month[:month]) + day


def _approx_sza(lat_deg, doy, hour_ut):
    """Approximate solar zenith angle (degrees) – accuracy ~1-2°."""
    lat  = math.radians(lat_deg)
    decl = math.radians(23.45 * math.sin(math.radians(360.0 / 365.0 * (doy - 81))))
    ha   = math.radians(15.0 * (hour_ut - 12.0))
    cos_sza = (math.sin(lat) * math.sin(decl) +
               math.cos(lat) * math.cos(decl) * math.cos(ha))
    cos_sza = max(-1.0, min(1.0, cos_sza))
    return math.degrees(math.acos(cos_sza))


# ── LUT file I/O (compatible with i.hyper.atcorr LUT format) ─────────────────

LUT_MAGIC   = 0x4C555400
LUT_VERSION = 1


def _save_lut(path, cfg, arrays):
    """Write an i.hyper.atcorr-compatible binary LUT file."""
    with open(path, "wb") as fh:
        fh.write(struct.pack("<II", LUT_MAGIC, LUT_VERSION))
        fh.write(struct.pack("<iii", len(cfg.aod), len(cfg.h2o), len(cfg.wl)))
        cfg.aod.astype("<f4").tofile(fh)
        cfg.h2o.astype("<f4").tofile(fh)
        cfg.wl.astype("<f4").tofile(fh)
        arrays.R_atm.ravel().astype("<f4").tofile(fh)
        arrays.T_down.ravel().astype("<f4").tofile(fh)
        arrays.T_up.ravel().astype("<f4").tofile(fh)
        arrays.s_alb.ravel().astype("<f4").tofile(fh)
    gs.verbose(f"LUT saved to {path}")


def _load_lut(path, atcorr_mod):
    """Load a binary LUT written by i.hyper.atcorr or _save_lut().

    Returns (cfg_like, LutArrays) where cfg_like carries .aod, .h2o, .wl.
    """
    LutConfig  = atcorr_mod.LutConfig
    LutArrays  = atcorr_mod.LutArrays
    with open(path, "rb") as fh:
        magic, version = struct.unpack("<II", fh.read(8))
        if magic != LUT_MAGIC:
            gs.fatal(f"Not a valid LUT file: {path}")
        n_aod, n_h2o, n_wl = struct.unpack("<iii", fh.read(12))
        aod  = np.frombuffer(fh.read(n_aod * 4), dtype="<f4").copy()
        h2o  = np.frombuffer(fh.read(n_h2o * 4), dtype="<f4").copy()
        wl   = np.frombuffer(fh.read(n_wl  * 4), dtype="<f4").copy()
        n    = n_aod * n_h2o * n_wl
        R_atm  = np.frombuffer(fh.read(n * 4), dtype="<f4").copy().reshape(n_aod, n_h2o, n_wl)
        T_down = np.frombuffer(fh.read(n * 4), dtype="<f4").copy().reshape(n_aod, n_h2o, n_wl)
        T_up   = np.frombuffer(fh.read(n * 4), dtype="<f4").copy().reshape(n_aod, n_h2o, n_wl)
        s_alb  = np.frombuffer(fh.read(n * 4), dtype="<f4").copy().reshape(n_aod, n_h2o, n_wl)
    gs.verbose(f"LUT loaded from {path}: AOD={list(aod)}, H2O={list(h2o)}, WL={list(wl)}")
    # Reconstruct a minimal cfg-like object
    cfg = LutConfig(wl=wl, aod=aod, h2o=h2o)
    lut = LutArrays(R_atm=R_atm, T_down=T_down, T_up=T_up, s_alb=s_alb)
    return cfg, lut


# ── Per-pixel correction helpers ──────────────────────────────────────────────

def _bilinear_interp_lut(arrays, cfg, aod_arr, h2o_arr, wl_idx):
    """Vectorised bilinear interpolation into the LUT for 2-D aod/h2o arrays.

    Returns R_atm, T_down, T_up, s_alb arrays of shape (nrows, ncols).
    """
    aod_grid = cfg.aod   # shape (n_aod,)
    h2o_grid = cfg.h2o   # shape (n_h2o,)

    # Clamp to grid bounds
    aod_arr = np.clip(aod_arr, aod_grid[0], aod_grid[-1])
    h2o_arr = np.clip(h2o_arr, h2o_grid[0], h2o_grid[-1])

    # Find bracketing indices
    ia = np.searchsorted(aod_grid, aod_arr, side="right").clip(1, len(aod_grid) - 1)
    ih = np.searchsorted(h2o_grid, h2o_arr, side="right").clip(1, len(h2o_grid) - 1)

    da = (aod_arr - aod_grid[ia - 1]) / (aod_grid[ia] - aod_grid[ia - 1] + 1e-30)
    dh = (h2o_arr - h2o_grid[ih - 1]) / (h2o_grid[ih] - h2o_grid[ih - 1] + 1e-30)

    def _interp(arr4d):
        # arr4d shape: (n_aod, n_h2o, n_wl) → extract wl_idx → (n_aod, n_h2o)
        v = arr4d[:, :, wl_idx]
        v00 = v[ia - 1, ih - 1]
        v10 = v[ia,     ih - 1]
        v01 = v[ia - 1, ih    ]
        v11 = v[ia,     ih    ]
        return (v00 * (1 - da) * (1 - dh) +
                v10 * da       * (1 - dh) +
                v01 * (1 - da) * dh       +
                v11 * da       * dh)

    return (_interp(arrays.R_atm),
            _interp(arrays.T_down),
            _interp(arrays.T_up),
            _interp(arrays.s_alb))


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    """GRASS module entry point for single-band atmospheric correction.

    Parses GRASS options/flags, locates the i.hyper.atcorr Python API
    (``atcorr.py``), builds a 6SV2.1 look-up table (LUT) over the requested
    AOD × H₂O grid, applies Lambertian surface reflectance inversion to the
    input raster, and writes the atmospherically corrected output raster.

    The function handles two correction modes:

    - **Scene-average**: A single AOD/H₂O pair is used for the whole scene.
      The LUT is sliced once and the inversion is applied vectorially across
      all pixels.
    - **Per-pixel**: Spatially varying AOD and/or H₂O rasters are read.
      LUT values are bilinearly interpolated at each pixel location before
      inversion.

    Radiance inputs are converted to TOA reflectance using:

    .. math::

        \\rho_\\text{toa} = \\frac{\\pi \\, L \\, d^2}{E_0 \\, \\cos\\theta_s}

    where *d²* is the Earth-Sun distance correction factor, *E₀* the
    Thuillier 2003 solar irradiance at the requested wavelength, and
    *θ_s* the solar zenith angle.

    Notes
    -----
    All GRASS module options and flags are consumed via ``gs.parser()``;
    see the ``#%`` header block at the top of this file for the full option
    specification.
    """
    opts, flags = gs.parser()

    # ── -l: list bands for sensor and exit ───────────────────────────────────
    if flags["l"]:
        sensor_key = opts["sensor"]
        if not sensor_key:
            # No sensor given — list all sensors
            print("Available sensors (use sensor= to filter):")
            for key in sorted(_sensors.SENSORS.keys()):
                print(f"  {key}")
            return
        try:
            specs = _sensors.load_band_specs(sensor_key)
        except (KeyError, FileNotFoundError) as exc:
            gs.fatal(f"Cannot load sensor '{sensor_key}': {exc}")
        print(f"Bands for sensor '{sensor_key}':")
        print(f"  {'Band name':<25}  {'Centre (µm)':>12}  {'FWHM (nm)':>10}")
        print(f"  {'-'*25}  {'-'*12}  {'-'*10}")
        for band, (centre_um, fwhm_nm) in specs.items():
            print(f"  {band:<25}  {centre_um:>12.4f}  {fwhm_nm:>10.1f}")
        return

    # ── Validate required options (not enforced by parser when -l is available)
    missing = []
    if not opts["input"]:
        missing.append("input")
    if not opts["output"]:
        missing.append("output")
    if not opts["doy"]:
        missing.append("doy")
    if missing:
        gs.fatal("Required parameter(s) not set: " + ", ".join(missing))

    # ── Locate and import the libsixsv Python API ────────────────────────────
    api_dir = _find_atcorr_api()
    if api_dir is None:
        gs.fatal(
            "Cannot locate atcorr.py from libsixsv.\n"
            "Set the LIBSIXSV_PYTHON environment variable to the "
            "directory containing atcorr.py, or place libsixsv/ "
            "next to i.atcorr2/ in the same parent directory."
        )
    sys.path.insert(0, api_dir)
    try:
        import atcorr as _atcorr
    except ImportError as exc:
        gs.fatal(f"Failed to import atcorr API from {api_dir}: {exc}")

    gs.verbose(f"Using atcorr API from {api_dir}  (libgrass_sixsv version: {_atcorr.version()})")

    # ── Parse 6S parameter file if supplied ──────────────────────────────────
    file_params = {}
    if opts["parameters"]:
        file_params = _parse_6s_param_file(opts["parameters"])
        gs.verbose(f"6S parameter file parsed: {file_params}")

    # ── Geometry parameters (options override file values) ────────────────────
    sza = float(opts["sza"]) if opts["sza"] else file_params.get("sza")
    if sza is None:
        gs.fatal("Solar zenith angle (sza=) is required when no parameters= file is given.")

    vza          = float(opts["vza"])
    raa          = float(opts["raa"])
    altitude_km  = float(opts["altitude"])
    # target_elevation= is in km (6S convention); elevation= raster is in
    # metres (G_OPT_R_ELEV standard) — scene-average derived below.
    target_elev  = float(opts["target_elevation"])

    doy = int(opts["doy"])

    # ── Atmospheric state ─────────────────────────────────────────────────────
    atmo_model    = _ATMO_NAME_TO_INT[opts["atmosphere"]]
    aerosol_model = _AEROSOL_NAME_TO_INT[opts["aerosol"]]

    # Override with file values if present and option is still at default
    if not opts["atmosphere"] and "atmo_model" in file_params:
        atmo_model = file_params["atmo_model"]
    if not opts["aerosol"] and "aerosol_model" in file_params:
        aerosol_model = file_params["aerosol_model"]

    ozone_du         = float(opts["ozone"])
    surface_pressure = float(opts["surface_pressure"])

    # ── LUT grid ──────────────────────────────────────────────────────────────
    aod_grid = np.array([float(v) for v in opts["aod"].split(",")], dtype=np.float32)
    h2o_grid = np.array([float(v) for v in opts["h2o"].split(",")], dtype=np.float32)

    # If the 6S file gave us a single AOD estimate, include it in the grid
    if "aod" in file_params:
        fa = file_params["aod"]
        if not any(abs(fa - v) < 1e-6 for v in aod_grid):
            aod_grid = np.sort(np.append(aod_grid, fa)).astype(np.float32)

    # ── Wavelength ────────────────────────────────────────────────────────────
    wavelength = float(opts["wavelength"]) if opts["wavelength"] else None

    if wavelength is None and opts["band"]:
        band_key = opts["band"]
        sensor_key = opts["sensor"] or None

        # Auto-detect sensor from band name when sensor= is not given
        if sensor_key is None:
            try:
                sensor_key = _sensors.find_sensor_for_band(band_key)
                gs.verbose(f"Auto-detected sensor '{sensor_key}' for band '{band_key}'")
            except ValueError as exc:
                gs.fatal(str(exc))

        try:
            band_centers = _sensors.load_band_centers(sensor_key)
        except (KeyError, FileNotFoundError) as exc:
            gs.fatal(f"Cannot load sensor '{sensor_key}': {exc}")

        if band_key not in band_centers:
            available = ", ".join(sorted(band_centers.keys()))
            gs.fatal(
                f"Band '{band_key}' not found for sensor '{sensor_key}'. "
                f"Available bands: {available}"
            )
        wavelength = band_centers[band_key]
        gs.verbose(
            f"Sensor '{sensor_key}' band '{band_key}' → "
            f"centre wavelength {wavelength:.4f} µm"
        )

    if wavelength is None:
        gs.fatal(
            "wavelength= (µm) is required, or specify band= "
            "(optionally sensor=) for automatic SRF-weighted lookup."
        )
    wl_arr = np.array([wavelength], dtype=np.float32)

    flag_polar  = flags["P"]   # polarization

    # ── Build LutConfig ───────────────────────────────────────────────────────
    cfg = _atcorr.LutConfig(
        wl            = wl_arr,
        aod           = aod_grid,
        h2o           = h2o_grid,
        sza           = sza,
        vza           = vza,
        raa           = raa,
        altitude_km   = altitude_km,
        atmo_model    = atmo_model,
        aerosol_model = aerosol_model,
        surface_pressure = surface_pressure,
        ozone_du      = ozone_du,
        enable_polar  = 1 if flag_polar else 0,
    )

    # ── Compute / load LUT ────────────────────────────────────────────────────
    lut_path = opts["lut"] if opts["lut"] else None

    if lut_path and os.path.exists(lut_path):
        gs.message(f"Loading pre-computed LUT from {lut_path} ...")
        cfg_lut, arrays = _load_lut(lut_path, _atcorr)
        # Verify wavelength is present in the loaded LUT
        wl_diffs = np.abs(cfg_lut.wl - wavelength)
        wl_idx = int(np.argmin(wl_diffs))
        if wl_diffs[wl_idx] > 0.005:
            gs.warning(
                f"Wavelength {wavelength} µm not found in LUT "
                f"(nearest: {cfg_lut.wl[wl_idx]:.4f} µm). "
                "Recompute the LUT with the correct wavelength."
            )
        cfg = cfg_lut
    else:
        gs.message(
            f"Computing LUT: {len(aod_grid)} AOD × {len(h2o_grid)} H2O × 1 WL ..."
        )
        arrays = _atcorr.compute_lut(cfg)
        wl_idx = 0
        if lut_path:
            _save_lut(lut_path, cfg, arrays)

    # ── Read input raster ─────────────────────────────────────────────────────
    input_map  = opts["input"]
    output_map = opts["output"]

    # Adjust region to input raster for the duration of this run, then restore.
    # gs.use_temp_region() creates a process-local WIND_OVERRIDE so the user's
    # current region is unchanged after the module exits (mirrors i.atcorr
    # behaviour: G_get_set_window / Rast_set_window).
    with gs.use_temp_region():
        gs.run_command("g.region", raster=input_map, quiet=True)
        region = gs.region()
        nrows  = region["rows"]
        ncols  = region["cols"]

        gs.message(f"Reading {input_map} ({nrows}×{ncols}) ...")
        data = gs.array.array()
        data.read(input_map)                     # shape (nrows, ncols), float64
        L = data.array.astype(np.float64)        # radiance in W m⁻² sr⁻¹ µm⁻¹

        null_mask = np.isnan(L)

        # Override target_elev with scene-average of elevation= map (metres → km)
        if opts["elevation"]:
            elev_data = gs.array.array()
            elev_data.read(opts["elevation"])
            elev_m = elev_data.array
            valid = elev_m[~np.isnan(elev_m)]
            if valid.size > 0:
                target_elev = float(np.mean(valid)) / 1000.0
                gs.verbose(f"Mean elevation {target_elev * 1000:.0f} m → target_elevation {target_elev:.3f} km")

        # ── Radiance → TOA reflectance conversion ─────────────────────────────
        # ρ_toa = (π × L × d²) / (E₀ × cos θs)
        E0   = _atcorr.solar_E0(wavelength)   # W m⁻² µm⁻¹
        d2   = _atcorr.earth_sun_dist2(doy)   # AU²
        mu_s = math.cos(math.radians(sza))
        if mu_s < 0.01:
            gs.warning("cos(SZA) < 0.01 (sun near/below horizon). Results may be unreliable.")
            mu_s = 0.01
        rho_toa = (math.pi * L * float(d2)) / (E0 * mu_s)
        rho_toa = np.clip(rho_toa, 0.0, 2.0)

        # ── Build per-pixel AOD / H2O fields ──────────────────────────────────
        use_perpixel = False

        if opts["aod_map"] or opts["visibility"] or opts["h2o_map"]:
            use_perpixel = True

            if opts["aod_map"]:
                aod_data = gs.array.array()
                aod_data.read(opts["aod_map"])
                aod_scene = aod_data.array.astype(np.float32)
            elif opts["visibility"]:
                # Convert per-pixel visibility (km) to AOD via Koschmieder
                vis_data = gs.array.array()
                vis_data.read(opts["visibility"])
                aod_scene = _vis_km_to_aod(vis_data.array.astype(np.float32))
            else:
                aod_val = (float(opts["aod_val"]) if opts["aod_val"]
                           else float(np.median(aod_grid)))
                aod_scene = np.full((nrows, ncols), aod_val, dtype=np.float32)

            if opts["h2o_map"]:
                h2o_data = gs.array.array()
                h2o_data.read(opts["h2o_map"])
                h2o_scene = h2o_data.array.astype(np.float32)
            else:
                h2o_val = (float(opts["h2o_val"]) if opts["h2o_val"]
                           else float(np.median(h2o_grid)))
                h2o_scene = np.full((nrows, ncols), h2o_val, dtype=np.float32)

        if not use_perpixel:
            # Scene-average AOD: aod_val= > visibility_val= > 6S file > mid-LUT
            if opts["aod_val"]:
                aod_val = float(opts["aod_val"])
            elif opts["visibility_val"]:
                aod_val = float(_vis_km_to_aod(float(opts["visibility_val"])))
                gs.verbose(f"Converted visibility {opts['visibility_val']} km → AOD {aod_val:.3f}")
            else:
                aod_val = float(file_params.get("aod", np.median(aod_grid)))
            h2o_val = (float(opts["h2o_val"]) if opts["h2o_val"]
                       else float(np.median(h2o_grid)))

        # ── Atmospheric correction ─────────────────────────────────────────────
        gs.message("Applying atmospheric correction ...")

        if use_perpixel:
            # Per-pixel correction via vectorised bilinear interpolation
            R_atm, T_down, T_up, s_alb = _bilinear_interp_lut(
                arrays, cfg, aod_scene, h2o_scene, wl_idx
            )
            refl_boa = _atcorr.invert(rho_toa, R_atm, T_down, T_up, s_alb)
        else:
            # Scene-average correction
            sl = _atcorr.lut_slice(cfg, arrays,
                                   aod_val=aod_val, h2o_val=h2o_val)
            refl_boa = _atcorr.invert(rho_toa,
                                      float(sl.R_atm[0]),
                                      float(sl.T_down[0]),
                                      float(sl.T_up[0]),
                                      float(sl.s_alb[0]))

        out = refl_boa
        out[null_mask] = np.nan

        # ── Write output raster ────────────────────────────────────────────────
        gs.message(f"Writing {output_map} ...")
        out_arr = gs.array.array(np.float64)
        out_arr.array = out.astype(np.float64)

        out_arr.write(output_map, overwrite=gs.overwrite())

    # Copy colour table from input
    gs.run_command("r.colors", map=output_map, raster=input_map, quiet=True)

    # Write metadata
    gs.run_command(
        "r.support", map=output_map,
        title=f"Atmospherically corrected {input_map}",
        description=(
            f"i.atcorr2: 6SV2.1 LUT-based correction via i.hyper.atcorr. "
            f"WL={wavelength} µm, SZA={sza:.1f}°, AOD={aod_val if not use_perpixel else 'per-pixel'}, "
            f"H2O={h2o_val if not use_perpixel else 'per-pixel'}, "
            f"aerosol={opts['aerosol']}, atmosphere={opts['atmosphere']}."
        ),
        quiet=True,
    )
    gs.run_command("r.timestamp", map=output_map, date="none", quiet=True)

    gs.message("Atmospheric correction complete.")


if __name__ == "__main__":
    main()
