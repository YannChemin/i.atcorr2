## DESCRIPTION

**i.atcorr2** performs atmospheric correction of a single-band raster
using the 6SV2.1 radiative-transfer engine exposed by the
[`grass_sixsv`](https://github.com/YannChemin/libsixsv) Python bindings.

It is a Python counterpart to [i.atcorr](../grass/imagery/i.atcorr/) that
replaces the internal C++ 6S pixel-by-pixel computation with a pre-computed
**look-up table (LUT)** approach.  The LUT is a 3-D array indexed by
[AOD × H₂O × wavelength]; once computed it can be saved to disk and reused
across multiple bands or time steps with no additional radiative-transfer cost.

> **Important:** current region settings are ignored.  The region is adjusted
> to cover the input raster map before correction is applied.  The previous
> settings are restored afterwards.

All atmospheric parameters (`sza`, `atmosphere`, `aerosol`, etc.) are
supplied directly as GRASS options — no separate 6S conditions text file is
required.  For backward compatibility the `parameters=` option accepts the
same conditions file used by *i.atcorr*; explicit options always take
precedence.

Spectral band centre wavelengths can be determined automatically from the
`sensor=` and `band=` options, which use SRF-weighted mean wavelengths
computed from the sensor CSV files in `sensors_csv/`.  Explicit `wavelength=`
overrides the automatic lookup.

---

### i.atcorr vs i.atcorr2 — feature comparison

| Feature | *i.atcorr* (classic) | *i.atcorr2* (this module) |
|---|---|---|
| **Language** | C++ | Python |
| **Engine** | 6S C++ port, pixel-by-pixel | 6SV2.1 via `grass_sixsv`, LUT-based |
| **Parameter input** | 6S conditions text file (`parameters=`, required) | Explicit GRASS options; 6S file optionally accepted for backward compatibility |
| **Spectral band** | Sensor code (iwave −2 to 33) or wl range in conditions file | `sensor=`/`band=` SRF lookup (24 sensors) **or** `wavelength=` µm |
| **Atmospheric model** | Numeric code 0–6 in conditions file | Named option: `us62`, `tropical`, `midsum`, … |
| **Aerosol concentration** | Visibility (km) or AOD in conditions file | `aod=` LUT grid + `aod_val=` scene value |
| **Water vapour** | Fixed by atmospheric model | `h2o=` LUT grid + `h2o_val=` scene value |
| **Per-pixel altitude** | `elevation=` raster (metres) | `elevation=` raster (metres); scene mean converted to km for 6S |
| **Per-pixel visibility** | `visibility=` raster (km) | `visibility=` raster (km) or `visibility_val=` scalar; converted to AOD via Koschmieder approximation |
| **Per-pixel AOD** | — | `aod_map=` raster |
| **Per-pixel H₂O** | — | `h2o_map=` raster |
| **LUT caching** | RB-tree keyed on (alt, vis) bins | Full 3-D [AOD × H₂O × WL] grid, saveable to `lut=` file |
| **Multi-band reuse** | Each band reruns the full 6S computation | Compute LUT once, reuse across bands via `lut=` |
| **Polarization** | No | `-P` flag (Stokes I, Q, U via 6SV2.1 vector RT) |
| **Input** | Radiance or reflectance (`-r` flag) | Radiance (W m⁻² sr⁻¹ µm⁻¹) only; `doy=` required |
| **Radiance → TOA** | Handled internally by 6S transform | Explicit: `ρ_toa = (π·L·d²) / (E₀·cos θs)` using Thuillier spectrum |
| **Output** | Float or integer (`-i` flag); rescalable | Float reflectance [0, 1] |
| **Dependencies** | GRASS GIS (C build) | GRASS GIS + `grass_sixsv` + NumPy |

---

### Processing pipeline

1. Parse GRASS options (and optionally a 6S conditions file for geometry).
2. If `sensor=` and `band=` are given, look up the SRF-weighted centre
   wavelength from `sensors_csv/`; `wavelength=` overrides this.
3. Build a `LutConfig`: geometry + atmospheric state + AOD/H₂O grid + centre
   wavelength.
4. If `lut=` file already exists: load it.  Otherwise call
   `atcorr.compute_lut()` via `grass_sixsv` and, if `lut=` is specified,
   save the result for future reuse.
5. Adjust the computational region to the input raster (restored on exit).
6. Read input raster (radiance in W m⁻² sr⁻¹ µm⁻¹).
7. If `elevation=` is given, read it and compute the scene-mean elevation
   (metres → km) as `target_elevation`.
8. Convert radiance to TOA reflectance:
   ```
   rho_toa = (pi * L * d2) / (E0 * cos(sza))
   ```
   where *d*² = Earth–Sun distance [AU²] from `doy=`, and *E*₀ = Thuillier
   solar irradiance at `wavelength=`.
9. Build per-pixel AOD/H₂O fields — priority chain:
   - AOD: `aod_map=` > `visibility=` (Koschmieder) > `aod_val=` > `visibility_val=` (Koschmieder) > conditions file > mid-LUT
   - H₂O: `h2o_map=` > `h2o_val=` > mid-LUT
10. Interpolate the LUT at each pixel's (AOD, H₂O) values (bilinear); or use
    the scene-average slice when no per-pixel maps are supplied.
11. Invert the Lambertian forward model to obtain surface (BOA) reflectance:
    ```
    rho_boa = (rho_toa - R_atm) / (T_down * T_up + s_alb * (rho_toa - R_atm))
    ```
12. Write output raster (float reflectance).

---

## PARAMETER CHOICES

### A. Atmospheric model (`atmosphere=`)

| Value | Meaning |
|---|---|
| `none` | No gaseous absorption |
| `tropical` | Tropical |
| `midsum` | Mid-latitude summer |
| `midwin` | Mid-latitude winter |
| `subarctsum` | Sub-arctic summer |
| `suarctwint` | Sub-arctic winter |
| `us62` | US Standard 1962 **(default)** |

Corresponds to codes 0–6 in the i.atcorr conditions file.

### B. Aerosol model (`aerosol=`)

| Value | Meaning |
|---|---|
| `none` | No aerosols |
| `continental` | Continental **(default)** |
| `maritime` | Maritime |
| `urban` | Urban |
| `desert` | Background desert (Shettle) |
| `biomass` | Biomass burning |
| `stratospheric` | Stratospheric |

Corresponds to codes 0–6 in the i.atcorr conditions file.
Custom aerosol models (codes 7–11 in i.atcorr) are not supported; use
*i.hyper.atcorr* directly for those cases.

### C. AOD LUT grid (`aod=`) and scene value (`aod_val=`)

`aod=` defines the aerosol optical depth grid points at 550 nm at which the
LUT is computed (comma-separated list, e.g. `aod=0.0,0.1,0.2,0.5`).  The
denser the grid, the more accurate the bilinear interpolation but the longer
the LUT computation.

`aod_val=` (or `aod_map=` for per-pixel values) selects the actual AOD used
for correction.  The priority chain is: `aod_map=` > `visibility=`
(Koschmieder) > `aod_val=` > `visibility_val=` (Koschmieder) > value from
`parameters=` file > median of `aod=` grid.

### D. H₂O LUT grid (`h2o=`) and scene value (`h2o_val=`)

`h2o=` defines the column water vapour grid in g/cm² (e.g. `h2o=1.0,2.0,3.5`).
`h2o_val=` (or `h2o_map=`) sets the value used for pixel-level correction.

### E. Sensor and band (`sensor=`, `band=`)

When `sensor=` and `band=` are both provided, the SRF-weighted centre
wavelength is looked up automatically from the sensor CSV files in
`sensors_csv/`.  The 24 supported sensor keys are:

`sentinel2a`, `sentinel2b`, `landsat7_etm`, `landsat8`, `spot6`, `spot7`,
`pleiades1a`, `pleiades1b`, `worldview2`, `worldview3`, `worldview4`,
`geoeye1`, `quickbird2`, `ikonos`, `rapideye`, `planetscope_0c_0d`,
`planetscope_0e`, `planetscope_0f_10`, `avnir`, `vgt1_spot4`, `vgt2_spot5`,
`prism_b`, `prism_f`, `prism_n`.

`wavelength=` overrides `sensor=`/`band=` when given explicitly.

### F. Target elevation and sensor altitude

`target_elevation=` is the scene-average ground elevation in km above sea
level (default 0).  When `elevation=` (a raster map in **metres**,
`G_OPT_R_ELEV` standard) is provided, its spatial mean is converted to km
and used as `target_elevation`.

`altitude=` is the sensor altitude in km.  Values > 900 are treated as
satellite orbit (default 1000).  Set to 0 for ground-level instruments.

### G. Visibility (`visibility=`, `visibility_val=`)

Per-pixel (`visibility=`) or scene-level (`visibility_val=`) meteorological
visibility in km is converted to AOD at 550 nm via the Koschmieder
approximation:

```
AOD_550 ≈ 3.912 / V_km − 0.01162
```

clipped to [0.001, 5.0].  These inputs are overridden by `aod_map=` /
`aod_val=` when provided.

### H. LUT file (`lut=`)

Path to a binary LUT file compatible with the i.hyper.atcorr LUT format
(magic `0x4C555400`, version 1).  If the file exists it is loaded; otherwise
the LUT is computed and written to that path.  This allows the same LUT to be
shared across multiple *i.atcorr2* runs (different bands of the same scene)
or with the *i.hyper.atcorr* C module.

---

## NOTES

### Locating the `grass_sixsv` Python API

The module searches for `atcorr.py` in the following locations, in order:

1. `$GISBASE/scripts/` — system or `g.extension` install
2. `$GRASS_ADDON_BASE/scripts/` — per-user `g.extension` install
3. `../libsixsv/python/` — sibling clone of [github.com/YannChemin/libsixsv](https://github.com/YannChemin/libsixsv)
4. Directory in the `LIBSIXSV_PYTHON` environment variable

Clone and build [libsixsv](https://github.com/YannChemin/libsixsv) with
`make && make install` before running this module.

### Radiance input

The input raster must contain calibrated radiance values in
W m⁻² sr⁻¹ µm⁻¹.  These are converted to TOA reflectance internally using
the Thuillier solar spectrum and the Earth–Sun distance for the given `doy=`.
`doy=` is therefore always required.

### Backward compatibility with i.atcorr parameter files

The `parameters=` option accepts the same 6S conditions text file used by
*i.atcorr*.  Solar geometry (SZA approximated from latitude, longitude, date
and UTC hour), atmospheric model code, aerosol model code and AOD/visibility
are extracted from the file.  Any value explicitly provided as a GRASS option
overrides the file value.

The spectral band code (*iwave* line in the file) is parsed but not used for
automatic wavelength lookup; use `sensor=`/`band=` or explicit `wavelength=`.

### LUT reuse across multiple bands

When correcting several bands of the same scene, compute and save the LUT
on the first call.  Subsequent calls load it instantly:

```sh
# Step 1: build the LUT and correct the red band
i.atcorr2 input=landsat8_B4_rad output=landsat8_B4_boa \
          sensor=landsat8 band=LC08_B4 sza=35.0 doy=180 \
          aod=0.0,0.1,0.2,0.5  h2o=1.0,2.0,3.5 \
          lut=/tmp/landsat8.lut

# Step 2: reuse the LUT for NIR (wavelength must be present in the LUT)
i.atcorr2 input=landsat8_B5_rad output=landsat8_B5_boa \
          sensor=landsat8 band=LC08_B5 sza=35.0 doy=180 \
          lut=/tmp/landsat8.lut
```

For full hyperspectral cubes (many bands at once), use *i.hyper.atcorr*
directly on a Raster3D input.

### Per-pixel atmospheric maps

When `aod_map=`, `visibility=`, or `h2o_map=` rasters are provided the
correction applies bilinear LUT interpolation independently for each pixel.
All maps must cover the same region as the input raster.  Common sources
include MODIS AOD (MOD04) and precipitable water (MOD05) reprojected with
*r.import*, or meteorological visibility fields from NWP model output.

---

## EXAMPLES

### Sentinel-2 with sensor/band auto-lookup

```sh
i.atcorr2 \
    input=sentinel2_B4_rad output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 vza=4.0 raa=95.0 doy=180 \
    aod=0.05,0.15,0.30  h2o=1.0,2.5 \
    aod_val=0.12        h2o_val=1.8 \
    atmosphere=us62 aerosol=continental
```

### Landsat 8 radiance input

```sh
i.atcorr2 \
    input=landsat8_B4_rad output=landsat8_B4_boa \
    sensor=landsat8 band=LC08_B4 doy=210 \
    sza=35.0 vza=0.0 raa=180.0 \
    aod=0.0,0.2,0.5  h2o=1.0,3.0
```

### Per-pixel visibility raster (converted to AOD)

```sh
i.atcorr2 \
    input=sentinel2_B4_rad output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 doy=180 visibility=vis_km_map \
    aod=0.0,0.1,0.2,0.5  h2o=1.0,2.5
```

### Per-pixel AOD and H₂O maps (HLS)

```sh
i.atcorr2 \
    input=hls_red_rad output=hls_red_boa \
    wavelength=0.640 sza=30.0 doy=200 \
    aod_map=aod_550nm  h2o_map=tcw \
    aod=0.0,0.1,0.2,0.4,0.8  h2o=0.5,1.5,3.0,5.0 \
    lut=/tmp/hls.lut
```

### Backward-compatible i.atcorr parameter file

```sh
i.atcorr2 \
    input=etm_band4_rad output=etm_band4_boa \
    parameters=ETM4_atmospheric_input.txt \
    wavelength=0.776 doy=150 \
    aod=0.0,0.1,0.15,0.3  h2o=1.0,2.0
```

### Polarization-enabled correction

```sh
i.atcorr2 -P \
    input=sentinel2_B2_rad output=sentinel2_B2_boa_polar \
    sensor=sentinel2a band=B2 \
    sza=30.0 vza=5.0 raa=90.0 doy=180 \
    aod=0.0,0.1,0.3  h2o=1.0,2.5
```

---

## SEE ALSO

- i.atcorr — classic C++ 6S module
- i.hyper.atcorr — 6SV2.1 hyperspectral module (Raster3D)
- i.landsat.toar — Landsat TOA reflectance conversion
- i.sentinel.import — Sentinel-2 scene import
- r.import — raster import with reprojection
- m.proj — coordinate reprojection

---

## REFERENCES

- Vermote, E.F., Tanré, D., Deuzé, J.L., Herman, M., & Morcrette, J.J.
  (1997). Second simulation of the satellite signal in the solar spectrum,
  6S: An overview. *IEEE Transactions on Geoscience and Remote Sensing*,
  35(3), 675–686.
- Kotchenova, S.Y., Vermote, E.F., Matarrese, R., & Klemm, F.J. (2006).
  Validation of a vector version of the 6S radiative transfer code for
  atmospheric correction of satellite data. *Applied Optics*, 45(26),
  6762–6774.
- Thuillier, G., et al. (2003). The solar spectral irradiance from 200 to
  2400 nm as measured by the SOLSPEC spectrometer from the Atlas and Eureca
  missions. *Solar Physics*, 214(1), 1–22.

---

## AUTHORS

Yann Chemin — based on
[i.atcorr](../grass/imagery/i.atcorr/) (Christo Zietsman / GRASS Development
Team) and *i.hyper.atcorr*.
