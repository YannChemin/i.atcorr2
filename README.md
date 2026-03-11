# i.atcorr2

Atmospheric correction for single-band rasters using the
[`grass_sixsv`](https://github.com/yannchemin/libsixsv) Python bindings
(6SV2.1 LUT-based engine).

---

## i.atcorr vs i.atcorr2

| | **i.atcorr** | **i.atcorr2** |
|---|---|---|
| **Language** | C++ | Python |
| **Engine** | 6S (C++ port of Fortran, pixel-by-pixel) | 6SV2.1 via `grass_sixsv` (LUT-based) |
| **Parameter input** | 6S conditions text file (`parameters=`, required) | Explicit GRASS options; 6S file optionally accepted for backward compatibility |
| **Spectral band** | Sensor codes (iwave −2 to 33) or wl range in conditions file | `sensor=`/`band=` SRF lookup (24 sensors) **or** explicit `wavelength=` µm |
| **Atmospheric model** | Codes 0–6 in the conditions file | Named options: `us62`, `tropical`, `midsum`, … |
| **Aerosol concentration** | Visibility (km) or AOD in conditions file | `aod=` LUT grid + `aod_val=` scene value |
| **Water vapour** | Fixed per-scene from atmospheric model | `h2o=` LUT grid + `h2o_val=` scene value |
| **Per-pixel altitude** | `elevation=` raster (metres) | `elevation=` raster (metres); scene mean converted to km for 6S |
| **Per-pixel visibility** | `visibility=` raster (km) | `visibility=` raster (km) or `visibility_val=` scalar; both converted to AOD via Koschmieder approximation |
| **Per-pixel AOD** | — | `aod_map=` raster |
| **Per-pixel H₂O** | — | `h2o_map=` raster |
| **LUT caching** | RB-tree keyed on (alt, vis) bins | Full 3-D `[AOD × H₂O × WL]` grid, saveable to `lut=` file |
| **Multi-band reuse** | Each band reruns the full 6S computation | Compute LUT once, reuse across bands via `lut=` file |
| **Polarization** | No | `-P` flag (Stokes I, Q, U via 6SV2.1 vector RT) |
| **Radiance → TOA** | Handled internally by 6S transform | Explicit: `ρ_toa = (π·L·d²)/(E₀·cos θs)`; requires `doy=` |
| **Input DN scaling** | `range=min,max` (integer) | `range=min,max` (integer) |
| **Output** | Float or integer (`-i` flag); rescalable | Float reflectance [0, 1] |
| **Dependencies** | GRASS GIS (C build) | GRASS GIS + `grass_sixsv` + NumPy |

---

## When to use which

**Use i.atcorr** when:
- You have a 6S conditions file already prepared.
- You need the classic sensor iwave codes (Landsat TM/ETM+, ASTER, IRS, etc.).
- You cannot install the `grass_sixsv` shared library.

**Use i.atcorr2** when:
- You want a named sensor/band lookup (`sensor=sentinel2a band=B4`) without writing a conditions file.
- You want to correct multiple bands sharing the same scene geometry without rerunning the full RT computation (save LUT once, reuse).
- You have spatially variable AOD or water vapour rasters (e.g. from MODIS, Sentinel-5P).
- You want explicit, scriptable parameters instead of a text conditions file.
- You need 6SV2.1 vector (polarized) radiative transfer.

---

## Quick-start

```sh
# Sensor/band auto-lookup (SRF-weighted centre wavelength)
i.atcorr2 -r \
    input=sentinel2_B4 output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 vza=4.0 raa=95.0 \
    aod=0.05,0.15,0.30  h2o=1.0,2.5 \
    aod_val=0.12        h2o_val=1.8 \
    atmosphere=us62 aerosol=continental

# Explicit wavelength (TOA reflectance input)
i.atcorr2 -r \
    input=sentinel2_B4 output=sentinel2_B4_boa \
    wavelength=0.665 \
    sza=28.5 vza=4.0 raa=95.0 \
    aod=0.05,0.15,0.30  h2o=1.0,2.5 \
    aod_val=0.12        h2o_val=1.8 \
    atmosphere=us62 aerosol=continental

# Radiance input – supply day of year for E₀/d² conversion
i.atcorr2 \
    input=landsat8_B4_dn output=landsat8_B4_boa \
    sensor=landsat8 band=LC08_B4 doy=210 \
    sza=35.0 range=0,65535 \
    aod=0.0,0.2,0.5 h2o=1.0,3.0

# Compute LUT once, reuse for a second band
i.atcorr2 -r \
    input=landsat8_B4 output=landsat8_B4_boa \
    sensor=landsat8 band=LC08_B4 sza=35.0 \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.0 \
    lut=/tmp/scene.lut

i.atcorr2 -r \
    input=landsat8_B5 output=landsat8_B5_boa \
    sensor=landsat8 band=LC08_B5 sza=35.0 \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.0 \
    lut=/tmp/scene.lut          # loaded, not recomputed

# Per-pixel visibility raster (converted to AOD via Koschmieder)
i.atcorr2 -r \
    input=sentinel2_B4 output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 visibility=vis_km_map \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.5

# Per-pixel AOD and H₂O maps
i.atcorr2 -r \
    input=hls_red output=hls_red_boa \
    wavelength=0.640 sza=30.0 \
    aod_map=aod_550nm h2o_map=tcw \
    aod=0.0,0.1,0.2,0.4,0.8 h2o=0.5,1.5,3.0,5.0 \
    lut=/tmp/hls.lut

# Backward-compatible 6S conditions file
i.atcorr2 \
    input=etm_band4 output=etm_band4_boa \
    parameters=ETM4_atmospheric_input.txt \
    wavelength=0.776 \
    aod=0.0,0.1,0.15,0.3 h2o=1.0,2.0
```

---

## Supported sensors (`sensor=`)

24 multispectral sensors are supported via SRF CSV files in `sensors_csv/`.
Band centre wavelengths are computed as the SRF-weighted mean; FWHM values
are derived from the half-maximum criterion or official datasheet overrides.

| Sensor key | Satellite |
|---|---|
| `sentinel2a`, `sentinel2b` | Sentinel-2A/B (MSI) |
| `landsat7_etm`, `landsat8` | Landsat 7 ETM+, Landsat 8 OLI |
| `spot6`, `spot7` | SPOT-6/7 |
| `pleiades1a`, `pleiades1b` | Pléiades-1A/B |
| `worldview2`, `worldview3`, `worldview4` | WorldView-2/3/4 |
| `geoeye1` | GeoEye-1 |
| `quickbird2` | QuickBird-2 |
| `ikonos` | IKONOS |
| `rapideye` | RapidEye |
| `planetscope_0c_0d`, `planetscope_0e`, `planetscope_0f_10` | PlanetScope |
| `avnir` | ALOS AVNIR-2 |
| `vgt1_spot4`, `vgt2_spot5` | SPOT-4/5 VEGETATION |
| `prism_b`, `prism_f`, `prism_n` | ALOS PRISM (backward/forward/nadir) |

Use `python3 -c "import sensors; print(sensors.list_bands('sentinel2a'))"` to
list band names for a given sensor.

---

## Dependencies

| Dependency | Required | Where to get it |
|---|---|---|
| GRASS GIS ≥ 8.0 | yes | https://grass.osgeo.org |
| `grass_sixsv` shared library | yes | https://github.com/yannchemin/libsixsv |
| `atcorr.py` Python bindings | yes | Installed by `libsixsv` `make install` |
| NumPy | yes | `pip install numpy` or system package |

The module searches for `atcorr.py` in this order:

1. `$GISBASE/scripts/` — system or `g.extension` install
2. `$GRASS_ADDON_BASE/scripts/` — per-user `g.extension` install
3. `../libsixsv/python/` — sibling clone of [github.com/yannchemin/libsixsv](https://github.com/yannchemin/libsixsv)
4. Directory in `$IHYPER_ATCORR_PYTHON` environment variable

Set `IHYPER_ATCORR_PYTHON` to override the search path.

---

## Running tests

```sh
cd testsuite
python3 -m grass.gunittest.main --grassdata /path/to/grassdata \
                                --location nc_spm_08_grass7 \
                                --location-type nc
```

The test suite has two classes:

- **`TestIAtcorr2Standalone`** (12 tests) — self-contained physics and
  self-consistency tests; no dependency on `i.atcorr`.
- **`TestIAtcorr2VsAtcorr`** (5 tests) — cross-comparison against `i.atcorr`;
  automatically skipped when `i.atcorr` is not installed.

---

## See also

- [i.atcorr](../grass/imagery/i.atcorr/) — classic C++ 6S module
- [i.hyper.atcorr](../i.hyper.atcorr/) — 6SV2.1 hyperspectral module (Raster3D)
