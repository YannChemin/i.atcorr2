# i.atcorr2

Atmospheric correction for single-band rasters using the
[`grass_sixsv`](https://github.com/YannChemin/libsixsv) Python bindings
(6SV2.1 LUT-based engine).

---

## i.atcorr vs i.atcorr2

| | **i.atcorr** | **i.atcorr2** |
|---|---|---|
| **Language** | C++ | Python |
| **Engine** | 6S (C++ port of Fortran, pixel-by-pixel) | 6SV2.1 via `grass_sixsv` (LUT-based) |
| **Parameter input** | 6S conditions text file (`parameters=`, required) | Explicit GRASS options; 6S file optionally accepted for backward compatibility |
| **Spectral band** | Sensor codes (iwave −2 to 33) or wl range in conditions file | `sensor=`/`band=` SRF lookup (44 sensors) **or** explicit `wavelength=` µm |
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
| **Input** | Radiance or reflectance (`-r` flag) | Radiance (W m⁻² sr⁻¹ µm⁻¹) only; `doy=` required |
| **Radiance → TOA** | Handled internally by 6S transform | Explicit: `ρ_toa = (π·L·d²)/(E₀·cos θs)` using Thuillier spectrum |
| **Output** | Float or integer (`-i` flag); rescalable | Float reflectance [0, 1] |
| **Dependencies** | GRASS GIS (C build) | GRASS GIS + `grass_sixsv` + NumPy |

---

## When to use which

**Use i.atcorr** when:
- You have a 6S conditions file already prepared.
- You need the classic sensor iwave codes (Landsat TM/ETM+, ASTER, IRS, etc.).
- You cannot install the `grass_sixsv` shared library.

**Use i.atcorr2** when:
- Your input is calibrated radiance (W m⁻² sr⁻¹ µm⁻¹).
- You want a named sensor/band lookup (`sensor=sentinel2a band=B4`) without writing a conditions file.
- You want to correct multiple bands sharing the same scene geometry without rerunning the full RT computation (save LUT once, reuse).
- You have spatially variable AOD or water vapour rasters (e.g. from MODIS, Sentinel-5P).
- You want explicit, scriptable parameters instead of a text conditions file.
- You need 6SV2.1 vector (polarized) radiative transfer.

---

## Quick-start

```sh
# Sensor/band auto-lookup (SRF-weighted centre wavelength)
i.atcorr2 \
    input=sentinel2_B4_rad output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 vza=4.0 raa=95.0 doy=180 \
    aod=0.05,0.15,0.30  h2o=1.0,2.5 \
    aod_val=0.12        h2o_val=1.8 \
    atmosphere=us62 aerosol=continental

# Explicit wavelength
i.atcorr2 \
    input=sentinel2_B4_rad output=sentinel2_B4_boa \
    wavelength=0.665 \
    sza=28.5 vza=4.0 raa=95.0 doy=180 \
    aod=0.05,0.15,0.30  h2o=1.0,2.5 \
    aod_val=0.12        h2o_val=1.8 \
    atmosphere=us62 aerosol=continental

# Compute LUT once, reuse for a second band
i.atcorr2 \
    input=landsat8_B4_rad output=landsat8_B4_boa \
    sensor=landsat8 band=LC08_B4 sza=35.0 doy=210 \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.0 \
    lut=/tmp/scene.lut

i.atcorr2 \
    input=landsat8_B5_rad output=landsat8_B5_boa \
    sensor=landsat8 band=LC08_B5 sza=35.0 doy=210 \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.0 \
    lut=/tmp/scene.lut          # loaded, not recomputed

# Per-pixel visibility raster (converted to AOD via Koschmieder)
i.atcorr2 \
    input=sentinel2_B4_rad output=sentinel2_B4_boa \
    sensor=sentinel2a band=B4 \
    sza=28.5 doy=180 visibility=vis_km_map \
    aod=0.0,0.1,0.2,0.5 h2o=1.0,2.5

# Per-pixel AOD and H₂O maps
i.atcorr2 \
    input=hls_red_rad output=hls_red_boa \
    wavelength=0.640 sza=30.0 doy=200 \
    aod_map=aod_550nm h2o_map=tcw \
    aod=0.0,0.1,0.2,0.4,0.8 h2o=0.5,1.5,3.0,5.0 \
    lut=/tmp/hls.lut

# Backward-compatible 6S conditions file
i.atcorr2 \
    input=etm_band4_rad output=etm_band4_boa \
    parameters=ETM4_atmospheric_input.txt \
    wavelength=0.776 doy=150 \
    aod=0.0,0.1,0.15,0.3 h2o=1.0,2.0
```

---

## Listing band names (`-l` flag)

To see the available band names, centre wavelengths, and FWHM for any sensor,
run with the `-l` flag:

```sh
# List bands for a specific sensor
i.atcorr2 -l sensor=sentinel2a

# List all sensor keys (omit sensor=)
i.atcorr2 -l
```

Example output for `sensor=landsat9_oli2`:

```
Bands for sensor 'landsat9_oli2':
  Band name                   Centre (µm)   FWHM (nm)
  -------------------------  ------------  ----------
  CoastalAerosol                   0.4428        20.0
  Blue                             0.4823        60.0
  Green                            0.5609        57.0
  Pan                              0.5939       172.0
  Red                              0.6543        37.0
  NIR                              0.8646        28.0
  Cirrus                           1.3740        21.0
  SWIR1                            1.6084        85.0
  SWIR2                            2.2011       187.0
```

---

## Supported sensors (`sensor=`)

44 multispectral sensors are supported via SRF CSV files in `sensors_csv/`.
Band centre wavelengths are computed as the SRF-weighted mean; FWHM values
are derived from the half-maximum criterion or official datasheet overrides.
Sensors marked † use trapezoidal SRF approximations derived from official
band-limit tables (no digitized per-wavelength RSR was published).

| Sensor key | Satellite / instrument |
|---|---|
| `sentinel2a`, `sentinel2b`, `sentinel2c` | Sentinel-2A/B/C (MSI) |
| `landsat9_oli2`, `landsat8` | Landsat 9 OLI-2, Landsat 8 OLI |
| `landsat7_etm` | Landsat 7 ETM+ |
| `landsat5_tm`, `landsat4_tm` | Landsat 4/5 TM |
| `landsat5_mss`, `landsat4_mss`, `landsat3_mss`, `landsat2_mss`, `landsat1_mss` | Landsat 1–5 MSS |
| `modis_terra` | MODIS Terra (bands 1–19, solar reflective) |
| `aster` | Terra ASTER (VNIR + SWIR, bands 1–9) |
| `eo1_ali` | EO-1 ALI |
| `pleiades_neo` † | Pléiades Neo (7 bands incl. coastal aerosol + red edge) |
| `amazonia1` † | Amazônia-1 WFI-2 (INPE, Brazil) |
| `cbers4a_mux` † | CBERS-4A MUX (China-Brazil, same bands as Amazônia-1) |
| `spot6`, `spot7` | SPOT-6/7 |
| `pleiades1a`, `pleiades1b` | Pléiades-1A/B |
| `worldview2`, `worldview3`, `worldview4` | WorldView-2/3/4 |
| `geoeye1` | GeoEye-1 |
| `quickbird2` | QuickBird-2 |
| `ikonos` | IKONOS |
| `rapideye` | RapidEye |
| `planetscope_0c_0d`, `planetscope_0e`, `planetscope_0f_10` | PlanetScope Dove Classic (4-band) |
| `planetscope_dove_r` | PlanetScope Dove-R / PS2.SD (4-band) |
| `planetscope_superdove` | PlanetScope SuperDove / PSB.SD (8-band) |
| `skysat_1_2`, `skysat_3_13`, `skysat_14_19` | Planet SkySat (5-band incl. Pan; by generation) |
| `avnir` | ALOS AVNIR-2 |
| `vgt1_spot4`, `vgt2_spot5` | SPOT-4/5 VEGETATION |
| `prism_b`, `prism_f`, `prism_n` | ALOS PRISM (backward/forward/nadir) |

Run `i.atcorr2 -l sensor=<key>` to list band names for any sensor.

---

## Dependencies

| Dependency | Required | Where to get it |
|---|---|---|
| GRASS GIS ≥ 8.0 | yes | https://grass.osgeo.org |
| `grass_sixsv` shared library | yes | https://github.com/YannChemin/libsixsv |
| `atcorr.py` Python bindings | yes | Installed by `libsixsv` `make install` |
| NumPy | yes | `pip install numpy` or system package |

The module searches for `atcorr.py` in this order:

1. `$GISBASE/scripts/` — system or `g.extension` install
2. `$GRASS_ADDON_BASE/scripts/` — per-user `g.extension` install
3. `../libsixsv/python/` — sibling clone of [github.com/YannChemin/libsixsv](https://github.com/YannChemin/libsixsv)
4. Directory in `$LIBSIXSV_PYTHON` environment variable

Set `LIBSIXSV_PYTHON` to override the search path.

---

## Running tests

```sh
cd testsuite
python3 -m grass.gunittest.main --grassdata /path/to/grassdata \
                                --location nc_spm_08_grass7 \
                                --location-type nc
```

The test suite has two classes:

- **`TestIAtcorr2Standalone`** (11 tests) — self-contained physics and
  self-consistency tests using synthetic radiance input; no dependency on `i.atcorr`.
- **`TestIAtcorr2VsAtcorr`** (5 tests) — cross-comparison against `i.atcorr`;
  automatically skipped when `i.atcorr` is not installed.

---

## See also

- [i.atcorr](../grass/imagery/i.atcorr/) — classic C++ 6S module
- [libsixsv](https://github.com/YannChemin/libsixsv) — 6SV2.1 shared library
