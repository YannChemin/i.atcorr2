"""Satellite sensor spectral band centre-wavelength and FWHM lookup for i.atcorr2.

Reads the spectral response function (SRF) CSV files shipped in
``sensors_csv/`` and computes the SRF-weighted centre wavelength and FWHM
for each band.  Official manufacturer/agency values override the SRF-computed
FWHM where the half-maximum criterion is unreliable (see :data:`_FWHM_OVERRIDES`).

Centre wavelength
-----------------
SRF-weighted mean:

.. math::

    \\lambda_c = \\frac{\\sum_i \\lambda_i \\, r_i}{\\sum_i r_i}

FWHM
----
Width between the two outermost wavelength points where ``r ≥ max(r) / 2``.
Official manufacturer/agency values are used instead of the half-maximum
computation for sensors whose SRF curve is unreliable for this purpose
(extended tails, multi-lobe shape — see :data:`_FWHM_OVERRIDES` for details).

CSV format
----------
All files in ``sensors_csv/`` follow the same layout:

- Row 0: header — first column is the wavelength label (any name), subsequent
  columns are band names.
- Rows 1+: first column = wavelength in nanometres; remaining columns =
  spectral response (0–1); empty cells are treated as 0.
- Duplicate column names (e.g. PRISM CSVs) are disambiguated by appending
  ``_1``, ``_2``, … suffixes.

FWHM sources
------------
Sensors using official values in :data:`_FWHM_OVERRIDES`:

- **Sentinel-2A/B/C**: ESA Spectral Response Functions document, 2021 edition.
- **Landsat 9 OLI-2**: USGS OLI-2 GLAMR-derived band summary (L9_OLI2_RSR.xlsx).
- **Landsat 8 OLI**: USGS OLI spectral band parameter tables.
- **Landsat 7 ETM+**: NASA/USGS band designation tables (range → FWHM).
- **Landsat 4/5 TM**: USGS TM band designation tables (range → FWHM).
- **Landsat 1–5 MSS**: USGS MSS band summary 50%-point ranges.
- **RapidEye**: Planet Labs Spectral Response Curves datasheet, Rev. C (2020).
- **AVNIR-2**: JAXA ALOS Mission Overview instrument specification.
- **QuickBird-2**: DigitalGlobe imagery user guide band limits.  The SRF CSV
  has extended tails far outside the true passband, making the half-maximum
  criterion unreliable.
- **VGT1/VGT2 MIR band**: ESA/CNES VEGETATION instrument specification
  (passband 1580–1750 nm → 170 nm).  The multi-lobe SRF shape causes the
  half-maximum criterion to underestimate the true FWHM (~88 nm vs 170 nm).

Sensors using SRF-computed FWHM (all others):

- SPOT-6/7, Pleiades-1A/B, WorldView-2/3/4, GeoEye-1, IKONOS, PlanetScope,
  VGT1/VGT2 (B0/B2/B3 bands), PRISM-B/F/N, MODIS Terra, ASTER VNIR/SWIR,
  EO-1 ALI — computed from the SRF CSV files in ``sensors_csv/``.
"""

import csv
import os

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_CSV_DIR = os.path.join(_HERE, "sensors_csv")

# Map from user-facing sensor key → CSV filename in sensors_csv/
SENSORS = {
    # ── Sentinel-2 ────────────────────────────────────────────────────────────
    "sentinel2a":        "sentinel_2A_msi.csv",
    "sentinel2b":        "sentinel_2B_msi.csv",
    "sentinel2c":        "sentinel2c.csv",
    # ── Landsat (OLI) ─────────────────────────────────────────────────────────
    "landsat9_oli2":     "landsat9_oli2.csv",
    "landsat8":          "landsat_8.csv",
    # ── Landsat (TM) ──────────────────────────────────────────────────────────
    "landsat5_tm":       "landsat5_tm.csv",
    "landsat4_tm":       "landsat4_tm.csv",
    # ── Landsat (ETM+) ────────────────────────────────────────────────────────
    "landsat7_etm":      "etmplus.csv",
    # ── Landsat (MSS) ─────────────────────────────────────────────────────────
    "landsat5_mss":      "landsat5_mss.csv",
    "landsat4_mss":      "landsat4_mss.csv",
    "landsat3_mss":      "landsat3_mss.csv",
    "landsat2_mss":      "landsat2_mss.csv",
    "landsat1_mss":      "landsat1_mss.csv",
    # ── MODIS / ASTER / EO-1 ──────────────────────────────────────────────────
    "modis_terra":       "modis_terra.csv",
    "aster":             "aster_vnir_swir.csv",
    "eo1_ali":           "eo1_ali.csv",
    # ── Pléiades Neo / Amazônia-1 / CBERS-4A ─────────────────────────────────
    # SRF approximated as flat-top trapezoids from official band-limit tables
    # (USGS OFR 2021-1030-P, -N, -J).  No digitized per-wavelength RSR published.
    "pleiades_neo":      "pleiades_neo.csv",
    "amazonia1":         "amazonia1.csv",
    "cbers4a_mux":       "cbers4a_mux.csv",
    # ── Commercial VHR ────────────────────────────────────────────────────────
    "worldview2":        "worldview2.csv",
    "worldview3":        "worldview3.csv",
    "worldview4":        "worldview4.csv",
    "spot6":             "spot6.csv",
    "spot7":             "spot7.csv",
    "rapideye":          "rapideye.csv",
    "ikonos":            "ikonos.csv",
    "geoeye1":           "geoeye1.csv",
    "quickbird2":        "quickbird2.csv",
    "pleiades1a":        "pleiades1a.csv",
    "pleiades1b":        "pleiades1b.csv",
    "avnir":             "AVNIR.csv",
    # ── SPOT-4/5 VGT ──────────────────────────────────────────────────────────
    "vgt1_spot4":        "VGT1_spot4.csv",
    "vgt2_spot5":        "VGT2_spot5.csv",
    # ── PlanetScope ───────────────────────────────────────────────────────────
    "planetscope_0c_0d":    "planetscope_0c_0d.csv",
    "planetscope_0e":       "planetscope_0e.csv",
    "planetscope_0f_10":    "planetscope_0f_10.csv",
    "planetscope_dove_r":   "planetscope_dove_r.csv",
    "planetscope_superdove": "planetscope_superdove.csv",
    # ── SkySat ─────────────────────────────────────────────────────────────────
    "skysat_1_2":   "skysat_1_2.csv",
    "skysat_3_13":  "skysat_3_13.csv",
    "skysat_14_19": "skysat_14_19.csv",
    # ── PRISM (ALOS-2) ────────────────────────────────────────────────────────
    "prism_b":           "PRISM-B.csv",
    "prism_f":           "PRISM-F.csv",
    "prism_n":           "PRISM-N.csv",
}

# Official FWHM overrides (nm) for sensors where the manufacturer / agency
# provides explicit values that supersede what the half-maximum computation
# gives on the SRF curve.  Format: {sensor_key: {band_name: fwhm_nm}}.
#
# Sentinel-2A/B: ESA S2 Spectral Response Functions document (2021 edition)
#   https://sentinels.copernicus.eu/web/sentinel/user-guides/sentinel-2-msi/document-library
# Landsat 7 ETM+: NASA/USGS band designation tables (centre ± FWHM/2 = wl range)
# Landsat 8 OLI:  USGS OLI spectral band parameters
# RapidEye: Planet Labs Spectral Response Curves datasheet Rev. C (2020)
# AVNIR-2: JAXA ALOS Mission Overview (band ranges → FWHM)
# QuickBird-2: DigitalGlobe imagery user guide (band limits → FWHM)
# VGT1/VGT2 MIR: ESA/CNES VEGETATION instrument spec (1580–1750 nm → 170 nm)

_FWHM_OVERRIDES: dict[str, dict[str, float]] = {
    # ── Sentinel-2A (ESA official) ─────────────────────────────────────────
    "sentinel2a": {
        "S2A_SR_AV_B1":  20.0,   # 443 nm – Coastal aerosol
        "S2A_SR_AV_B2":  65.0,   # 490 nm – Blue
        "S2A_SR_AV_B3":  35.0,   # 560 nm – Green
        "S2A_SR_AV_B4":  30.0,   # 665 nm – Red
        "S2A_SR_AV_B5":  15.0,   # 705 nm – Red Edge 1
        "S2A_SR_AV_B6":  15.0,   # 740 nm – Red Edge 2
        "S2A_SR_AV_B7":  20.0,   # 783 nm – Red Edge 3
        "S2A_SR_AV_B8": 115.0,   # 842 nm – NIR broad
        "S2A_SR_AV_B8A": 20.0,   # 865 nm – NIR narrow
        "S2A_SR_AV_B9":  20.0,   # 945 nm – Water vapour
        "S2A_SR_AV_B10": 30.0,   # 1375 nm – Cirrus
        "S2A_SR_AV_B11": 90.0,   # 1610 nm – SWIR 1
        "S2A_SR_AV_B12":180.0,   # 2190 nm – SWIR 2
    },
    # ── Sentinel-2C (ESA official — from same SRF document as 2A/2B) ──────
    "sentinel2c": {
        "S2C_SR_AV_B1":  20.0,   # 443 nm – Coastal aerosol
        "S2C_SR_AV_B2":  65.0,   # 492 nm – Blue
        "S2C_SR_AV_B3":  35.0,   # 560 nm – Green
        "S2C_SR_AV_B4":  30.0,   # 665 nm – Red
        "S2C_SR_AV_B5":  15.0,   # 704 nm – Red Edge 1
        "S2C_SR_AV_B6":  15.0,   # 739 nm – Red Edge 2
        "S2C_SR_AV_B7":  20.0,   # 780 nm – Red Edge 3
        "S2C_SR_AV_B8": 115.0,   # 842 nm – NIR broad
        "S2C_SR_AV_B8A": 20.0,   # 864 nm – NIR narrow
        "S2C_SR_AV_B9":  20.0,   # 943 nm – Water vapour
        "S2C_SR_AV_B10": 30.0,   # 1375 nm – Cirrus
        "S2C_SR_AV_B11": 90.0,   # 1610 nm – SWIR 1
        "S2C_SR_AV_B12":180.0,   # 2190 nm – SWIR 2
    },
    # ── Sentinel-2B (ESA official — slightly different from 2A) ───────────
    "sentinel2b": {
        "SR_AV_B1":  20.0,   # 442 nm – Coastal aerosol
        "SR_AV_B2":  66.0,   # 492 nm – Blue
        "SR_AV_B3":  36.0,   # 559 nm – Green
        "SR_AV_B4":  31.0,   # 665 nm – Red
        "SR_AV_B5":  16.0,   # 704 nm – Red Edge 1
        "SR_AV_B6":  15.0,   # 739 nm – Red Edge 2
        "SR_AV_B7":  20.0,   # 780 nm – Red Edge 3
        "SR_AV_B8": 106.0,   # 833 nm – NIR broad
        "SR_AV_B8A": 21.0,   # 864 nm – NIR narrow
        "SR_AV_B9":  20.0,   # 943 nm – Water vapour
        "SR_AV_B10": 30.0,   # 1377 nm – Cirrus
        "SR_AV_B11": 94.0,   # 1610 nm – SWIR 1
        "SR_AV_B12":185.0,   # 2186 nm – SWIR 2
    },
    # ── Landsat 7 ETM+ (NASA/USGS band ranges → FWHM) ─────────────────────
    "landsat7_etm": {
        "1 blue":        65.0,   # 450–515 nm
        "2 green":       80.0,   # 525–605 nm
        "3 red":         60.0,   # 630–690 nm
        "4 nir":        125.0,   # 775–900 nm
        "5 swir 1":     200.0,   # 1550–1750 nm
        "7 swir 2":     270.0,   # 2080–2350 nm
        "pan":          380.0,   # 520–900 nm
    },
    # ── Landsat 9 OLI-2 (USGS GLAMR band summary) ─────────────────────────
    "landsat9_oli2": {
        "CoastalAerosol": 20.0,  # 435–451 nm
        "Blue":           60.0,  # 452–512 nm
        "Green":          57.0,  # 533–590 nm
        "Red":            37.0,  # 636–673 nm
        "NIR":            28.0,  # 851–879 nm
        "Cirrus":         21.0,  # 1363–1384 nm
        "SWIR1":          85.0,  # 1566–1651 nm
        "SWIR2":         187.0,  # 2107–2294 nm
        "Pan":           172.0,  # 503–676 nm
    },
    # ── Landsat 8 OLI (USGS band parameters) ──────────────────────────────
    "landsat8": {
        "Coastal Aerosol": 20.0,  # 430–450 nm
        "Blue":            60.0,  # 450–510 nm
        "Green":           60.0,  # 530–590 nm
        "Red":             30.0,  # 640–670 nm
        "Nir":             30.0,  # 850–880 nm
        "Cirrus":          20.0,  # 1360–1380 nm
        "SWIR1":           80.0,  # 1570–1650 nm
        "SWIR2":          180.0,  # 2110–2290 nm
        "Pan":            180.0,  # 500–680 nm
    },
    # ── Landsat 5 TM (USGS band ranges) ───────────────────────────────────
    "landsat5_tm": {
        "Blue":   70.0,   # 450–520 nm
        "Green":  80.0,   # 520–600 nm
        "Red":    60.0,   # 630–690 nm
        "NIR":   130.0,   # 760–900 nm
        "SWIR1": 200.0,   # 1550–1750 nm
        "SWIR2": 270.0,   # 2080–2350 nm
    },
    # ── Landsat 4 TM (same band spec as L5 TM) ────────────────────────────
    "landsat4_tm": {
        "Blue":   70.0,
        "Green":  80.0,
        "Red":    60.0,
        "NIR":   130.0,
        "SWIR1": 200.0,
        "SWIR2": 270.0,
    },
    # ── Landsat 1–3 MSS (USGS band summary 50%-point ranges) ──────────────
    # L1: Green 504–602=98, Red 605–701=96, NIR1 703–793=90, NIR2 802–1094=292
    "landsat1_mss": {
        "Green":  98.0,
        "Red":    96.0,
        "NIR1":   90.0,
        "NIR2":  292.0,
    },
    # L2: Green 498–600=102, Red 607–711=104, NIR1 706–784=78, NIR2 804–1095=291
    "landsat2_mss": {
        "Green": 102.0,
        "Red":   104.0,
        "NIR1":   78.0,
        "NIR2":  291.0,
    },
    # L3: Green 495–593=98, Red 607–704=97, NIR1 702–793=91, NIR2 804–1088=284
    "landsat3_mss": {
        "Green":  98.0,
        "Red":    97.0,
        "NIR1":   91.0,
        "NIR2":  284.0,
    },
    # ── Landsat 4–5 MSS (USGS band summary 50%-point ranges) ──────────────
    # L4 MSS Band 1 Green 500–600, Band 2 Red 600–700, Band 3 NIR1 700–800, Band 4 NIR2 800–1100
    "landsat4_mss": {
        "Green": 100.0,
        "Red":   100.0,
        "NIR1":  100.0,
        "NIR2":  300.0,
    },
    "landsat5_mss": {
        "Green": 100.0,
        "Red":   100.0,
        "NIR1":  100.0,
        "NIR2":  300.0,
    },
    # ── Pléiades Neo (USGS OFR 2021-1030-P, Table 2 — band limits → FWHM) ───
    # SRF is trapezoidal approximation; FWHM = upper – lower band limit.
    "pleiades_neo": {
        "Pan":      350.0,   # 450–800 nm
        "DeepBlue":  60.0,   # 400–460 nm  (coastal aerosol)
        "Blue":      60.0,   # 460–520 nm
        "Green":     80.0,   # 520–600 nm
        "Red":       80.0,   # 600–680 nm
        "RedEdge":   40.0,   # 680–720 nm
        "NIR":      170.0,   # 730–900 nm
    },
    # ── Amazônia-1 WFI-2 (USGS OFR 2021-1030-N, Table 2) ─────────────────
    "amazonia1": {
        "Blue":   70.0,   # 450–520 nm
        "Green":  70.0,   # 520–590 nm
        "Red":    60.0,   # 630–690 nm
        "NIR":   120.0,   # 770–890 nm
    },
    # ── CBERS-4A MUX (USGS OFR 2021-1030-J, Table 2) ─────────────────────
    # Same band limits as Amazônia-1 WFI-2 — identical instrument design.
    "cbers4a_mux": {
        "Band5_Blue":   70.0,   # 450–520 nm
        "Band6_Green":  70.0,   # 520–590 nm
        "Band7_Red":    60.0,   # 630–690 nm
        "Band8_NIR":   120.0,   # 770–890 nm
    },
    # ── RapidEye (Planet Labs datasheet 2020) ─────────────────────────────
    "rapideye": {
        "Blue":     66.0,   # 440–510 nm
        "Green":    71.0,   # 520–590 nm
        "Red":      62.0,   # 630–685 nm  (sometimes quoted 55 nm)
        "Red Edge": 39.0,   # 690–730 nm
        "NIR":      98.0,   # 760–850 nm
    },
    # ── QuickBird-2 (DigitalGlobe imagery user guide — band limits) ───────────
    # ── IKONOS (Space Imaging / GeoEye spectral characterisation report) ──────
    # Trapezoidal CSV generated from official band limits.
    # Official FWHM values taken directly from the characterisation table.
    "ikonos": {
        "Panchromatic": 403.0,  # 525.8–928.5 nm
        "Blue":          71.3,  # 444.7–516.0 nm
        "Green":         88.6,  # 506.4–595.0 nm
        "Red":           65.8,  # 631.9–697.7 nm
        "NIR":           95.4,  # 757.3–852.7 nm
    },
    # SRF CSV tails extend far outside the true passband, making the half-max
    # criterion unreliable.  Official limits used instead:
    #   Pan 450–900 nm, Blue 450–520, Green 520–600, Red 630–690, NIR 760–900.
    "quickbird2": {
        "pan":   450.0,
        "blue":   70.0,
        "green":  80.0,
        "red":    60.0,
        "nir":   140.0,
    },
    # ── VGT1/SPOT-4 MIR (ESA/CNES VEGETATION spec: 1580–1750 nm) ──────────
    # SRF has a complex multi-lobe shape; half-max yields ~88 nm, but the
    # instrument passband is 1580–1750 nm → FWHM ≈ 170 nm.
    "vgt1_spot4": {
        "MIR": 170.0,
    },
    # ── VGT2/SPOT-5 MIR (same instrument design as VGT1) ──────────────────
    "vgt2_spot5": {
        "MIR": 170.0,
    },
    # ── AVNIR-2 / ALOS (JAXA instrument spec) ─────────────────────────────
    "avnir": {
        "band 1": 80.0,    # 420–500 nm
        "band 2": 80.0,    # 520–600 nm
        "band 3": 80.0,    # 610–690 nm
        "band 4":130.0,    # 760–890 nm
    },
}

def _load_srf_csv(sensor_key):
    """Load the SRF CSV for *sensor_key* and return ``(wl_nm, bands)``.

    Parameters
    ----------
    sensor_key : str
        Key from :data:`SENSORS`.

    Returns
    -------
    wl_nm : numpy.ndarray, shape (N,)
        Wavelength axis in nanometres.
    bands : dict[str, numpy.ndarray]
        Mapping of band name → SRF array (same length as *wl_nm*).

    Raises
    ------
    KeyError
        If *sensor_key* is not in :data:`SENSORS`.
    FileNotFoundError
        If the corresponding CSV file is missing.
    """
    csv_path = os.path.join(_CSV_DIR, SENSORS[sensor_key])
    wl_list = []
    band_lists = {}

    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        # Strip surrounding whitespace/quotes from header tokens
        header = [h.strip().strip('"') for h in header]
        # Deduplicate column names (e.g. PRISM CSVs reuse the same name)
        raw_names = header[1:]
        band_names = []
        seen: dict[str, int] = {}
        for name in raw_names:
            if name in seen:
                seen[name] += 1
                name = f"{name}_{seen[name]}"
            else:
                seen[name] = 0
            band_names.append(name)
        for name in band_names:
            band_lists[name] = []

        for row in reader:
            if not row or not row[0].strip():
                continue
            try:
                wl = float(row[0])
            except ValueError:
                continue
            wl_list.append(wl)
            for i, name in enumerate(band_names):
                cell = row[i + 1].strip() if i + 1 < len(row) else ""
                try:
                    band_lists[name].append(float(cell))
                except ValueError:
                    band_lists[name].append(0.0)

    wl_nm = np.array(wl_list, dtype=np.float64)
    bands = {name: np.array(vals, dtype=np.float64)
             for name, vals in band_lists.items()}
    return wl_nm, bands

def _srf_fwhm(wl_nm, srf):
    """Compute FWHM (nm) from an SRF array using the half-maximum criterion.

    Parameters
    ----------
    wl_nm : numpy.ndarray
        Wavelength axis in nm.
    srf : numpy.ndarray
        Spectral response values (same length).

    Returns
    -------
    float or None
        FWHM in nm, or ``None`` if the SRF has no significant response.
    """
    peak = srf.max()
    if peak <= 0:
        return None
    above = wl_nm[srf >= peak / 2.0]
    if len(above) < 2:
        return 0.0  # degenerate single-point response
    return float(above[-1] - above[0])

def load_band_centers(sensor_key):
    """Return SRF-weighted centre wavelengths (µm) for all bands of a sensor.

    The centre wavelength is computed as the SRF-weighted mean:

    .. math::

        \\lambda_c = \\frac{\\sum_i \\lambda_i \\, r_i}{\\sum_i r_i}

    where :math:`r_i` is the spectral response at wavelength
    :math:`\\lambda_i`.  Bands whose total response is zero are omitted.

    Parameters
    ----------
    sensor_key : str
        Key from :data:`SENSORS` (e.g. ``"sentinel2a"``).

    Returns
    -------
    dict[str, float]
        Band name → centre wavelength in µm.
    """
    wl_nm, bands = _load_srf_csv(sensor_key)
    centers = {}
    for name, srf in bands.items():
        total = srf.sum()
        if total > 0:
            centers[name] = float(np.dot(wl_nm, srf) / total) / 1000.0
    return centers

def load_band_fwhm(sensor_key):
    """Return FWHM (nm) for all bands of a sensor.

    Priority order for each band:

    1. Entry in :data:`_FWHM_OVERRIDES` — official manufacturer/agency value.
    2. Half-maximum width computed from the SRF curve in the CSV.

    Parameters
    ----------
    sensor_key : str
        Key from :data:`SENSORS` (e.g. ``"sentinel2a"``).

    Returns
    -------
    dict[str, float]
        Band name → FWHM in nm.  Bands with zero total response are omitted.
    """
    wl_nm, bands = _load_srf_csv(sensor_key)
    overrides = _FWHM_OVERRIDES.get(sensor_key, {})
    result = {}
    for name, srf in bands.items():
        if srf.sum() <= 0:
            continue
        if name in overrides:
            result[name] = overrides[name]
        else:
            fw = _srf_fwhm(wl_nm, srf)
            if fw is not None:
                result[name] = fw
    return result

def load_band_specs(sensor_key):
    """Return centre wavelength (µm) and FWHM (nm) for all bands of a sensor.

    Combines :func:`load_band_centers` and :func:`load_band_fwhm` into a
    single call.

    Parameters
    ----------
    sensor_key : str
        Key from :data:`SENSORS`.

    Returns
    -------
    dict[str, tuple[float, float]]
        Band name → ``(centre_um, fwhm_nm)``, sorted by centre wavelength.
    """
    centers = load_band_centers(sensor_key)
    fwhms = load_band_fwhm(sensor_key)
    specs = {name: (centers[name], fwhms[name])
             for name in centers if name in fwhms}
    return dict(sorted(specs.items(), key=lambda kv: kv[1][0]))

def find_sensor_for_band(band_name):
    """Find the sensor key(s) that contain *band_name* in their SRF table.

    Searches all sensor CSVs and returns the unique sensor key whose band list
    includes *band_name*.  Raises :exc:`ValueError` when the name is ambiguous
    (found in more than one sensor) or not found in any sensor.

    Parameters
    ----------
    band_name : str
        Exact band column name as it appears in the sensor SRF CSV header.

    Returns
    -------
    str
        The unique sensor key (from :data:`SENSORS`) that contains *band_name*.

    Raises
    ------
    ValueError
        If *band_name* is not found in any sensor, or if it is found in more
        than one sensor (use ``sensor=`` explicitly in that case).
    """
    matches = []
    for key in SENSORS:
        try:
            centers = load_band_centers(key)
        except Exception:
            continue
        if band_name in centers:
            matches.append(key)
    if len(matches) == 1:
        return matches[0]
    if len(matches) == 0:
        known = sorted(SENSORS.keys())
        raise ValueError(
            f"Band '{band_name}' not found in any sensor SRF table. "
            f"Available sensors: {', '.join(known)}. "
            "Use list_bands(sensor_key) to see band names for a sensor."
        )
    raise ValueError(
        f"Band '{band_name}' is ambiguous — found in sensors: "
        f"{', '.join(sorted(matches))}. "
        "Specify sensor= explicitly."
    )


def list_sensors():
    """Return the list of available sensor keys.

    Returns
    -------
    list[str]
    """
    return sorted(SENSORS.keys())

def list_bands(sensor_key):
    """Return band names with centre wavelength and FWHM for a sensor.

    Parameters
    ----------
    sensor_key : str
        Key from :data:`SENSORS`.

    Returns
    -------
    dict[str, tuple[float, float]]
        Band name → ``(centre_um, fwhm_nm)``, sorted by centre wavelength.
    """
    return load_band_specs(sensor_key)
