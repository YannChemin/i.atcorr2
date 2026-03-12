#!/usr/bin/env python3
"""Process new SRF source files and generate output CSVs."""

import os
import numpy as np
import pandas as pd

CSV_DIR = os.path.dirname(os.path.abspath(__file__))


def resample_1nm(wl_nm, df_bands):
    """Resample band columns to 1 nm integer steps."""
    wl_new = np.arange(int(np.ceil(wl_nm.min())), int(np.floor(wl_nm.max())) + 1)
    result = {"WL_nm": wl_new}
    for col in df_bands.columns:
        result[col] = np.interp(wl_new, wl_nm, df_bands[col].values)
    return pd.DataFrame(result)


def write_csv(df, out_path):
    """Write DataFrame to CSV: WL_nm as int, band values as %.4f."""
    with open(out_path, "w") as fh:
        fh.write(",".join(df.columns) + "\n")
        for _, row in df.iterrows():
            wl = int(row["WL_nm"])
            vals = ",".join(f"{row[c]:.4f}" for c in df.columns if c != "WL_nm")
            fh.write(f"{wl},{vals}\n")
    print(f"  Written: {out_path}")


# ── PlanetScope 0C/0D, 0E, 0F/10 ────────────────────────────────────────────
for src, out in [
    ("PlanetScope_RSR_SatID_0c_0d.csv", "planetscope_0c_0d.csv"),
    ("PlanetScope_RSR_SatID_0e.csv",    "planetscope_0e.csv"),
    ("PlanetScope_RSR_SatID_0f_10.csv", "planetscope_0f_10.csv"),
]:
    print(f"\nProcessing {src} → {out}")
    df = pd.read_csv(os.path.join(CSV_DIR, src), encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    # First column is wavelength in µm → multiply by 1000 to get nm
    wl_col = df.columns[0]
    wl_nm = df[wl_col].values * 1000.0
    # Remaining columns: Blue, Green, Red, NIR
    band_cols = df.columns[1:]
    bands_df = df[band_cols].copy()
    bands_df.columns = ["Blue", "Green", "Red", "NIR"]
    # Round wl_nm to integers (already on 10nm grid but need exact integers)
    wl_nm = np.round(wl_nm).astype(float)
    out_df = pd.DataFrame({"WL_nm": wl_nm.astype(int)})
    for col in bands_df.columns:
        out_df[col] = bands_df[col].values
    write_csv(out_df, os.path.join(CSV_DIR, out))

# ── Dove-R ───────────────────────────────────────────────────────────────────
print("\nProcessing PlanetScope_RSR_dove_r.csv → planetscope_dove_r.csv")
df = pd.read_csv(os.path.join(CSV_DIR, "PlanetScope_RSR_dove_r.csv"))
df.columns = [c.strip() for c in df.columns]
wl_col = df.columns[0]
wl_nm = df[wl_col].values.astype(float)
band_cols = df.columns[1:]
bands_df = df[band_cols].copy()
bands_df.columns = ["Blue", "Green", "Red", "NIR"]
# Already in nm — build output directly (already integer steps)
out_df = pd.DataFrame({"WL_nm": wl_nm.astype(int)})
for col in bands_df.columns:
    out_df[col] = bands_df[col].values
write_csv(out_df, os.path.join(CSV_DIR, "planetscope_dove_r.csv"))

# ── SuperDove ────────────────────────────────────────────────────────────────
print("\nProcessing PlanetScope_RSR_Superdove.csv → planetscope_superdove.csv")
df = pd.read_csv(os.path.join(CSV_DIR, "PlanetScope_RSR_Superdove.csv"))
df.columns = [c.strip() for c in df.columns]
wl_col = df.columns[0]
wl_nm = df[wl_col].values.astype(float)
# Band name mapping
band_map = {
    "Coastal-Blue response": "CoastalBlue",
    "Blue response": "Blue",
    "Green_i response": "Green_i",
    "Green_ii response": "Green_ii",
    "Yellow response": "Yellow",
    "Red response": "Red",
    "Red-edge response": "RedEdge",
    "NIR response": "NIR",
}
out_df = pd.DataFrame({"WL_nm": wl_nm.astype(int)})
for src_name, dst_name in band_map.items():
    out_df[dst_name] = df[src_name].values
write_csv(out_df, os.path.join(CSV_DIR, "planetscope_superdove.csv"))

# ── RapidEye ─────────────────────────────────────────────────────────────────
print("\nProcessing RapidEye_RSR.csv → rapideye.csv")
df = pd.read_csv(os.path.join(CSV_DIR, "RapidEye_RSR.csv"))
df.columns = [c.strip() for c in df.columns]
wl_col = df.columns[0]
wl_nm = df[wl_col].values * 1000.0  # µm → nm
# Rename "Red Edge" → "RedEdge"
cols = list(df.columns[1:])
cols = [c if c != "Red Edge" else "RedEdge" for c in cols]
bands_df = df[df.columns[1:]].copy()
bands_df.columns = cols
# Resample to 1 nm integer steps
resampled = resample_1nm(wl_nm, bands_df)
write_csv(resampled, os.path.join(CSV_DIR, "rapideye.csv"))

# ── SkySat ────────────────────────────────────────────────────────────────────

def load_skysat(filepath):
    """Load a single SkySat CSV (µm wavelength, BOM possible)."""
    df = pd.read_csv(filepath, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    wl_col = df.columns[0]
    wl_nm = df[wl_col].values * 1000.0  # µm → nm
    bands_df = df[df.columns[1:]].copy()
    bands_df.columns = ["Blue", "Green", "Red", "NIR", "Pan"]
    return wl_nm, bands_df


def average_skysats(file_list):
    """Average multiple SkySat files onto a common 1nm grid."""
    # Load all files
    all_data = []
    for fp in file_list:
        wl_nm, bands_df = load_skysat(fp)
        # Resample each to 1nm grid
        resampled = resample_1nm(wl_nm, bands_df)
        all_data.append(resampled)

    # Find common wavelength range (intersection)
    wl_min = max(d["WL_nm"].min() for d in all_data)
    wl_max = min(d["WL_nm"].max() for d in all_data)
    wl_common = np.arange(wl_min, wl_max + 1)

    # Filter each to common range and stack
    filtered = []
    for d in all_data:
        mask = (d["WL_nm"] >= wl_min) & (d["WL_nm"] <= wl_max)
        filtered.append(d[mask].reset_index(drop=True))

    # Average band columns
    band_cols = ["Blue", "Green", "Red", "NIR", "Pan"]
    avg = pd.DataFrame({"WL_nm": wl_common})
    for col in band_cols:
        stacked = np.stack([d[col].values for d in filtered], axis=0)
        avg[col] = np.mean(stacked, axis=0)
    return avg


# SkySat 1-2
print("\nProcessing Skysat 1+2 → skysat_1_2.csv")
files_1_2 = [
    os.path.join(CSV_DIR, "Skysat_RSR_Skysat1.csv"),
    os.path.join(CSV_DIR, "Skysat_RSR_Skysat2.csv"),
]
avg_1_2 = average_skysats(files_1_2)
write_csv(avg_1_2, os.path.join(CSV_DIR, "skysat_1_2.csv"))

# SkySat 3-13
print("\nProcessing Skysat 3-13 → skysat_3_13.csv")
files_3_13 = [
    os.path.join(CSV_DIR, f"Skysat_RSR_Skysat{i}.csv")
    for i in range(3, 14)
]
avg_3_13 = average_skysats(files_3_13)
write_csv(avg_3_13, os.path.join(CSV_DIR, "skysat_3_13.csv"))

# SkySat 14-19
print("\nProcessing Skysat 14-19 → skysat_14_19.csv")
fp_14_19 = os.path.join(CSV_DIR, "Skysat_RSR_Skysat14-SkySat19.csv")
wl_nm, bands_df = load_skysat(fp_14_19)
resampled = resample_1nm(wl_nm, bands_df)
write_csv(resampled, os.path.join(CSV_DIR, "skysat_14_19.csv"))

print("\nAll CSVs written successfully.")
