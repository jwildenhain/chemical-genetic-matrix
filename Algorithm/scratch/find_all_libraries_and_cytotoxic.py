import os
import glob
import pandas as pd

print("=================================================================")
print("SEARCHING FOR BIOACTIVE & CYTOTOXIC LIBRARIES IN CHEMGRID")
print("=================================================================")

chemgrid_dir = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/"
all_files = sorted(glob.glob(os.path.join(chemgrid_dir, "*")))

print("Files in ChemGRID data directory:")
for f in all_files:
    size_mb = os.path.getsize(f) / (1024 * 1024)
    print(f"  - {os.path.basename(f):45s} | Size: {size_mb:7.2f} MB")

# Check for cytotoxicity or bioactivity in CSV files
csv_files = glob.glob(os.path.join(chemgrid_dir, "*.csv"))
print("\nScanning CSV column headers for cytotoxicity or strain screening data:")
for cf in csv_files:
    try:
        df = pd.read_csv(cf, sep='\t', nrows=2)
        if df.shape[1] == 1:
            df = pd.read_csv(cf, nrows=2)
        cols = list(df.columns)
        cyto_cols = [c for c in cols if 'cyto' in c.lower() or 'bioact' in c.lower() or 'active' in c.lower() or 'screen' in c.lower()]
        print(f"\nFile: {os.path.basename(cf)}")
        print(f"  Total Columns: {len(cols)} | Total Rows Sample: {df.shape[0]}")
        if cyto_cols:
            print(f"  Relevant columns: {cyto_cols[:10]}")
    except Exception as e:
        print(f"  Error reading {os.path.basename(cf)}: {e}")
