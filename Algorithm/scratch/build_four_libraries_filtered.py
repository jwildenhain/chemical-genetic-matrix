import os
import gzip
import pandas as pd
import numpy as np

print("=================================================================")
print("CHECKING SCREENED STRAIN COUNTS PER COMPOUND (> 4 STRAINS SCREENED)")
print("=================================================================")

files = {
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Maybridge1000_sdf_chemgrid.csv',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Lopac_sdf_chemgrid.csv',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTMTS3_sdf_chemgrid.csv',
    'Spectrum128': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Spectrum128_sdf_chemgrid.csv'
}

df_sdf_all = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_sdf_all.columns if c not in {'smiles', 'compound_id', 'library'}]

print(f"Total tasks evaluated: {len(tasks)}")

filtered_records = []

for lib_name, filepath in files.items():
    if os.path.exists(filepath):
        df = pd.read_csv(filepath, sep='\t')
        cid_col = 'code' if 'code' in df.columns else ('ID' if 'ID' in df.columns else 'cid')
        smiles_col = 'parent_smiles' if 'parent_smiles' in df.columns else 'smiles'
        
        # Match strain columns in this CSV
        strain_cols = [c for c in df.columns if c in tasks]
        
        print(f"\nLibrary: {lib_name:15s} | File: {os.path.basename(filepath)}")
        print(f"  - Total Rows: {len(df)} | Matched Strain Columns: {len(strain_cols)}")
        
        # Calculate number of screened strains per compound (non-null measurements)
        if strain_cols:
            df['screened_count'] = df[strain_cols].notna().sum(axis=1)
            valid_gt4 = df[df['screened_count'] > 4]
            print(f"  - Compounds with > 4 Strains Screened: {len(valid_gt4)} / {len(df)} ({len(valid_gt4)/len(df)*100:.1f}%)")
            
            # Show summary of screened counts
            print(f"  - Screened Strains Distribution: Min={df['screened_count'].min()}, Median={df['screened_count'].median()}, Max={df['screened_count'].max()}")

# Check SQL dumps for Maybridge Bioactive & Cytotoxic tables
dump_path = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/daily_all-databases_2014-08-04.sql.gz"
if os.path.exists(dump_path):
    print(f"\nScanning {os.path.basename(dump_path)} for Maybridge Bioactive & Cytotoxic tables...")
    with gzip.open(dump_path, 'rt', errors='ignore') as f:
        tables = set()
        for line in f:
            if 'CREATE TABLE' in line:
                tname = line.split()[2].strip('`(')
                if any(k in tname.lower() for k in ['may', 'bioact', 'cyto', 'toxic']):
                    tables.add(tname)
        print("  - Found Tables:", sorted(list(tables)))
