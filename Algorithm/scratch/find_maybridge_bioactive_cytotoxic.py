import os
import glob
import gzip
import pandas as pd

print("=================================================================")
print("SEARCHING FOR MAYBRIDGE BIOACTIVE & CYTOTOXIC LIBRARIES (>4 STRAINS)")
print("=================================================================")

chemgrid_dir = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/"
csv_files = glob.glob(os.path.join(chemgrid_dir, "*.csv"))

for cf in sorted(csv_files):
    try:
        df = pd.read_csv(cf, sep='\t', nrows=5)
        if df.shape[1] == 1:
            df = pd.read_csv(cf, nrows=5)
        cols = list(df.columns)
        
        # Check strain columns (columns that aren't metadata)
        non_strain = {'sdf_structure', 'ID', 'mol_formular', 'mol_mass', 'product_name', 'db', 'c_log_p', 
                      'PM_tpsa', 'parent_smiles', 'code', 'cid', 'xsd', 'xmean', 'xmedian', 'xmin', 'xmax', 
                      'hg', 'active', 'resistent', 'nonactive', 'screened', 'ratio'}
        strains = [c for c in cols if c not in non_strain and not c.startswith('BCUT') and not c.startswith('PPSA')]
        
        # Check Maybridge identifiers or bioactive/cytotoxic names
        has_maybridge = any('may' in c.lower() for c in cols) or 'may' in os.path.basename(cf).lower()
        has_cyto_bio = any('cyto' in c.lower() or 'bioact' in c.lower() for c in cols) or 'cyto' in os.path.basename(cf).lower() or 'bio' in os.path.basename(cf).lower()
        
        print(f"\nFile: {os.path.basename(cf)}")
        print(f"  - Total Columns: {len(cols)} | Strains Count: {len(strains)}")
        print(f"  - Maybridge Relevant: {has_maybridge} | Cyto/Bioactive Relevant: {has_cyto_bio}")
        if strains:
            print(f"  - Strain Columns Sample: {strains[:10]}")
    except Exception as e:
        print(f"Error reading {os.path.basename(cf)}: {e}")

# Check ChemGrid.latest.sql and daily dumps for Maybridge tables
sql_files = [
    '/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql'
]

for sf in sql_files:
    if os.path.exists(sf):
        print(f"\nScanning SQL File: {os.path.basename(sf)}...")
        with open(sf, 'r', errors='ignore') as f:
            may_tables = set()
            cyto_tables = set()
            for line in f:
                if 'CREATE TABLE' in line:
                    tname = line.split()[2].strip('`(')
                    if 'may' in tname.lower() or 'bio' in tname.lower() or 'cyto' in tname.lower():
                        may_tables.add(tname)
            print("  - Relevant Tables Found:", sorted(list(may_tables)))
