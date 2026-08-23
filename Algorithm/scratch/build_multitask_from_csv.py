import os
import glob
import re
import numpy as np
import pandas as pd
from collections import defaultdict
from rdkit import Chem

print("=================================================================")
print("BUILDING CLEAN MULTI-TASK YEAST BIOACTIVITY DATASET")
print("=================================================================")

target_libs = {"Spectrum", "SPECMTS3", "Maybridge1000", "Lopac"}

# Load 100% valid SMILES structures from cgm_ChemGRID2011.sdf
sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
cmpd_smiles = {}
cmpd_libs = {}

supplier = Chem.SDMolSupplier(sdf_path)
for mol in supplier:
    if mol is None:
        continue
    code = mol.GetProp('code') if mol.HasProp('code') else (mol.GetProp('id') if mol.HasProp('id') else None)
    smiles = mol.GetProp('parent_smiles') if mol.HasProp('parent_smiles') else Chem.MolToSmiles(mol)
    db = mol.GetProp('db') if mol.HasProp('db') else ''
    
    if code and smiles:
        clean_c = code.strip()
        # Verify valid RDKit SMILES
        if Chem.MolFromSmiles(smiles) is not None:
            cmpd_smiles[clean_c] = smiles
            cmpd_libs[clean_c] = db

print(f"Loaded {len(cmpd_smiles)} verified compound SMILES from SDF.")

# Search for pre-joined CSV files in zfspool
csv_files = [
    "/mnt/zfspool/Backups/Jan/2019_ugreen/usbsticks/usbstick_kingston_pwprotected/CGM/csv/chemGRID_Z_norm_Spectrum_4_PP.csv",
    "/mnt/zfspool/Backups/Jan/2019_ugreen/usbsticks/usbstick_kingston_pwprotected/CGM/csv/stats_properties_chemGRID_Spectrum_4932.csv"
]

dfs = []
for f in csv_files:
    if os.path.exists(f):
        print(f"Reading CSV: {f}")
        try:
            df = pd.read_csv(f, low_memory=False)
            dfs.append(df)
        except Exception as e:
            print("Error reading CSV:", e)

combined_raw_df = pd.concat(dfs, ignore_index=True)
print(f"Combined CSV Rows: {len(combined_raw_df)}")

# Process measurements
measurements = []

for idx, r in combined_raw_df.iterrows():
    cid = str(r.get('supplier_obj_id', '')).strip()
    smiles_raw = str(r.get('parent_smiles', '')).strip()
    strain = str(r.get('strain', '')).strip()
    lib = str(r.get('library', '')).strip()
    z = r.get('z_score', np.nan)
    val = r.get('value', np.nan)
    
    # Priority: Lookup in verified cmpd_smiles SDF dictionary first!
    valid_smiles = cmpd_smiles.get(cid)
    if not valid_smiles and Chem.MolFromSmiles(smiles_raw) is not None:
        valid_smiles = smiles_raw
        
    if valid_smiles and strain and strain != 'nan' and (lib in target_libs or not lib):
        try:
            z_f = float(z)
        except (ValueError, TypeError):
            z_f = np.nan
        try:
            val_f = float(val)
        except (ValueError, TypeError):
            val_f = np.nan
            
        measurements.append({
            'compound_id': cid,
            'smiles': valid_smiles,
            'strain': strain,
            'library': lib,
            'z_score': z_f,
            'value': val_f
        })

meas_df = pd.DataFrame(measurements)
print(f"Extracted {len(meas_df)} measurement records with verified SMILES.")

# Parse time-points & resolve earliest time-point per strain
# e.g., Alg12-18h vs Alg12-20h vs Alg12-22h
all_strains = meas_df['strain'].unique()
base_strain_groups = defaultdict(list)
time_pattern = re.compile(r"^(.*?)[_-]?(\d+)\s*h?$", re.IGNORECASE)

for s in all_strains:
    match = time_pattern.match(s)
    if match:
        base_name = match.group(1).rstrip("_-")
        tp = int(match.group(2))
        base_strain_groups[base_name].append((tp, s))
    else:
        base_strain_groups[s].append((999, s))

selected_strains = set()
strain_rename_map = {}

for base_name, tp_list in base_strain_groups.items():
    tp_list.sort(key=lambda x: x[0])
    earliest_tp, earliest_strain = tp_list[0]
    selected_strains.add(earliest_strain)
    strain_rename_map[earliest_strain] = base_name if earliest_tp != 999 else earliest_strain

print(f"Selected {len(selected_strains)} earliest time-point strains out of {len(all_strains)} total.")

meas_df = meas_df[meas_df['strain'].isin(selected_strains)].copy()
meas_df['clean_strain'] = meas_df['strain'].map(strain_rename_map)

# Flag extreme outliers
outliers = []
clean_rows = []

for strain, group in meas_df.groupby('clean_strain'):
    mean_z = group['z_score'].mean()
    std_z = group['z_score'].std()
    
    for _, row in group.iterrows():
        z_v = row['z_score']
        v_v = row['value']
        
        is_outlier = False
        if not np.isnan(z_v) and abs(z_v) > 5.0:
            is_outlier = True
        elif not np.isnan(v_v) and (v_v < -0.5 or v_v > 3.5):
            is_outlier = True
            
        if is_outlier:
            outliers.append(row.to_dict())
        else:
            clean_rows.append(row.to_dict())

outliers_df = pd.DataFrame(outliers)
outliers_df.to_csv("reports/extreme_outliers.csv", index=False)
print(f"Flagged {len(outliers_df)} extreme outliers -> Saved to reports/extreme_outliers.csv.")

clean_df = pd.DataFrame(clean_rows)

# Median aggregation across replicates for each (compound_id, clean_strain)
aggregated = clean_df.groupby(['compound_id', 'smiles', 'library', 'clean_strain']).agg({
    'z_score': 'median',
    'value': 'median'
}).reset_index()

print(f"Aggregated replicate measurements: {len(aggregated)} unique compound-strain pairs.")

# Create Pivot Table: rows = SMILES / compound_id, columns = clean_strain
pivot_z = aggregated.pivot_table(index=['smiles', 'compound_id', 'library'], columns='clean_strain', values='z_score')
pivot_val = aggregated.pivot_table(index=['smiles', 'compound_id', 'library'], columns='clean_strain', values='value')

# Binarization rule:
# Active (1) if z_score <= -2.5 or value < 0.75 (sensitive/inhibitory hit)
# Inactive (0) if |z_score| < 1.0 or value >= 0.85
binary_matrix = pd.DataFrame(index=pivot_z.index, columns=pivot_z.columns)

for col in pivot_z.columns:
    z_col = pivot_z[col]
    val_col = pivot_val[col]
    
    # Active: z <= -2.5 or val < 0.75
    active_mask = (z_col <= -2.5) | (val_col < 0.75)
    # Inactive: abs(z) < 1.0 or val >= 0.85
    inactive_mask = (z_col.abs() < 1.0) | (val_col >= 0.85)
    
    binary_matrix[col] = np.nan
    binary_matrix.loc[inactive_mask, col] = 0
    binary_matrix.loc[active_mask, col] = 1

# Filter strains with at least 15 total screens
min_screens = 15
valid_strain_cols = [c for c in binary_matrix.columns if binary_matrix[c].notna().sum() >= min_screens]
binary_matrix = binary_matrix[valid_strain_cols].reset_index()

output_csv = "data/yeast_bioactivity_multitask.csv"
binary_matrix.to_csv(output_csv, index=False)

print("\n=================================================================")
print("FINAL DATASET CREATION SUMMARY")
print("=================================================================")
print(f"Saved dataset: {output_csv}")
print(f"Total Compounds: {len(binary_matrix)}")
print(f"Total Yeast Strain Tasks: {len(valid_strain_cols)}")
print(f"Sample Strain Tasks: {valid_strain_cols[:5]}")
