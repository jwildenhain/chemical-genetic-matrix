import os
import sys
import re
import math
import numpy as np
import pandas as pd
from collections import defaultdict
from rdkit import Chem

# Ensure output directories exist
os.makedirs("data", exist_ok=True)
os.makedirs("reports", exist_ok=True)

sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
sql_path = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

target_libs = {"Spectrum", "SPECMTS3", "Maybridge1000", "Lopac"}

print("=================================================================")
print("STEP 1: Loading Compound SMILES from SDF")
print("=================================================================")

cmpd_smiles = {}
cmpd_names = {}
cmpd_libs = {}

supplier = Chem.SDMolSupplier(sdf_path)
for i, mol in enumerate(supplier):
    if mol is None:
        continue
    code = mol.GetProp('code') if mol.HasProp('code') else (mol.GetProp('id') if mol.HasProp('id') else None)
    smiles = mol.GetProp('parent_smiles') if mol.HasProp('parent_smiles') else Chem.MolToSmiles(mol)
    pname = mol.GetProp('product_name') if mol.HasProp('product_name') else ''
    db = mol.GetProp('db') if mol.HasProp('db') else ''
    
    if code and smiles:
        c_code = code.strip()
        cmpd_smiles[c_code] = smiles
        cmpd_names[c_code] = pname
        cmpd_libs[c_code] = db

print(f"Loaded {len(cmpd_smiles)} valid compound structures from SDF.")

print("\n=================================================================")
print("STEP 2: Parsing Screening Data & Resolving Strain Time-Points")
print("=================================================================")

# Raw parsed measurements: strain -> compound_id -> list of dicts
raw_screen_data = defaultdict(lambda: defaultdict(list))
all_parsed_strains = set()

with open(sql_path, 'r', encoding='latin1') as f:
    for line in f:
        if line.startswith("INSERT INTO `exp_info`") or line.startswith("INSERT INTO `exp_raw`"):
            idx = line.find("VALUES ")
            if idx != -1:
                content = line[idx + 7:]
                items = content.split("),(")
                for item in items:
                    item = item.lstrip("(").rstrip(");\n")
                    parts = item.split(",")
                    if len(parts) >= 12:
                        row = parts[0].strip("' ")
                        col = parts[1].strip("' ")
                        val_str = parts[2].strip("' ")
                        expid = parts[6].strip("' ") # Strain / ORF ID
                        lib = parts[11].strip("' ")
                        
                        if lib in target_libs and expid and val_str:
                            try:
                                val = float(val_str)
                                all_parsed_strains.add(expid)
                                raw_screen_data[expid][lib].append({
                                    'value': val,
                                    'row': row,
                                    'col': col,
                                    'library': lib
                                })
                            except ValueError:
                                pass

print(f"Parsed screening measurements across {len(all_parsed_strains)} unique raw strain strings.")

# Time-point resolution logic: group strains by base name and pick earliest time point
# e.g., Alg12-18h, Alg12-20h, Alg12-22h -> Alg12-18h
strain_timepoint_map = {}
base_strain_groups = defaultdict(list)

time_pattern = re.compile(r"^(.*?)[_-]?(\d+)\s*h?$", re.IGNORECASE)

for strain in all_parsed_strains:
    match = time_pattern.match(strain)
    if match:
        base_name = match.group(1).rstrip("_-")
        tp = int(match.group(2))
        base_strain_groups[base_name].append((tp, strain))
    else:
        base_strain_groups[strain].append((999, strain))

selected_strains = set()
strain_rename_map = {}

for base_name, tp_list in base_strain_groups.items():
    tp_list.sort(key=lambda x: x[0]) # sort by time point ascending (earliest first)
    earliest_tp, earliest_strain = tp_list[0]
    selected_strains.add(earliest_strain)
    strain_rename_map[earliest_strain] = base_name if earliest_tp != 999 else earliest_strain

print(f"Resolved strain time-points: Selected {len(selected_strains)} earliest-time-point strains out of {len(all_parsed_strains)} total.")

print("\n=================================================================")
print("STEP 3: Detecting & Flagging Extreme Outliers")
print("=================================================================")

outlier_records = []
clean_strain_cmpd_vals = defaultdict(lambda: defaultdict(list))

for strain, lib_dict in raw_screen_data.items():
    if strain not in selected_strains:
        continue
    final_strain_name = strain_rename_map[strain]
    
    for lib, recs in lib_dict.items():
        vals = [r['value'] for r in recs]
        if len(vals) == 0:
            continue
        
        mean_v = np.mean(vals)
        std_v = np.std(vals) if len(vals) > 1 else 0.0
        
        for r in recs:
            v = r['value']
            # Flag extreme outliers (value < 0 or > 3.5 or |Z| > 5.0 relative to strain mean)
            is_outlier = False
            z_val = ((v - mean_v) / std_v) if std_v > 1e-5 else 0.0
            
            if v < -0.2 or v > 4.0 or abs(z_val) > 5.0:
                is_outlier = True
                outlier_records.append({
                    'strain_raw': strain,
                    'strain_cleaned': final_strain_name,
                    'library': lib,
                    'row': r['row'],
                    'col': r['col'],
                    'value': v,
                    'strain_mean': round(mean_v, 4),
                    'strain_std': round(std_v, 4),
                    'z_score': round(z_val, 2),
                    'reason': 'Extreme reading' if (v < -0.2 or v > 4.0) else '|Z| > 5.0'
                })
            
            if not is_outlier:
                clean_strain_cmpd_vals[final_strain_name][lib].append(v)

outlier_df = pd.DataFrame(outlier_records)
outlier_df.to_csv("reports/extreme_outliers.csv", index=False)
print(f"Flagged {len(outlier_records)} extreme outlier measurements. Saved to reports/extreme_outliers.csv.")

print("\n=================================================================")
print("STEP 4: Aggregating Replicates via Median Calculation & Building Matrix")
print("=================================================================")

plate_well_cmpd = {}

with open(sql_path, 'r', encoding='latin1') as f:
    for line in f:
        if line.startswith("INSERT INTO `info_cmpd2plate`"):
            idx = line.find("VALUES ")
            if idx != -1:
                content = line[idx + 7:]
                items = content.split("),(")
                for item in items:
                    item = item.lstrip("(").rstrip(");\n")
                    parts = item.split(",")
                    if len(parts) >= 29:
                        obj_id = parts[0].strip("' ")
                        supplier_obj_id = parts[12].strip("' ") if len(parts) > 12 else parts[0].strip("' ")
                        plate_num = parts[25].strip("' ") if len(parts) > 25 else ""
                        p_row = parts[26].strip("' ") if len(parts) > 26 else ""
                        p_col = parts[27].strip("' ") if len(parts) > 27 else ""
                        db = parts[28].strip("' ") if len(parts) > 28 else ""
                        
                        cid = supplier_obj_id if supplier_obj_id in cmpd_smiles else (obj_id if obj_id in cmpd_smiles else None)
                        if cid and db in target_libs and plate_num and p_row and p_col:
                            key = (db, plate_num, p_row, p_col)
                            plate_well_cmpd[key] = cid

print(f"Mapped {len(plate_well_cmpd)} plate/well coordinates to compound SMILES.")

# Build strain -> compound -> list of measurements
strain_cmpd_reps = defaultdict(lambda: defaultdict(list))

with open(sql_path, 'r', encoding='latin1') as f:
    for line in f:
        if line.startswith("INSERT INTO `exp_info`") or line.startswith("INSERT INTO `exp_raw`"):
            idx = line.find("VALUES ")
            if idx != -1:
                content = line[idx + 7:]
                items = content.split("),(")
                for item in items:
                    item = item.lstrip("(").rstrip(");\n")
                    parts = item.split(",")
                    if len(parts) >= 12:
                        p_row = parts[0].strip("' ")
                        p_col = parts[1].strip("' ")
                        val_str = parts[2].strip("' ")
                        expid = parts[6].strip("' ")
                        plate_num = parts[9].strip("' ") if len(parts) > 9 else ""
                        lib = parts[11].strip("' ")
                        
                        if lib in target_libs and expid in selected_strains:
                            final_strain = strain_rename_map[expid]
                            key = (lib, plate_num, p_row, p_col)
                            cid = plate_well_cmpd.get(key)
                            if cid:
                                try:
                                    val = float(val_str)
                                    strain_cmpd_reps[final_strain][cid].append(val)
                                except ValueError:
                                    pass

print(f"Successfully aligned compound-strain measurements across {len(strain_cmpd_reps)} strains!")

# Compute Medians and Binarize:
# Active (1) if median growth < 0.75 (inhibitory/sensitive hit)
# Inactive (0) if median growth >= 0.85
all_compounds_set = set()
valid_strains = []

for strain, cmpd_dict in strain_cmpd_reps.items():
    if len(cmpd_dict) >= 20: # Keep strains with at least 20 compound screens
        valid_strains.append(strain)
        for cid in cmpd_dict.keys():
            all_compounds_set.add(cid)

valid_strains.sort()
sorted_compounds = sorted(all_compounds_set)

print(f"Retained {len(valid_strains)} valid strains and {len(sorted_compounds)} unique compounds for multi-task dataset.")

# Construct DataFrame
dataset_rows = []
for cid in sorted_compounds:
    smiles = cmpd_smiles[cid]
    lib = cmpd_libs.get(cid, '')
    row_dict = {
        'smiles': smiles,
        'compound_id': cid,
        'library': lib
    }
    for strain in valid_strains:
        vals = strain_cmpd_reps[strain].get(cid, [])
        if len(vals) > 0:
            med_val = float(np.median(vals))
            if med_val < 0.75:
                row_dict[strain] = 1
            elif med_val >= 0.85:
                row_dict[strain] = 0
            else:
                row_dict[strain] = np.nan
        else:
            row_dict[strain] = np.nan
    dataset_rows.append(row_dict)

dataset_df = pd.DataFrame(dataset_rows)
csv_output_path = "data/yeast_bioactivity_multitask.csv"
dataset_df.to_csv(csv_output_path, index=False)

print("\n=================================================================")
print("EXTRACTION SUMMARY & DATASET STATS")
print("=================================================================")
print(f"Saved dataset: {csv_output_path}")
print(f"Total Compounds: {len(dataset_df)}")
print(f"Total Yeast Strain Tasks: {len(valid_strains)}")
print(f"Target Strains: {valid_strains[:10]}...")
