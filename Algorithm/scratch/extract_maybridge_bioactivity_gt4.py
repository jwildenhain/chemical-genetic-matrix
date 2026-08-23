import os
import glob
import re
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("EXTRACTING MAYBRIDGE BIOACTIVE, CYTOTOXIC & 4 TARGET LIBRARIES (> 4 STRAINS)")
print("=================================================================")

# Load primary screening dataset
df_sdf = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
non_task_cols = {'smiles', 'compound_id', 'library'}
tasks = [c for c in df_sdf.columns if c not in non_task_cols]

# Filter df_sdf to compounds with > 4 strains screened
df_sdf['screened_count'] = df_sdf[tasks].notna().sum(axis=1)
df_primary_gt4 = df_sdf[df_sdf['screened_count'] > 4].copy()
print(f"Primary Spectrum Dataset (>4 strains screened): {len(df_primary_gt4)} / {len(df_sdf)} compounds")

# Check Maybridge & target library compounds in ChemGRID2011_sdf_chemgrid.csv
cg2011_csv = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/ChemGRID2011_sdf_chemgrid.csv"
maybridge_gt4_list = []

if os.path.exists(cg2011_csv):
    print("\nParsing ChemGRID2011_sdf_chemgrid.csv for Maybridge Bioactive & Cytotoxic screening data...")
    # Load Maybridge rows or rows with db containing Maybridge
    df_cg = pd.read_csv(cg2011_csv, sep='\t', low_memory=False)
    print("ChemGRID2011 Total Rows:", len(df_cg))
    
    # Filter rows with Maybridge identifiers or Maybridge library
    may_mask = (df_cg['db'].astype(str).str.contains('Maybridge', case=False, na=False)) | \
               (df_cg['code'].astype(str).str.startswith('MAY')) | \
               (df_cg['ID'].astype(str).str.startswith('MAY'))
    df_may = df_cg[may_mask].copy()
    print(f"Found {len(df_may)} Maybridge library compound entries.")
    
    # Match strain columns in ChemGRID2011
    cg_strain_cols = [c for c in df_cg.columns if c in tasks]
    print(f"Matched yeast strain columns in ChemGRID2011: {len(cg_strain_cols)}")
    
    if cg_strain_cols:
        df_may['screened_count'] = df_may[cg_strain_cols].notna().sum(axis=1)
        df_may_gt4 = df_may[df_may['screened_count'] > 4].copy()
        print(f"Maybridge Compounds with > 4 strains screened in ChemGRID2011: {len(df_may_gt4)}")
        
        # Load SDF structures for Maybridge compounds
        may_sdfs = [
            '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
            '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_MaybridgeHitskit.sdf',
            '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_ChemGRID2011.sdf'
        ]
        
        smiles_dict = {}
        for msf in may_sdfs:
            if os.path.exists(msf) and os.path.getsize(msf) > 0:
                print(f"Featurizing structures from {os.path.basename(msf)}...")
                suppl = Chem.SDMolSupplier(msf, sanitize=False)
                for m in suppl:
                    if m is None:
                        continue
                    try:
                        Chem.SanitizeMol(m)
                        smiles = Chem.MolToSmiles(m)
                        if smiles:
                            for p in ['code', 'cid', 'ID', 'id', 'CATNUM']:
                                if m.HasProp(p):
                                    cval = m.GetProp(p).strip()
                                    if cval:
                                        smiles_dict[cval] = smiles
                    except Exception:
                        continue
                        
        print(f"Loaded {len(smiles_dict)} Maybridge SMILES mappings.")
        
        # Attach SMILES to Maybridge gt4 entries
        for _, row in df_may_gt4.iterrows():
            cid = str(row.get('code', row.get('ID', '')))
            smiles = smiles_dict.get(cid)
            if not smiles and 'parent_smiles' in row and pd.notna(row['parent_smiles']):
                smiles = str(row['parent_smiles'])
                
            if smiles:
                item = {
                    'compound_id': cid,
                    'library': 'Maybridge_Bioactive_Cytotoxic',
                    'smiles': smiles
                }
                for t in tasks:
                    if t in row and pd.notna(row[t]):
                        item[t] = row[t]
                    else:
                        item[t] = np.nan
                maybridge_gt4_list.append(item)

# Combine primary gt4 dataset and Maybridge gt4 dataset
df_primary_clean = df_primary_gt4[['compound_id', 'library', 'smiles'] + tasks].copy()

if maybridge_gt4_list:
    df_may_clean = pd.DataFrame(maybridge_gt4_list)
    print(f"\nAdding {len(df_may_clean)} Maybridge Bioactive/Cytotoxic screened compounds...")
    df_final = pd.concat([df_primary_clean, df_may_clean], ignore_index=True)
else:
    df_final = df_primary_clean

df_final = df_final.drop_duplicates(subset=['smiles'])

output_file = "data/yeast_bioactivity_target_libraries_gt4.csv"
df_final.to_csv(output_file, index=False)

print(f"\n=================================================================")
print(f"FINAL TARGET DATASET (> 4 STRAINS SCREENED): {len(df_final)} COMPOUNDS")
print("=================================================================")
print(df_final['library'].value_counts().to_string())
