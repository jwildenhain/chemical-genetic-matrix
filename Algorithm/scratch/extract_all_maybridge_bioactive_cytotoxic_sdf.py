import os
import re
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("EXTRACTING MAYBRIDGE BIOACTIVE & CYTOTOXIC + TARGET LIBRARIES DATASET (>4 STRAINS)")
print("=================================================================")

# Load Maybridge1000 dataset
may_csv = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Maybridge1000_sdf_chemgrid.csv"
may_sdf = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf"

df_may = pd.read_csv(may_csv, sep='\t', on_bad_lines='skip', low_memory=False)
df_may['screened'] = pd.to_numeric(df_may['screened'], errors='coerce')
df_may_gt4 = df_may[df_may['screened'] > 4].copy()

print(f"Maybridge Compounds with > 4 strains screened: {len(df_may_gt4)}")

# Map Maybridge compound codes to SMILES from SDF
may_smiles = {}
if os.path.exists(may_sdf):
    suppl = Chem.SDMolSupplier(may_sdf, sanitize=False)
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
                            may_smiles[cval] = smiles
        except Exception:
            continue

print(f"Mapped {len(may_smiles)} Maybridge SMILES structures.")

# Load target screening dataset
df_target = pd.read_csv("data/yeast_bioactivity_target_libraries_gt4.csv")
tasks = [c for c in df_target.columns if c not in {'smiles', 'compound_id', 'library', 'screened_count'}]

# Create records for Maybridge gt4 compounds
may_records = []
for _, row in df_may_gt4.iterrows():
    cid = str(row.get('code', row.get('ID', '')))
    smiles = may_smiles.get(cid)
    if not smiles and 'parent_smiles' in row and pd.notna(row['parent_smiles']):
        smiles = str(row['parent_smiles'])
        
    if smiles:
        item = {
            'compound_id': cid if cid else f"MAYBRIDGE_{len(may_records)+1}",
            'library': 'Maybridge_Bioactive_Cytotoxic',
            'smiles': smiles
        }
        for t in tasks:
            if t in row and pd.notna(row[t]):
                item[t] = row[t]
            else:
                item[t] = np.nan
        may_records.append(item)

df_may_records = pd.DataFrame(may_records)
print(f"Created {len(df_may_records)} Maybridge Bioactive & Cytotoxic screened compound records.")

# Combine with existing target library compounds
df_combined = pd.concat([df_target[['compound_id', 'library', 'smiles'] + tasks], df_may_records], ignore_index=True)
df_combined = df_combined.drop_duplicates(subset=['smiles'])

output_file = "data/yeast_bioactivity_target_libraries_gt4.csv"
df_combined.to_csv(output_file, index=False)

print(f"\n=================================================================")
print(f"FINAL COMBINED TARGET DATASET (> 4 STRAINS SCREENED): {len(df_combined)} COMPOUNDS")
print("=================================================================")
print(df_combined['library'].value_counts().to_string())
