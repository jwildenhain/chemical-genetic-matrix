import os
import glob
import re
import pandas as pd
from rdkit import Chem

print("=================================================================")
print("EXTENSIVE SEARCH FOR ALL MOLECULAR IDENTIFIERS & STRUCTURES")
print("=================================================================")

chemgrid_dir = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/"
all_sdfs = glob.glob(os.path.join(chemgrid_dir, "*.sdf"))

print(f"Found {len(all_sdfs)} SDF structure files in {chemgrid_dir}:")
for sf in sorted(all_sdfs):
    size_mb = os.path.getsize(sf) / (1024 * 1024)
    print(f"  - {os.path.basename(sf):35s} | {size_mb:7.2f} MB")

master_records = []
seen_smiles = set()

def extract_cid(m, lib_name, fallback_id):
    for p in ['cid', 'code', 'id', 'ID', 'CATNUM', 'product_name', 'CAS']:
        if m.HasProp(p):
            val = m.GetProp(p).strip()
            if val:
                return val
    return fallback_id

for sf in sorted(all_sdfs):
    size_b = os.path.getsize(sf)
    if size_b == 0:
        print(f"\nSkipping 0-byte file: {os.path.basename(sf)}")
        continue
        
    lib_name = os.path.basename(sf).replace('.sdf', '')
    print(f"\nProcessing {os.path.basename(sf)}...")
    suppl = Chem.SDMolSupplier(sf, sanitize=False)
    added = 0
    for idx, m in enumerate(suppl):
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
            smiles = Chem.MolToSmiles(m)
            if smiles and smiles not in seen_smiles:
                seen_smiles.add(smiles)
                cid = extract_cid(m, lib_name, f"{lib_name}_{idx+1}")
                master_records.append({
                    'compound_id': cid,
                    'library': lib_name,
                    'smiles': smiles
                })
                added += 1
        except Exception:
            continue
    print(f"  -> Added {added} new unique molecular structures.")

df_compounds = pd.DataFrame(master_records)
print(f"\n=================================================================")
print(f"TOTAL UNIQUE MOLECULAR STRUCTURES EXTRACTED: {len(df_compounds)}")
print("=================================================================")
print(df_compounds['library'].value_counts().to_string())

# Merge with screening matrix across 102 yeast strain tasks
df_screen = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_screen.columns if c not in {'smiles', 'compound_id', 'library'}]

df_extended = pd.merge(df_compounds, df_screen[['smiles'] + tasks], on='smiles', how='left')

output_path = "data/yeast_bioactivity_extended_master.csv"
df_extended.to_csv(output_path, index=False)
print(f"\nSaved extended dataset to {output_path}: {len(df_extended)} compounds x {len(tasks)} tasks.")
