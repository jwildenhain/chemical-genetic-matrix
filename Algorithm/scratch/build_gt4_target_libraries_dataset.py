import os
import glob
import re
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("BUILDING DATASET FOR TARGET LIBRARIES WITH > 4 STRAINS SCREENED")
print("=================================================================")

# Load compound structures for target libraries (Spectrum, SPECMTS3, Maybridge1000, Lopac)
sdfs = {
    'Spectrum': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTRUM_ED.sdf',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_SPECMTS3.sdf',
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Lopac.sdf'
}

all_records = []
seen_smiles = set()

def extract_cid(m, lib_name):
    for p in ['cid', 'code', 'id', 'ID', 'CATNUM', 'product_name']:
        if m.HasProp(p):
            val = m.GetProp(p).strip()
            if val:
                num = re.sub(r'[^0-9]', '', val)
                if num:
                    if 'LOP' in val or lib_name == 'Lopac':
                        return 'LOP' + num.zfill(8)
                    elif 'SPE' in val or lib_name in ('Spectrum', 'SPECMTS3'):
                        return 'SPE' + num.zfill(8)
                    elif 'MAY' in val or lib_name == 'Maybridge1000':
                        return 'MAY' + num.zfill(8)
                    return val
    return None

for lib_name, sf in sdfs.items():
    if os.path.exists(sf):
        suppl = Chem.SDMolSupplier(sf, sanitize=False)
        count = 0
        for m in suppl:
            if m is None:
                continue
            try:
                Chem.SanitizeMol(m)
                smiles = Chem.MolToSmiles(m)
                if smiles and smiles not in seen_smiles:
                    seen_smiles.add(smiles)
                    cid = extract_cid(m, lib_name)
                    if not cid:
                        cid = f"{lib_name}_{count+1}"
                    all_records.append({
                        'compound_id': cid,
                        'library': lib_name,
                        'smiles': smiles
                    })
                    count += 1
            except Exception:
                continue
        print(f"Library {lib_name:15s}: Loaded {count} structures.")

df_compounds = pd.DataFrame(all_records)
print(f"\nTotal Target Library Compounds: {len(df_compounds)}")

# Merge with screening matrix across 102 yeast strain tasks
df_screen = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_screen.columns if c not in {'smiles', 'compound_id', 'library'}]

df_merged = pd.merge(df_compounds, df_screen[['smiles'] + tasks], on='smiles', how='left')

# Calculate number of screened yeast strain tasks per compound (> 4 filter)
df_merged['screened_count'] = df_merged[tasks].notna().sum(axis=1)
df_gt4 = df_merged[df_merged['screened_count'] > 4].copy()

print(f"\nCompounds with > 4 Yeast Strains Screened: {len(df_gt4)} / {len(df_merged)} ({len(df_gt4)/len(df_merged)*100:.1f}%)")
print("Library breakdown of > 4 strains screened:\n", df_gt4['library'].value_counts().to_string())

output_path = "data/yeast_bioactivity_target_libraries_gt4.csv"
df_gt4.to_csv(output_path, index=False)
print(f"\nSaved filtered dataset to {output_path}: {len(df_gt4)} compounds x {len(tasks)} tasks.")
