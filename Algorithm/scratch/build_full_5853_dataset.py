import os
import re
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("BUILDING FULL 5,853-COMPOUND DATASET ACROSS 4 TARGET LIBRARIES")
print("=================================================================")

# Load compound structures from library SDFs
sdfs = {
    'Spectrum': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTRUM_ED.sdf',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_SPECMTS3.sdf',
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Lopac.sdf'
}

all_records = []

def extract_cid(m, lib_name):
    for p in ['cid', 'code', 'id', 'ID', 'CATNUM']:
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
        valid = 0
        for m in suppl:
            if m is None:
                continue
            try:
                Chem.SanitizeMol(m)
                smiles = Chem.MolToSmiles(m)
                cid = extract_cid(m, lib_name)
                if smiles and cid:
                    all_records.append({
                        'compound_id': cid,
                        'library': lib_name,
                        'smiles': smiles
                    })
                    valid += 1
            except Exception:
                continue
        print(f"Library {lib_name:15s}: Loaded {valid} valid structures.")

df_all = pd.DataFrame(all_records).drop_duplicates(subset=['compound_id'])
print(f"\nTotal Unique Library Compounds: {len(df_all)}")

# Merge experimental screening data from yeast_bioactivity_multitask_sdf.csv
df_existing = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_existing.columns if c not in {'smiles', 'compound_id', 'library'}]

df_merged = pd.merge(df_all, df_existing[['compound_id'] + tasks], on='compound_id', how='left')

output_path = "data/yeast_bioactivity_all_libraries.csv"
df_merged.to_csv(output_path, index=False)
print(f"Saved complete 4-library dataset ({len(df_merged)} compounds x {len(tasks)} tasks) to {output_path}.")
