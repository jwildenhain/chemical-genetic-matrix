import os
import re
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("EXTRACTING ALL LIBRARIES & STRAIN BIOACTIVITY/CYTOTOXIC DATA")
print("=================================================================")

# Load compound structures from all available SDF files
sdf_files = {
    'Spectrum': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTRUM_ED.sdf',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_SPECMTS3.sdf',
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Lopac.sdf',
    'ChemGRID2011': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_ChemGRID2011.sdf'
}

records = []
seen_smiles = set()

for lib_name, sf in sdf_files.items():
    if not os.path.exists(sf):
        continue
    print(f"Reading SDF file: {os.path.basename(sf)}...")
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
                cid = None
                for p in ['cid', 'code', 'id', 'ID', 'CATNUM']:
                    if m.HasProp(p):
                        cid = m.GetProp(p).strip()
                        break
                if not cid:
                    cid = f"{lib_name}_{count+1}"
                    
                records.append({
                    'compound_id': cid,
                    'library': lib_name,
                    'smiles': smiles
                })
                count += 1
        except Exception:
            continue
    print(f"  -> Added {count} unique structures from {lib_name}.")

df_compounds = pd.DataFrame(records)
print(f"\nTotal Unique Compounds Extracted across ALL Libraries: {len(df_compounds)}")

# Load screening matrix from yeast_bioactivity_multitask_sdf.csv or daily dumps
df_screen = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_screen.columns if c not in {'smiles', 'compound_id', 'library'}]

df_master = pd.merge(df_compounds, df_screen[['smiles'] + tasks], on='smiles', how='left')

output_path = "data/yeast_bioactivity_master_all_libraries.csv"
df_master.to_csv(output_path, index=False)
print(f"\nMaster Dataset Saved to {output_path}: {len(df_master)} compounds x {len(tasks)} strain tasks.")
print("Library breakdown:\n", df_master['library'].value_counts().to_string())
