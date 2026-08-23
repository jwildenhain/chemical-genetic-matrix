import os
import re
import pandas as pd
from rdkit import Chem

print("=================================================================")
print("AUDITING TOTAL COMPOUND COUNTS ACROSS TARGET LIBRARIES")
print("=================================================================")

sdfs = {
    'Spectrum': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTRUM_ED.sdf',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_SPECMTS3.sdf',
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Lopac.sdf'
}

lib_counts = {}
all_lib_ids = set()

for lib_name, sf in sdfs.items():
    if os.path.exists(sf):
        suppl = Chem.SDMolSupplier(sf, sanitize=False)
        mols = [m for m in suppl if m is not None]
        lib_counts[lib_name] = len(mols)
        print(f"Library: {lib_name:15s} | Raw SDF Structures: {len(mols)}")
        for m in mols:
            for p in ['cid', 'id', 'ID', 'code', 'CATNUM']:
                if m.HasProp(p):
                    all_lib_ids.add(m.GetProp(p).strip())

total_raw_compounds = sum(lib_counts.values())
print(f"\nTotal Raw Compounds in Target Libraries (SDF files): {total_raw_compounds}")
print(f"Total Unique Compounds in Target Libraries:          {len(all_lib_ids)}")

df_screen = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
screened_cids = set(df_screen['compound_id'])
print(f"\nCompounds in 102-Strain Bioactivity Dataset:       {len(screened_cids)}")
print(f"Unscreened / Partial Library Compounds:              {len(all_lib_ids) - len(screened_cids)}")

lib_breakdown = df_screen['library'].value_counts()
print("\nDataset Compounds Breakdown by Library:")
print(lib_breakdown.to_string())
