import os
import re
import json
import pandas as pd
from rdkit import Chem

print("=================================================================")
print("EXTRACTING MOLECULAR STRUCTURES FROM CHEMGRID SDF DATABASE")
print("=================================================================")

multitask_file = "data/yeast_bioactivity_multitask.csv"
df = pd.read_csv(multitask_file)
target_ids = set(df['compound_id'])
print(f"Loaded dataset {multitask_file} with {len(df)} compounds.")

sdfs = [
    '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Lopac.sdf',
    '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_Maybridge1000.sdf',
    '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_SPECMTS3.sdf',
    '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTRUM_ED.sdf',
    '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/cgm_ChemGRID2011.sdf'
]

sdf_smiles_map = {}
sdf_mol_map = {}

def extract_id_variants(val):
    if not val:
        return []
    s = str(val).strip()
    res = {s}
    s_clean = re.sub(r'\s+', '', s)
    res.add(s_clean)
    num = re.sub(r'[^0-9]', '', s_clean)
    if num:
        res.add(num)
        res.add(num.zfill(8))
        if 'LOP' in s_clean:
            res.add('LOP' + num)
            res.add('LOP' + num.zfill(8))
            res.add('LOPAC ' + num.zfill(5))
        elif 'SPE' in s_clean:
            res.add('SPE' + num)
            res.add('SPE' + num.zfill(8))
            res.add('SPECTRUM ' + num.zfill(8))
        elif 'MAY' in s_clean:
            res.add('MAY' + num)
            res.add('MAY' + num.zfill(8))
    return res

for sf in sdfs:
    if not os.path.exists(sf):
        continue
    print(f"Reading SDF: {os.path.basename(sf)}...")
    suppl = Chem.SDMolSupplier(sf, sanitize=False)
    for m in suppl:
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
            smiles = Chem.MolToSmiles(m)
        except Exception:
            continue
            
        if not smiles:
            continue
            
        prop_vals = []
        for p in ['cid', 'id', 'ID', 'code', 'CATNUM', 'product_name']:
            if m.HasProp(p):
                prop_vals.append(m.GetProp(p))
                
        for p_val in prop_vals:
            for variant in extract_id_variants(p_val):
                if variant in target_ids and variant not in sdf_smiles_map:
                    sdf_smiles_map[variant] = smiles
                    sdf_mol_map[variant] = m

matched_count = len(sdf_smiles_map)
print(f"\nTotal compounds successfully matched from SDF structures: {matched_count} / {len(df)} ({matched_count/len(df)*100:.2f}%)")

# Update SMILES column in dataframe with SDF-derived SMILES where available
updated_count = 0
for idx, row in df.iterrows():
    cid = row['compound_id']
    if cid in sdf_smiles_map:
        df.at[idx, 'smiles'] = sdf_smiles_map[cid]
        updated_count += 1

output_file = "data/yeast_bioactivity_multitask_sdf.csv"
df.to_csv(output_file, index=False)
print(f"Saved SDF-validated dataset to {output_file} with {updated_count} SDF structure updates.")
