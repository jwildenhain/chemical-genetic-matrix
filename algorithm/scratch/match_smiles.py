import pandas as pd
from rdkit import Chem

df = pd.read_csv('data/yeast_bioactivity_extended_master.csv')

def get_canonical(smiles):
    if pd.isna(smiles): return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol: return Chem.MolToSmiles(mol)
    except: pass
    return None

print("Canonicalizing target SMILES...")
meb_smi = "COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CC=C3"
cyc_smi = "C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C"

meb_can = get_canonical(meb_smi)
cyc_can = get_canonical(cyc_smi)

print(f"Mebendazole Canonical: {meb_can}")
print(f"Cyclosporine Canonical: {cyc_can}")

print("Checking against database...")
matches = []
for idx, row in df.iterrows():
    can = get_canonical(row['smiles'])
    if can == meb_can:
        print(f"FOUND MEBENDAZOLE! ID: {row['compound_id']}")
        matches.append(('Mebendazole', row))
    elif can == cyc_can:
        print(f"FOUND CYCLOSPORINE! ID: {row['compound_id']}")
        matches.append(('Cyclosporine', row))

for name, row in matches:
    print(f"\n{name} Bioactivity:")
    active = []
    for col in row.index:
        if col not in ('compound_id', 'library', 'smiles') and row[col] == 1.0:
            active.append(col)
    print(active)
