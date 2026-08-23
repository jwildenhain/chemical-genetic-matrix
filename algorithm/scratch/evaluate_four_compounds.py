import os
import sys
import json
import pandas as pd
import numpy as np
import deepchem as dc
from rdkit import Chem

sys.path.append(".")
from src.predict_bioactivity import predict_bioactivity

compounds = {
    'Ciclosporine': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
    'Mebendazole': 'COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CC=C3',
    'Fluconazole': 'C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O',
    'Chrysin': 'C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O'
}

print("=================================================================")
print("RUNNING BIOACTIVITY MODEL FOR: Ciclosporine, Mebendazole, Fluconazole, Chrysin")
print("=================================================================")

smiles_list = list(compounds.values())
names_list = list(compounds.keys())

# Run predictions
results = predict_bioactivity(smiles_list, top_k=10)

df_dataset = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_dataset.columns if c not in {'smiles', 'compound_id', 'library'}]

out_dict = {}

for name, res in zip(names_list, results):
    s = res['input_smiles']
    top_targets = [(strain, float(prob)) for strain, prob in res['top_active_targets']]
    all_scores = {k: float(v) for k, v in res['all_scores'].items()}
    
    # Check if compound is in experimental dataset
    exp_matches = df_dataset[df_dataset['smiles'] == s]
    exp_data = {}
    cid = "N/A (External Candidate)"
    lib = "N/A"
    
    if len(exp_matches) > 0:
        row = exp_matches.iloc[0]
        cid = str(row['compound_id'])
        lib = str(row['library'])
        for t in tasks:
            if not pd.isna(row[t]):
                exp_data[t] = int(row[t])
                
    out_dict[name] = {
        'name': name,
        'smiles': s,
        'compound_id': cid,
        'library': lib,
        'top_predicted_strains': top_targets,
        'experimental_screened_strains': len(exp_data),
        'experimental_active_strains': [t for t, v in exp_data.items() if v == 1]
    }

with open("reports/four_compounds_profile.json", "w") as f:
    json.dump(out_dict, f, indent=2)

print("\nSaved 4-compound evaluation profiles to reports/four_compounds_profile.json.")
