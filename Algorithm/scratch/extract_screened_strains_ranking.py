import os
import sys
import json
import pandas as pd
import numpy as np

sys.path.append(".")
from src.predict_bioactivity import predict_bioactivity

compounds = {
    'Fluconazole': 'C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O',
    'Chrysin': 'C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O',
    'Ciclosporine': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
    'Mebendazole': 'COC(=O)Nc1nc2ccc(C(=O)c3ccccc3)cc2[nH]1'
}

smiles_list = list(compounds.values())
names_list = list(compounds.keys())

# Get model predictions
results = predict_bioactivity(smiles_list, top_k=102)

df_dataset = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_dataset.columns if c not in {'smiles', 'compound_id', 'library'}]

screened_report = {}

for name, res in zip(names_list, results):
    s = res['input_smiles']
    all_scores = res['all_scores'] # task -> model predicted prob
    
    exp_matches = df_dataset[df_dataset['smiles'] == s]
    
    if len(exp_matches) == 0:
        screened_report[name] = {
            'status': 'Not Present in 4-Library Experimental Screen (Spectrum, SPECMTS3, Maybridge1000, Lopac)',
            'compound_id': 'N/A',
            'active_strains': [],
            'inactive_strains_sample': []
        }
        continue
        
    row = exp_matches.iloc[0]
    cid = str(row['compound_id'])
    lib = str(row['library'])
    
    active_strains = []
    inactive_strains = []
    
    for t in tasks:
        val = row[t]
        if not pd.isna(val):
            exp_val = int(val)
            pred_p = float(all_scores.get(t, 0.0))
            item = {
                'strain': t,
                'experimental_label': exp_val,
                'model_predicted_probability': round(pred_p, 4)
            }
            if exp_val == 1:
                active_strains.append(item)
            else:
                inactive_strains.append(item)
                
    # Sort experimental active strains by model prediction
    active_strains = sorted(active_strains, key=lambda x: x['model_predicted_probability'], reverse=True)
    # Sort experimental inactive strains by model prediction (lowest first)
    inactive_strains = sorted(inactive_strains, key=lambda x: x['model_predicted_probability'])
    
    screened_report[name] = {
        'status': 'Present in Screen',
        'compound_id': cid,
        'library': lib,
        'total_screened_strains': len(active_strains) + len(inactive_strains),
        'total_active_strains': len(active_strains),
        'total_inactive_strains': len(inactive_strains),
        'actual_active_strains': active_strains,
        'sample_bottom_inactive_strains': inactive_strains[:5],
        'sample_top_inactive_strains': inactive_strains[-5:]
    }

with open("reports/screened_data_ranking.json", "w") as f:
    json.dump(screened_report, f, indent=2)

print("Saved screened data ranking report to reports/screened_data_ranking.json.")
