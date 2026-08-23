import os
import sys
import json
import pandas as pd
import numpy as np

sys.path.append(".")
from src.predict_bioactivity import predict_bioactivity

compounds = {
    'Ciclosporine': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
    'Mebendazole': 'COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CC=C3',
    'Fluconazole': 'C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O',
    'Chrysin': 'C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O'
}

smiles_list = list(compounds.values())
names_list = list(compounds.keys())

# Run model predictions across all 102 yeast strain tasks
results = predict_bioactivity(smiles_list, top_k=102)

df_dataset = pd.read_csv("data/yeast_bioactivity_multitask_sdf.csv")
tasks = [c for c in df_dataset.columns if c not in {'smiles', 'compound_id', 'library'}]

report_data = {}

for name, res in zip(names_list, results):
    s = res['input_smiles']
    all_scores = res['all_scores'] # strain_task -> prob
    
    # Sort all 102 strains by predicted probability descending
    sorted_strains = sorted(all_scores.items(), key=lambda x: x[1], reverse=True)
    top5_pred = sorted_strains[:5]
    bottom5_pred = sorted_strains[-5:]
    
    # Look up experimental data from dataset if present
    exp_matches = df_dataset[df_dataset['smiles'] == s]
    exp_labels = {}
    cid = "External Candidate"
    lib = "N/A"
    
    if len(exp_matches) > 0:
        row = exp_matches.iloc[0]
        cid = str(row['compound_id'])
        lib = str(row['library'])
        for t in tasks:
            if not pd.isna(row[t]):
                exp_labels[t] = int(row[t])
                
    exp_active_strains = [t for t, v in exp_labels.items() if v == 1]
    
    # For top5 predicted, attach experimental status
    top5_formatted = []
    for strain, prob in top5_pred:
        exp_val = exp_labels.get(strain, "Unscreened")
        top5_formatted.append({
            'strain': strain,
            'predicted_probability': round(float(prob), 4),
            'experimental_label': exp_val
        })
        
    bottom5_formatted = []
    for strain, prob in bottom5_pred:
        exp_val = exp_labels.get(strain, "Unscreened")
        bottom5_formatted.append({
            'strain': strain,
            'predicted_probability': round(float(prob), 4),
            'experimental_label': exp_val
        })
        
    exp_active_formatted = []
    for strain in exp_active_strains:
        prob = all_scores.get(strain, 0.0)
        exp_active_formatted.append({
            'strain': strain,
            'experimental_label': 1,
            'predicted_probability': round(float(prob), 4)
        })
    # Sort actual actives by predicted probability
    exp_active_formatted = sorted(exp_active_formatted, key=lambda x: x['predicted_probability'], reverse=True)
        
    report_data[name] = {
        'name': name,
        'compound_id': cid,
        'library': lib,
        'smiles': s,
        'top5_predicted': top5_formatted,
        'bottom5_predicted': bottom5_formatted,
        'actual_experimental_actives': exp_active_formatted
    }

with open("reports/top_bottom_strains_report.json", "w") as f:
    json.dump(report_data, f, indent=2)

print("Saved top/bottom strains report to reports/top_bottom_strains_report.json.")
