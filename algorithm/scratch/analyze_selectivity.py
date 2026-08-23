import os
import sys
import json
import pandas as pd
import numpy as np

sys.path.append(".")
from src.predict_bioactivity import predict_bioactivity

print("=================================================================")
print("SELECTIVITY ANALYSIS ON FULL 4,408 COMPOUND MODEL (4 LIBRARIES)")
print("=================================================================")

target_compounds = {
    'Ciclosporine': 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C',
    'Mebendazole': 'COC(=O)Nc1nc2ccc(C(=O)c3ccccc3)cc2[nH]1',
    'Fluconazole': 'C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O',
    'Chrysin': 'C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)O)O'
}

smiles_list = list(target_compounds.values())
names_list = list(target_compounds.keys())

# Get model predictions for 4 target compounds
results = predict_bioactivity(smiles_list, top_k=102)

print("\n--- 4 TARGET COMPOUNDS SELECTIVITY CHECK ---")
for name, res in zip(names_list, results):
    all_scores = res['all_scores']
    active_strains = [strain for strain, prob in all_scores.items() if prob >= 0.50]
    num_active = len(active_strains)
    is_selective = (4 <= num_active <= 68)
    
    print(f"\nCompound: {name}")
    print(f"  - Predicted Active Strains (P >= 0.50): {num_active} / 102 ({num_active/102*100:.1f}%)")
    print(f"  - Satisfies Selectivity (4 <= k <= 68): {'YES (Selective)' if is_selective else 'NO'}")

# Evaluate selectivity across ALL 4,408 dataset compounds
df_dataset = pd.read_csv("data/yeast_bioactivity_all_libraries.csv")
all_smiles = df_dataset['smiles'].tolist()
all_cids = df_dataset['compound_id'].tolist()
all_libs = df_dataset['library'].tolist()

batch_size = 200
all_selectivity_results = []

for i in range(0, len(all_smiles), batch_size):
    batch_smiles = all_smiles[i:i+batch_size]
    batch_cids = all_cids[i:i+batch_size]
    batch_libs = all_libs[i:i+batch_size]
    
    batch_res = predict_bioactivity(batch_smiles, top_k=102)
    if not batch_res:
        continue
        
    for cid, lib, res in zip(batch_cids, batch_libs, batch_res):
        s = res['input_smiles']
        scores = res['all_scores']
        actives = [strain for strain, p in scores.items() if p >= 0.50]
        num_act = len(actives)
        is_sel = (4 <= num_act <= 68)
        
        all_selectivity_results.append({
            'compound_id': cid,
            'library': lib,
            'smiles': s,
            'num_active_strains': num_act,
            'percent_active': round(num_act / 102 * 100, 1),
            'is_selective': is_sel
        })

df_sel = pd.DataFrame(all_selectivity_results)
df_sel.to_csv("reports/compound_selectivity_analysis_full.csv", index=False)

selective_df = df_sel[df_sel['is_selective'] == True]
print(f"\n=================================================================")
print(f"Total Dataset Compounds Evaluated: {len(df_sel)}")
print(f"Selective Compounds (4 <= Active Strains <= 68): {len(selective_df)} / {len(df_sel)} ({len(selective_df)/len(df_sel)*100:.1f}%)")
print(f"Non-Selective (< 4 Active Strains): {(df_sel['num_active_strains'] < 4).sum()}")
print(f"Pan-Bioactive (> 68 Active Strains): {(df_sel['num_active_strains'] > 68).sum()}")
print("=================================================================")
