import os
import sys
import json
import argparse
import numpy as np
import pandas as pd
import deepchem as dc
from rdkit import Chem

def predict_bioactivity(smiles_list, model_dir="models/deepchem_multitask", dataset_file="data/yeast_bioactivity_multitask.csv", top_k=5):
    if not os.path.exists(dataset_file):
        print(f"Error: Dataset {dataset_file} not found.")
        return None
    
    df = pd.read_csv(dataset_file)
    non_task_cols = {'smiles', 'compound_id', 'library'}
    tasks = [c for c in df.columns if c not in non_task_cols]
    
    valid_smiles = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m:
            valid_smiles.append(Chem.MolToSmiles(m))
        else:
            print(f"Warning: Invalid SMILES ignored: {s}")
            
    if not valid_smiles:
        print("No valid SMILES provided.")
        return None
    
    print(f"Featurizing {len(valid_smiles)} input molecule(s)...")
    featurizer = dc.feat.CircularFingerprint(size=2048, radius=2)
    features = featurizer.featurize(valid_smiles)
    
    # Initialize dummy targets y with proper multi-task shape (num_samples, num_tasks)
    dummy_y = np.zeros((len(valid_smiles), len(tasks)))
    dataset = dc.data.NumpyDataset(
        X=features,
        y=dummy_y,
        ids=valid_smiles
    )
    
    print("Loading trained MultitaskClassifier model...")
    model = dc.models.MultitaskClassifier(
        n_tasks=len(tasks),
        n_features=2048,
        layer_sizes=[1024, 512, 256],
        model_dir=model_dir
    )
    
    print("Predicting bioactivity probabilities across yeast strains...")
    preds = model.predict(dataset) # shape: (num_samples, num_tasks, 2)
    
    results = []
    for i, s in enumerate(valid_smiles):
        active_probs = preds[i, :, 1] # Probability of active class (1)
        strain_scores = dict(zip(tasks, active_probs))
        sorted_strains = sorted(strain_scores.items(), key=lambda x: x[1], reverse=True)
        
        results.append({
            'input_smiles': s,
            'top_active_targets': sorted_strains[:top_k],
            'all_scores': strain_scores
        })
        
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict Bioactivity Across Yeast Strains for Input Molecules")
    parser.add_argument("smiles", nargs="*", help="One or more SMILES strings to evaluate")
    parser.add_argument("--top_k", type=int, default=5, help="Number of top sensitive yeast strain targets to report")
    args = parser.parse_args()

    smiles_to_test = args.smiles if args.smiles else [
        "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34", # Testosterone
        "CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC", # Cocaine
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" # Caffeine
    ]
    
    res = predict_bioactivity(smiles_to_test, top_k=args.top_k)
    if res:
        for r in res:
            print("\n===================================================")
            print(f"SMILES: {r['input_smiles']}")
            print(f"Top {args.top_k} Predicted Bioactive Yeast Target Strains:")
            for strain, prob in r['top_active_targets']:
                print(f"  - {strain:20s}: Bioactivity Probability = {prob:.4f}")
