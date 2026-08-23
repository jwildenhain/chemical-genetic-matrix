import os
import sys
import argparse
import numpy as np

# Apply NumPy patches for older DeepChem/TensorFlow compatibility
np.object = object
np.bool = bool
np.int = int
np.float = float
np.typeDict = np.sctypeDict

import tensorflow as tf
if not hasattr(tf.keras.optimizers, 'legacy'):
    tf.keras.optimizers.legacy = tf.keras.optimizers

import pandas as pd
import deepchem as dc

def predict_molecule(smiles, model_dir, mapping_csv, dataset_csv, show_all=False):
    print(f"Loading highly predictive strains from base dataset...")
    df = pd.read_csv(dataset_csv, nrows=1)
    tasks = [c for c in df.columns if c not in ('compound_id', 'smiles')]
    
    print(f"Initializing MultitaskClassifier for {len(tasks)} target tasks...")
    model = dc.models.MultitaskClassifier(
        n_tasks=len(tasks),
        n_features=2048,
        layer_sizes=[1024, 1024, 512],
        dropouts=[0.3, 0.3, 0.2],
        model_dir=model_dir
    )
    model.restore()
    
    print(f"Featurizing input SMILES: {smiles}")
    featurizer = dc.feat.CircularFingerprint(size=2048, radius=2)
    features = featurizer.featurize([smiles])
    
    if len(features) == 0 or features[0] is None:
        print("Error: Could not featurize the given SMILES string. It may be invalid.")
        sys.exit(1)
        
    dataset = dc.data.NumpyDataset(X=features, y=np.zeros((1, len(tasks))))
    
    print("Running inference...")
    # MultitaskClassifier classification output shape: (num_samples, num_tasks, 2)
    preds = model.predict(dataset)
    if isinstance(preds, list):
        preds = np.array(preds)
        
    if preds.ndim == 3 and preds.shape[2] == 2:
        prob_active = preds[0, :, 1]
    else:
        prob_active = preds[0]
    
    print("\nLoading Strain -> ORF mapping...")
    mapping_df = pd.read_csv(mapping_csv)
    # create a dictionary mapping strain -> ORF
    strain_to_orf = dict(zip(mapping_df['Strain'], mapping_df['ORF']))
    
    results = []
    for idx, strain in enumerate(tasks):
        prob = prob_active[idx]
        if prob >= 0.50 or show_all:
            orf = strain_to_orf.get(strain, "No ORF (WildType/Control)")
            results.append((strain, orf, prob))
            
    # Sort by highest probability
    results.sort(key=lambda x: x[2], reverse=True)
    
    print(f"\n=======================================================")
    print(f"PREDICTION RESULTS FOR: {smiles}")
    print(f"=======================================================")
    if len(results) == 0:
        print("No active targets predicted (Probability < 50% across all 105 strains).")
    else:
        print(f"Showing results for {len(results)} strains:")
        print(f"{'Strain':<25} | {'ORF':<15} | {'Probability':<10}")
        print("-" * 55)
        for strain, orf, prob in results:
            print(f"{strain:<25} | {orf:<15} | {prob:.2%}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict yeast strain bioactivity from SMILES using MultitaskClassifier.")
    parser.add_argument("--smiles", type=str, required=True, help="SMILES string of the molecule")
    parser.add_argument("--model_dir", type=str, default="model", help="Path to DeepChem model directory")
    parser.add_argument("--mapping_csv", type=str, default="../data/strain_orf_mapping.csv", help="Path to strain_orf_mapping.csv")
    parser.add_argument("--dataset_csv", type=str, default="../data/predictive_matrix_cleaned.csv", help="Path to base dataset to get task order")
    parser.add_argument("--all", action="store_true", help="Show probabilities for all strains, not just > 50%")
    
    args = parser.parse_args()
    
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    predict_molecule(args.smiles, args.model_dir, args.mapping_csv, args.dataset_csv, args.all)
