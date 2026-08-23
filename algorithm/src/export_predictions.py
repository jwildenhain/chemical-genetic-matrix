import os
import pandas as pd
import numpy as np

# Apply NumPy patches for older DeepChem/TensorFlow compatibility
np.object = object
np.bool = bool
np.int = int
np.float = float
np.typeDict = np.sctypeDict

import deepchem as dc
import tensorflow as tf

def export_probability_matrix():
    print("Loading filtered dataset...")
    df = pd.read_csv("data/multitask_dataset_filtered_binary.csv")
    
    strains = [c for c in df.columns if c not in ('compound_id', 'smiles')]
    
    print("Featurizing SMILES with ECFP4...")
    featurizer = dc.feat.CircularFingerprint(size=2048, radius=2)
    features = featurizer.featurize(df['smiles'].values)
    
    # Create DeepChem dataset (labels are dummy)
    labels = np.zeros((len(df), len(strains)))
    dataset = dc.data.NumpyDataset(X=features, y=labels)
    
    print(f"Loading MultitaskClassifier from models/deepchem_multitask...")
    model = dc.models.MultitaskClassifier(
        n_tasks=len(strains),
        n_features=2048,
        layer_sizes=[1024, 512, 256],
        dropouts=0.2,
        learning_rate=0.0005,
        model_dir="models/deepchem_multitask"
    )
    
    # Ensure model is restored
    model.restore()
    
    print("Running predictions...")
    # Predict outputs shape: (num_samples, num_tasks, 2)
    preds = model.predict(dataset)
    
    # Extract probability of being active (class 1)
    # preds has shape (num_samples, num_tasks, 2)
    prob_active = preds[:, :, 1]
    
    print("Creating probability matrix...")
    prob_df = pd.DataFrame(prob_active, columns=strains)
    prob_df.insert(0, 'compound_id', df['compound_id'].values)
    prob_df.insert(1, 'smiles', df['smiles'].values)
    
    out_path = "data/probability_matrix.csv"
    prob_df.to_csv(out_path, index=False)
    print(f"Exported probability matrix to {out_path}")

if __name__ == "__main__":
    # Suppress TF warnings
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    export_probability_matrix()
