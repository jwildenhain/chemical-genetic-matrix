import os
import numpy as np
import pandas as pd
import torch
import deepchem as dc
from sklearn.metrics import roc_auc_score, average_precision_score

# Add root directory to sys.path to allow importing from cgm_predictor
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cgm_predictor.models import CGMPredictor

def main():
    print("Evaluating Ensemble Model (Late Fusion of ECFP + GCN)...")
    
    # Force CPU to match previous setups
    device = torch.device('cpu')

    # Load data
    df = pd.read_csv('data/predictive_matrix_cleaned.csv')
    
    # We only care about evaluation, we'll use the DeepChem featurizer to split the data deterministically
    # just as we did before.
    
    tasks = list(df.columns)
    if 'compound_id' in tasks:
        tasks.remove('compound_id')
    if 'smiles' in tasks:
        tasks.remove('smiles')
    
    # 1. We need the same deterministic ScaffoldSplitter that was used for the other models
    # Wait, the easiest way is to use DeepChem's featurizer to re-create the split
    print("Featurizing data to get identical splits...")
    featurizer = dc.feat.CircularFingerprint(size=2048)
    loader = dc.data.CSVLoader(
        tasks=tasks,
        feature_field="smiles",
        id_field="smiles",
        featurizer=featurizer
    )
    dataset = loader.create_dataset('data/predictive_matrix_cleaned.csv')
    
    splitter = dc.splits.ScaffoldSplitter()
    train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(
        dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1
    )
    
    # We will extract the SMILES and labels from the datasets
    # Actually, dc datasets keep the IDs (which are SMILES strings for CSVLoader)
    
    train_smiles = train_dataset.ids
    train_y = train_dataset.y
    valid_smiles = valid_dataset.ids
    valid_y = valid_dataset.y
    test_smiles = test_dataset.ids
    test_y = test_dataset.y
    
    # Load Ensemble Model
    print("Loading ensemble predictor...")
    predictor = CGMPredictor(model_type='ensemble', n_tasks=len(tasks))
    
    def evaluate_split(smiles_list, y_true, split_name):
        print(f"Predicting on {split_name} split ({len(smiles_list)} samples)...")
        # predictor.predict() takes a single SMILES, so we iterate
        y_pred = []
        valid_y_true = []
        for i, s in enumerate(smiles_list):
            try:
                probs = predictor.predict(s)
                y_pred.append(probs)
                valid_y_true.append(y_true[i])
            except Exception:
                pass
                
        y_pred = np.array(y_pred)
        valid_y_true = np.array(valid_y_true)
        
        # Calculate metrics for each task, handling cases where only one class is present
        roc_aucs = []
        pr_aucs = []
        for i in range(valid_y_true.shape[1]):
            y_t = valid_y_true[:, i]
            y_p = y_pred[:, i]
            if len(np.unique(y_t)) > 1:
                roc_aucs.append(roc_auc_score(y_t, y_p))
                pr_aucs.append(average_precision_score(y_t, y_p))
                
        mean_roc = np.mean(roc_aucs) if roc_aucs else 0.0
        mean_pr = np.mean(pr_aucs) if pr_aucs else 0.0
        
        print(f"{split_name} ROC-AUC: {mean_roc:.4f} | PR-AUC: {mean_pr:.4f}")
        return mean_roc, mean_pr

    evaluate_split(train_smiles, train_y, "Train")
    evaluate_split(valid_smiles, valid_y, "Valid")
    evaluate_split(test_smiles, test_y, "Test")

if __name__ == "__main__":
    main()
