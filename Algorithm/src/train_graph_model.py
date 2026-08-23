import os
import sys
import json
import torch
import numpy as np

# Apply NumPy patches just in case
np.object = object
np.bool = bool
np.int = int
np.float = float
np.typeDict = np.sctypeDict

import pandas as pd
import deepchem as dc
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

print("=================================================================")
print("DEEPCHEM PYTORCH GCN MODEL (MODERN ENVIRONMENT)")
print("=================================================================")

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"PyTorch CUDA Available: {torch.cuda.is_available()}")

dataset_file = "data/predictive_matrix_cleaned.csv"
if not os.path.exists(dataset_file):
    print(f"Error: {dataset_file} not found.")
    sys.exit(1)

df = pd.read_csv(dataset_file)
tasks = [c for c in df.columns if c not in ('smiles', 'compound_id')]

print(f"Total Filtered Compounds: {len(df)}")
print(f"Total Highly Predictive Yeast Strains ({len(tasks)}): {tasks[:5]}...")

os.makedirs("models/deepchem_modern_gcn", exist_ok=True)
os.makedirs("reports", exist_ok=True)

print("\n=================================================================")
print("STEP 1: Featurizing Target Compounds with MolGraphConvFeaturizer")
print("=================================================================")

# MolGraphConvFeaturizer creates the PyTorch Geometric graphs
featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
loader = dc.data.CSVLoader(
    tasks=tasks,
    featurizer=featurizer,
    feature_field="smiles",
    id_field="smiles"
)

dataset = loader.create_dataset(dataset_file)
print(f"Featurized dataset shape: {dataset.get_shape()}")

print("\n=================================================================")
print("STEP 2: Scaffold Split & Class Balancing Weights")
print("=================================================================")

splitter = dc.splits.ScaffoldSplitter()
train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(
    dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1, seed=42
)

# Use BalancingTransformer to weight positive samples higher
transformer = dc.trans.BalancingTransformer(dataset=train_dataset)
train_dataset = transformer.transform(train_dataset)
valid_dataset = transformer.transform(valid_dataset)
test_dataset = transformer.transform(test_dataset)

print(f"Train Dataset size: {len(train_dataset)}")
print(f"Valid Dataset size: {len(valid_dataset)}")
print(f"Test Dataset size:  {len(test_dataset)}")

print("\n=================================================================")
print("STEP 3: Initializing & Training DeepChem GCNModel on RTX 5080")
print("=================================================================")

# GCNModel utilizes PyTorch Geometric directly on the GPU
model = dc.models.GCNModel(
    n_tasks=len(tasks),
    mode='classification',
    dropout=0.2,
    batch_size=64,
    learning_rate=0.001,
    model_dir="models/deepchem_modern_gcn",
    device=torch.device('cuda') if torch.cuda.is_available() else None
)

nb_epoch = 30
print(f"Training PyTorch Geometric GCN for {nb_epoch} epochs...")

loss_history = []
for epoch in range(1, nb_epoch + 1):
    loss = model.fit(train_dataset, nb_epoch=1)
    loss_history.append(float(loss))
    if epoch % 5 == 0 or epoch == 1:
        print(f"Epoch {epoch:2d}/{nb_epoch:2d} | Training Loss: {loss:.4f}")

print("\n=================================================================")
print("STEP 4: Evaluating Modern GCN Metrics")
print("=================================================================")

def compute_multitask_metrics(dset):
    y_true = dset.y
    w = dset.w
    raw_pred = model.predict(dset)
    if isinstance(raw_pred, list): raw_pred = np.array(raw_pred)
    y_pred_prob = raw_pred[:, :, 1] if raw_pred.ndim == 3 else raw_pred
        
    roc_aucs, prc_aucs = [], []
    for i in range(y_true.shape[1]):
        valid_mask = w[:, i] > 0
        yt = y_true[valid_mask, i]
        yp = y_pred_prob[valid_mask, i]
        if len(np.unique(yt)) > 1:
            roc_aucs.append(roc_auc_score(yt, yp))
            p, r, _ = precision_recall_curve(yt, yp)
            prc_aucs.append(auc(r, p))
            
    return {
        'roc_auc': float(np.nanmean(roc_aucs)) if roc_aucs else 0.0,
        'prc_auc': float(np.nanmean(prc_aucs)) if prc_aucs else 0.0
    }

train_scores = compute_multitask_metrics(train_dataset)
valid_scores = compute_multitask_metrics(valid_dataset)
test_scores  = compute_multitask_metrics(test_dataset)

print("\n--- PYTORCH GCN MULTI-TASK EVALUATION METRICS ---")
print(f"Train Set -> ROC-AUC: {train_scores['roc_auc']:.4f} | PR-AUC: {train_scores['prc_auc']:.4f}")
print(f"Valid Set -> ROC-AUC: {valid_scores['roc_auc']:.4f} | PR-AUC: {valid_scores['prc_auc']:.4f}")
print(f"Test Set  -> ROC-AUC: {test_scores['roc_auc']:.4f} | PR-AUC: {test_scores['prc_auc']:.4f}")

summary_dict = {
    'device': device,
    'total_compounds': len(df),
    'num_tasks': len(tasks),
    'train_roc_auc': round(train_scores['roc_auc'], 4),
    'valid_roc_auc': round(valid_scores['roc_auc'], 4),
    'test_roc_auc': round(test_scores['roc_auc'], 4),
    'test_prc_auc': round(test_scores['prc_auc'], 4)
}

with open("reports/modern_gcn_model_summary.json", "w") as sf:
    json.dump(summary_dict, sf, indent=2)

print("\nModern PyTorch GCN training completed successfully!")
