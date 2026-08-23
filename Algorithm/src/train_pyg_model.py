import os
import sys
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
import pandas as pd
import numpy as np
import deepchem as dc
from tqdm import tqdm

print("=================================================================")
print("PYTORCH GEOMETRIC GCN MODEL (NATIVE PyG ON RTX 5080)")
print("=================================================================")

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using {device} for training.")

dataset_file = "data/predictive_matrix_cleaned.csv"
if not os.path.exists(dataset_file):
    print(f"Error: {dataset_file} not found.")
    sys.exit(1)

df = pd.read_csv(dataset_file)
tasks = [c for c in df.columns if c not in ('smiles', 'compound_id')]

print(f"Total Compounds: {len(df)}")
print(f"Total Tasks (Strains): {len(tasks)}")

print("Featurizing SMILES...")
featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
features = featurizer.featurize(df['smiles'].values)

# Convert to PyG Data objects
dataset = []
for i, feat in enumerate(features):
    if not hasattr(feat, 'node_features') or feat.node_features is None: continue
    x = torch.tensor(feat.node_features, dtype=torch.float)
    if hasattr(feat, 'edge_index') and feat.edge_index is not None:
        edge_index = torch.tensor(feat.edge_index, dtype=torch.long)
    else:
        # For single atom molecules with no edges
        edge_index = torch.empty((2, 0), dtype=torch.long)
        
    y = torch.tensor(df[tasks].iloc[i].values, dtype=torch.float).view(1, -1)
    
    # Mask missing values
    mask = ~torch.isnan(y)
    y = torch.nan_to_num(y, 0.0)
    
    data = Data(x=x, edge_index=edge_index, y=y, mask=mask)
    dataset.append(data)

print(f"Successfully converted {len(dataset)} molecules to graphs.")

# Split dataset
np.random.seed(42)
indices = np.random.permutation(len(dataset))
train_idx = indices[:int(0.8*len(indices))]
val_idx = indices[int(0.8*len(indices)):int(0.9*len(indices))]
test_idx = indices[int(0.9*len(indices)):]

train_loader = DataLoader([dataset[i] for i in train_idx], batch_size=64, shuffle=True)
val_loader = DataLoader([dataset[i] for i in val_idx], batch_size=64)
test_loader = DataLoader([dataset[i] for i in test_idx], batch_size=64)

from torch_geometric.nn import GCNConv, global_max_pool

class GCN(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GCN, self).__init__()
        # DeepChem GCN typically uses 2 GCN layers by default
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        
        # Gating function for weighted sum
        self.gating = nn.Sequential(
            nn.Linear(hidden_channels, 1),
            nn.Sigmoid()
        )
        
        # MLP for final prediction (concat of weighted sum and max pool -> hidden -> output)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_channels * 2, hidden_channels),
            nn.ReLU(),
            nn.Dropout(p=0.2),
            nn.Linear(hidden_channels, out_channels)
        )

    def forward(self, x, edge_index, batch):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        
        # 1. Gated Weighted Sum Pooling
        gate_weights = self.gating(x)
        weighted_x = x * gate_weights
        from torch_geometric.nn import global_add_pool
        sum_pool = global_add_pool(weighted_x, batch)
        
        # 2. Max Pooling
        max_pool = global_max_pool(x, batch)
        
        # 3. Concatenate
        graph_repr = torch.cat([sum_pool, max_pool], dim=-1)
        
        # 4. MLP Readout
        return self.mlp(graph_repr)

model = GCN(in_channels=30, hidden_channels=128, out_channels=len(tasks)).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
# Use pos_weight to handle class imbalance (mostly inactive)
# The ratio of inactive to active is roughly 10:1 in screening data
pos_weight = torch.ones([len(tasks)], device=device) * 5.0
criterion = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pos_weight)

print(f"Training native PyG GCN model for 40 epochs on {device}...")
for epoch in range(1, 41):
    model.train()
    total_loss = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.batch)
        loss = criterion(out, data.y)
        loss = (loss * data.mask).sum() / (data.mask.sum() + 1e-8)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * data.num_graphs
    
    if epoch % 5 == 0 or epoch == 1:
        print(f"Epoch {epoch:2d}/40 | Loss: {total_loss / len(train_loader.dataset):.4f}")

def evaluate(loader):
    model.eval()
    y_true, y_pred, y_mask = [], [], []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            out = model(data.x, data.edge_index, data.batch)
            y_true.append(data.y.cpu())
            y_pred.append(torch.sigmoid(out).cpu())
            y_mask.append(data.mask.cpu())
            
    y_true = torch.cat(y_true, dim=0).numpy()
    y_pred = torch.cat(y_pred, dim=0).numpy()
    y_mask = torch.cat(y_mask, dim=0).numpy()
    
    roc_aucs, prc_aucs = [], []
    for i in range(y_true.shape[1]):
        mask_i = y_mask[:, i]
        if mask_i.sum() > 0:
            yt, yp = y_true[mask_i, i], y_pred[mask_i, i]
            if len(np.unique(yt)) > 1:
                roc_aucs.append(roc_auc_score(yt, yp))
                p, r, _ = precision_recall_curve(yt, yp)
                prc_aucs.append(auc(r, p))
    
    return np.mean(roc_aucs) if roc_aucs else 0.0, np.mean(prc_aucs) if prc_aucs else 0.0

train_roc, train_prc = evaluate(train_loader)
val_roc, val_prc = evaluate(val_loader)
test_roc, test_prc = evaluate(test_loader)

print("\n--- PYTORCH GEOMETRIC EVALUATION ---")
print(f"Train ROC-AUC: {train_roc:.4f} | PR-AUC: {train_prc:.4f}")
print(f"Valid ROC-AUC: {val_roc:.4f} | PR-AUC: {val_prc:.4f}")
print(f"Test  ROC-AUC: {test_roc:.4f} | PR-AUC: {test_prc:.4f}")

os.makedirs("models/deepchem_modern_gcn", exist_ok=True)
torch.save(model.state_dict(), "models/deepchem_modern_gcn/pyg_model.pt")

summary = {
    'device': 'cuda' if torch.cuda.is_available() else 'cpu',
    'train_roc_auc': round(train_roc, 4),
    'valid_roc_auc': round(val_roc, 4),
    'test_roc_auc': round(test_roc, 4),
    'test_prc_auc': round(test_prc, 4)
}
with open("reports/modern_gcn_model_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\nGraph model successfully trained and saved!")
