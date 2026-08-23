import os
import sys
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_max_pool, global_add_pool
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
import pandas as pd
import numpy as np
import deepchem as dc
from tqdm import tqdm

print("=================================================================")
print("HYBRID GCN-ECFP MODEL (EARLY FUSION)")
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

print("Featurizing SMILES (MolGraph and ECFP)...")
gcn_featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
gcn_features = gcn_featurizer.featurize(df['smiles'].values)

ecfp_featurizer = dc.feat.CircularFingerprint(size=2048)
ecfp_features = ecfp_featurizer.featurize(df['smiles'].values)

dataset = []
for i, feat in enumerate(gcn_features):
    if not hasattr(feat, 'node_features') or feat.node_features is None: continue
    
    x = torch.tensor(feat.node_features, dtype=torch.float)
    if hasattr(feat, 'edge_index') and feat.edge_index is not None:
        edge_index = torch.tensor(feat.edge_index, dtype=torch.long)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        
    ecfp = torch.tensor(ecfp_features[i], dtype=torch.float).view(1, -1)
    
    y = torch.tensor(df[tasks].iloc[i].values, dtype=torch.float).view(1, -1)
    
    # Mask missing values
    mask = ~torch.isnan(y)
    y = torch.nan_to_num(y, 0.0)
    
    data = Data(x=x, edge_index=edge_index, ecfp=ecfp, y=y, mask=mask)
    dataset.append(data)

print(f"Successfully converted {len(dataset)} molecules to Hybrid objects.")

# Split dataset
np.random.seed(42)
indices = np.random.permutation(len(dataset))
train_idx = indices[:int(0.8*len(indices))]
val_idx = indices[int(0.8*len(indices)):int(0.9*len(indices))]
test_idx = indices[int(0.9*len(indices)):]

train_loader = DataLoader([dataset[i] for i in train_idx], batch_size=64, shuffle=True)
val_loader = DataLoader([dataset[i] for i in val_idx], batch_size=64)
test_loader = DataLoader([dataset[i] for i in test_idx], batch_size=64)

class HybridModel(nn.Module):
    def __init__(self, node_channels, gcn_hidden, ecfp_size, ecfp_hidden, out_channels):
        super(HybridModel, self).__init__()
        
        # Branch 1: Graph Convolution (DeepChem GCN architecture)
        self.conv1 = GCNConv(node_channels, gcn_hidden)
        self.conv2 = GCNConv(gcn_hidden, gcn_hidden)
        
        self.gating = nn.Sequential(
            nn.Linear(gcn_hidden, 1),
            nn.Sigmoid()
        )
        
        # Branch 2: ECFP Dense Network (Extremely Reduced to prevent overfitting)
        self.ecfp_dense = nn.Sequential(
            nn.Linear(ecfp_size, 128),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(128, ecfp_hidden),
            nn.ReLU(),
            nn.Dropout(p=0.5)
        )
        
        # Joint MLP Readout
        joint_size = (2 * gcn_hidden) + ecfp_hidden
        self.mlp = nn.Sequential(
            nn.Linear(joint_size, 64),
            nn.ReLU(),
            nn.Dropout(p=0.5),
            nn.Linear(64, out_channels)
        )

    def forward(self, x, edge_index, ecfp, batch):
        # 1. Graph Forward
        x_gcn = self.conv1(x, edge_index)
        x_gcn = F.relu(x_gcn)
        x_gcn = F.dropout(x_gcn, p=0.3, training=self.training)
        x_gcn = self.conv2(x_gcn, edge_index)
        x_gcn = F.relu(x_gcn)
        x_gcn = F.dropout(x_gcn, p=0.3, training=self.training)
        
        gate_weights = self.gating(x_gcn)
        weighted_x = x_gcn * gate_weights
        sum_pool = global_add_pool(weighted_x, batch)
        max_pool = global_max_pool(x_gcn, batch)
        graph_repr = torch.cat([sum_pool, max_pool], dim=-1) # Shape: [batch, 2 * gcn_hidden]
        
        # 2. ECFP Forward
        ecfp_repr = self.ecfp_dense(ecfp) # Shape: [batch, ecfp_hidden]
        
        # 3. Concatenate (Early Fusion)
        joint_repr = torch.cat([graph_repr, ecfp_repr], dim=-1)
        
        # 4. Joint MLP
        return self.mlp(joint_repr)

# Severely reduced capacity: GCN 32, ECFP 32
model = HybridModel(node_channels=30, gcn_hidden=32, ecfp_size=2048, ecfp_hidden=32, out_channels=len(tasks)).to(device)
# Added strong weight decay for regularization
optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-3)

pos_weight = torch.ones([len(tasks)], device=device) * 5.0
criterion = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pos_weight)

print(f"Training Hybrid PyG+ECFP model for 60 epochs on {device}...")
for epoch in range(1, 61):
    model.train()
    total_loss = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        out = model(data.x, data.edge_index, data.ecfp, data.batch)
        loss = criterion(out, data.y)
        loss = (loss * data.mask).sum() / (data.mask.sum() + 1e-8)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * data.num_graphs
    
    if epoch % 5 == 0 or epoch == 1:
        print(f"Epoch {epoch:2d}/50 | Loss: {total_loss / len(train_loader.dataset):.4f}")

def evaluate(loader):
    model.eval()
    y_true, y_pred, y_mask = [], [], []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            out = model(data.x, data.edge_index, data.ecfp, data.batch)
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

print("\n--- HYBRID MODEL EVALUATION ---")
print(f"Train ROC-AUC: {train_roc:.4f} | PR-AUC: {train_prc:.4f}")
print(f"Valid ROC-AUC: {val_roc:.4f} | PR-AUC: {val_prc:.4f}")
print(f"Test  ROC-AUC: {test_roc:.4f} | PR-AUC: {test_prc:.4f}")

os.makedirs("models/hybrid_model", exist_ok=True)
torch.save(model.state_dict(), "models/hybrid_model/hybrid.pt")

print("\nHybrid model successfully trained and saved!")
