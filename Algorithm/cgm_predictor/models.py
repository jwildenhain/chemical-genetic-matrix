import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_max_pool, global_add_pool
from torch_geometric.data import Data
import deepchem as dc
import numpy as np

# ---------------------------------------------------------
# PyTorch Geometric (PyG) Model Definition
# ---------------------------------------------------------
class GCN(nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        
        self.gating = nn.Sequential(
            nn.Linear(hidden_channels, 1),
            nn.Sigmoid()
        )
        
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
        
        gate_weights = self.gating(x)
        weighted_x = x * gate_weights
        sum_pool = global_add_pool(weighted_x, batch)
        max_pool = global_max_pool(x, batch)
        
        graph_repr = torch.cat([sum_pool, max_pool], dim=-1)
        return self.mlp(graph_repr)

class HybridModel(nn.Module):
    def __init__(self, node_channels, gcn_hidden, ecfp_size, ecfp_hidden, out_channels):
        super(HybridModel, self).__init__()
        
        # Branch 1: Graph Convolution (DeepChem GCN architecture)
        from torch_geometric.nn import GCNConv
        self.conv1 = GCNConv(node_channels, gcn_hidden)
        self.conv2 = GCNConv(gcn_hidden, gcn_hidden)
        
        self.gating = nn.Sequential(
            nn.Linear(gcn_hidden, 1),
            nn.Sigmoid()
        )
        
        # Branch 2: ECFP Dense Network (Reduced to prevent overfitting)
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
        from torch_geometric.nn import global_max_pool, global_add_pool
        import torch.nn.functional as F
        
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
        graph_repr = torch.cat([sum_pool, max_pool], dim=-1)
        
        # 2. ECFP Forward
        ecfp_repr = self.ecfp_dense(ecfp)
        
        # 3. Concatenate (Early Fusion)
        joint_repr = torch.cat([graph_repr, ecfp_repr], dim=-1)
        
        # 4. Joint MLP
        return self.mlp(joint_repr)

class CGMPredictor:
    """Unified class for inferencing with ECFP, GCN, Hybrid, or Ensemble models."""
    def __init__(self, model_type='ecfp', n_tasks=105, use_gpu=False):
        self.model_type = model_type
        
        # Determine device
        if use_gpu and torch.cuda.is_available():
            self.device = torch.device('cuda')
        else:
            self.device = torch.device('cpu')
            
        if self.model_type == 'ecfp':
            self.featurizer = dc.feat.CircularFingerprint(size=2048)
            # Create dummy multitask classifier
            self.model = dc.models.MultitaskClassifier(
                n_tasks=n_tasks,
                n_features=2048,
                layer_sizes=[1024, 1024, 512],
                dropouts=[0.3, 0.3, 0.2],
                model_dir="models/deepchem_refined",
                device=self.device
            )
            self.model.restore()
        elif self.model_type == 'gcn':
            self.featurizer = dc.feat.MolGraphConvFeaturizer(use_edges=True)
            self.model = GCN(in_channels=30, hidden_channels=128, out_channels=n_tasks).to(self.device)
            self.model.load_state_dict(torch.load("models/deepchem_modern_gcn/pyg_model.pt", map_location=self.device))
            self.model.eval()
        elif self.model_type == 'hybrid':
            self.featurizer_gcn = dc.feat.MolGraphConvFeaturizer(use_edges=True)
            self.featurizer_ecfp = dc.feat.CircularFingerprint(size=2048)
            self.model = HybridModel(node_channels=30, gcn_hidden=32, ecfp_size=2048, ecfp_hidden=32, out_channels=n_tasks).to(self.device)
            self.model.load_state_dict(torch.load("models/hybrid_model/hybrid.pt", map_location=self.device))
            self.model.eval()
        elif self.model_type == 'ensemble':
            self.model_ecfp = CGMPredictor(model_type='ecfp', n_tasks=n_tasks, use_gpu=use_gpu)
            self.model_gcn = CGMPredictor(model_type='gcn', n_tasks=n_tasks, use_gpu=use_gpu)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")

    def predict(self, smiles):
        """Predicts bioactivity probabilities for a given SMILES string."""
        if self.model_type == 'ensemble':
            probs_ecfp = self.model_ecfp.predict(smiles)
            probs_gcn = self.model_gcn.predict(smiles)
            return (probs_ecfp + probs_gcn) / 2.0

        if self.model_type == 'ecfp':
            features = self.featurizer.featurize([smiles])
            if features.size == 0 or features[0].size == 0:
                raise ValueError("Featurization failed.")
            dummy_labels = np.zeros((1, self.model.n_tasks))
            dataset = dc.data.NumpyDataset(X=features, y=dummy_labels)
            preds = self.model.predict(dataset)
            return preds[0, :, 1] # Probability of active class
            
        elif self.model_type == 'gcn':
            features = self.featurizer.featurize([smiles])
            feat = features[0]
            if not hasattr(feat, 'node_features') or feat.node_features is None:
                raise ValueError("Graph featurization failed.")
            
            x = torch.tensor(feat.node_features, dtype=torch.float)
            if hasattr(feat, 'edge_index') and feat.edge_index is not None:
                edge_index = torch.tensor(feat.edge_index, dtype=torch.long)
            else:
                edge_index = torch.empty((2, 0), dtype=torch.long)
                
            batch = torch.zeros(x.size(0), dtype=torch.long)
            
            with torch.no_grad():
                out = self.model(x, edge_index, batch)
                probs = torch.sigmoid(out).numpy()[0]
            return probs
            
        elif self.model_type == 'hybrid':
            feat_gcn = self.featurizer_gcn.featurize([smiles])[0]
            feat_ecfp = self.featurizer_ecfp.featurize([smiles])[0]
            
            x = torch.tensor(feat_gcn.node_features, dtype=torch.float)
            if hasattr(feat_gcn, 'edge_index') and feat_gcn.edge_index is not None:
                edge_index = torch.tensor(feat_gcn.edge_index, dtype=torch.long)
            else:
                edge_index = torch.empty((2, 0), dtype=torch.long)
                
            batch = torch.zeros(x.size(0), dtype=torch.long)
            ecfp = torch.tensor(feat_ecfp, dtype=torch.float).view(1, -1)
            
            with torch.no_grad():
                out = self.model(x, edge_index, ecfp, batch)
                probs = torch.sigmoid(out).numpy()[0]
            return probs
