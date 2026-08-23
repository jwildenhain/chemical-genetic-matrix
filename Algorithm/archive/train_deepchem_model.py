import os
import sys
import json
import torch
import numpy as np
np.object = object
np.bool = bool
np.int = int
np.typeDict = np.sctypeDict
import pandas as pd
import deepchem as dc
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, accuracy_score

print("=================================================================")
print("DEEPCHEM GPU MULTI-TASK MODEL: TARGET LIBRARIES (>4 STRAINS SCREENED, 1,554 COMPS)")
print("=================================================================")

# GPU Verification
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"PyTorch CUDA Available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Using GPU: {torch.cuda.get_device_name(0)} (VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB)")

dataset_file = "data/multitask_dataset_filtered_binary.csv"
if not os.path.exists(dataset_file):
    print(f"Error: {dataset_file} not found. Did you run filter_and_assess_bias.py?")
    sys.exit(1)
print(f"Loading dataset: {dataset_file}")

df = pd.read_csv(dataset_file)
non_task_cols = {'smiles', 'compound_id'}
tasks = [c for c in df.columns if c not in non_task_cols]

print(f"Total Filtered Compounds: {len(df)}")
print(f"Total Yeast Strain Tasks ({len(tasks)}): {tasks[:5]}...")

os.makedirs("models/deepchem_multitask", exist_ok=True)
os.makedirs("reports", exist_ok=True)

print("\n=================================================================")
print("STEP 1: Featurizing Target Compounds with ECFP4 Fingerprints")
print("=================================================================")

featurizer = dc.feat.CircularFingerprint(size=2048, radius=2)
loader = dc.data.CSVLoader(
    tasks=tasks,
    featurizer=featurizer,
    feature_field="smiles",
    id_field="smiles"
)

dataset = loader.create_dataset(dataset_file)
print(f"Featurized dataset shape: {dataset.get_shape()}")

print("\n=================================================================")
print("STEP 2: Scaffold Split (Train 80% / Val 10% / Test 10%)")
print("=================================================================")

splitter = dc.splits.ScaffoldSplitter()
train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(
    dataset,
    frac_train=0.8,
    frac_valid=0.1,
    frac_test=0.1,
    seed=42
)

print(f"Original Train Dataset size: {len(train_dataset)}")
print(f"Valid Dataset size:          {len(valid_dataset)}")
print(f"Test Dataset size:           {len(test_dataset)}")

print("\n=================================================================")
print("STEP 3: Applying Vectorized MLSMOTE Class Imbalance Resampling")
print("=================================================================")

def apply_mlsmote_vectorized(X, Y, W, ids, oversample_count=1000, k=5):
    pos_indices = np.where(Y.sum(axis=1) > 0)[0]
    print(f"Bioactive compounds in train set: {len(pos_indices)} / {len(X)}")
    
    if len(pos_indices) < k + 1:
        print("Not enough minority samples for KNN interpolation. Returning original.")
        return X, Y, W, ids
        
    print(f"Generating {oversample_count} synthetic bioactive samples via Vectorized MLSMOTE...")
    knn = NearestNeighbors(n_neighbors=k+1, metric='euclidean').fit(X[pos_indices])
    distances, indices = knn.kneighbors(X[pos_indices])
    
    rand_base_idx = np.random.choice(len(pos_indices), size=oversample_count)
    rand_neighbor_subidx = np.random.randint(1, k+1, size=oversample_count)
    
    base_x = X[pos_indices[rand_base_idx]]
    base_y = Y[pos_indices[rand_base_idx]]
    base_w = W[pos_indices[rand_base_idx]]
    
    neighbor_pos_idx = indices[rand_base_idx, rand_neighbor_subidx]
    neighbor_x = X[pos_indices[neighbor_pos_idx]]
    neighbor_y = Y[pos_indices[neighbor_pos_idx]]
    neighbor_w = W[pos_indices[neighbor_pos_idx]]
    
    alphas = np.random.uniform(0.1, 0.9, size=(oversample_count, 1))
    syn_x_cont = base_x + alphas * (neighbor_x - base_x)
    syn_x_bin = (syn_x_cont >= 0.5).astype(np.float32)
    
    syn_y = np.logical_or(base_y, neighbor_y).astype(np.float32)
    syn_w = np.maximum(base_w, neighbor_w)
    syn_ids = np.array([f"SYNTHETIC_SMOTE_{i+1}" for i in range(oversample_count)])
    
    res_X = np.vstack([X, syn_x_bin])
    res_Y = np.vstack([Y, syn_y])
    res_W = np.vstack([W, syn_w])
    res_ids = np.concatenate([ids, syn_ids])
    
    return res_X, res_Y, res_W, res_ids

X_res, Y_res, W_res, ids_res = apply_mlsmote_vectorized(
    train_dataset.X,
    train_dataset.y,
    train_dataset.w,
    train_dataset.ids,
    oversample_count=1000,
    k=5
)

train_dataset_smote = dc.data.NumpyDataset(
    X=X_res,
    y=Y_res,
    w=W_res,
    ids=ids_res
)

print(f"MLSMOTE-Augmented Train Dataset size: {len(train_dataset_smote)} (Original: {len(train_dataset)})")

print("\n=================================================================")
print("STEP 4: Initializing & Training DeepChem MultitaskClassifier on GPU")
print("=================================================================")

model = dc.models.MultitaskClassifier(
    n_tasks=len(tasks),
    n_features=2048,
    layer_sizes=[1024, 512, 256],
    dropouts=0.2,
    learning_rate=0.0005,
    model_dir="models/deepchem_multitask"
)

nb_epoch = 30
print(f"Training MultitaskClassifier on GPU for {nb_epoch} epochs...")

loss_history = []
for epoch in range(1, nb_epoch + 1):
    loss = model.fit(train_dataset_smote, nb_epoch=1)
    loss_history.append(float(loss))
    if epoch % 5 == 0 or epoch == 1:
        print(f"Epoch {epoch:2d}/{nb_epoch:2d} | Training Loss: {loss:.4f}")

print("\n=================================================================")
print("STEP 5: Evaluating Model Metrics Across Yeast Strains")
print("=================================================================")

def compute_multitask_metrics(dset):
    y_true = dset.y
    w = dset.w
    raw_pred = model.predict(dset)
    
    if isinstance(raw_pred, list):
        raw_pred = np.array(raw_pred)
        
    if raw_pred.ndim == 3 and raw_pred.shape[2] == 2:
        y_pred_prob = raw_pred[:, :, 1]
    else:
        y_pred_prob = raw_pred
        
    roc_aucs = []
    prc_aucs = []
    accs = []
    per_task_rocs = []
    
    for i in range(y_true.shape[1]):
        valid_mask = w[:, i] > 0
        yt = y_true[valid_mask, i]
        yp = y_pred_prob[valid_mask, i]
        
        if len(np.unique(yt)) > 1:
            roc = roc_auc_score(yt, yp)
            p, r, _ = precision_recall_curve(yt, yp)
            prc = auc(r, p)
            acc = accuracy_score(yt, (yp >= 0.5).astype(int))
            
            roc_aucs.append(roc)
            prc_aucs.append(prc)
            accs.append(acc)
            per_task_rocs.append(roc)
        else:
            per_task_rocs.append(np.nan)
            
    return {
        'roc_auc': float(np.nanmean(roc_aucs)) if roc_aucs else 0.0,
        'prc_auc': float(np.nanmean(prc_aucs)) if prc_aucs else 0.0,
        'accuracy': float(np.nanmean(accs)) if accs else 0.0,
        'per_task_roc_auc': per_task_rocs
    }

train_scores = compute_multitask_metrics(train_dataset_smote)
valid_scores = compute_multitask_metrics(valid_dataset)
test_scores  = compute_multitask_metrics(test_dataset)

print("\n--- OVERALL MULTI-TASK EVALUATION METRICS (FILTERED >4 STRAINS DATASET) ---")
print(f"Train Set -> ROC-AUC: {train_scores['roc_auc']:.4f} | PR-AUC: {train_scores['prc_auc']:.4f} | Acc: {train_scores['accuracy']:.4f}")
print(f"Valid Set -> ROC-AUC: {valid_scores['roc_auc']:.4f} | PR-AUC: {valid_scores['prc_auc']:.4f} | Acc: {valid_scores['accuracy']:.4f}")
print(f"Test Set  -> ROC-AUC: {test_scores['roc_auc']:.4f} | PR-AUC: {test_scores['prc_auc']:.4f} | Acc: {test_scores['accuracy']:.4f}")

per_task_results = []
for idx, task_name in enumerate(tasks):
    roc = test_scores['per_task_roc_auc'][idx]
    per_task_results.append({
        'strain_task': task_name,
        'test_roc_auc': round(float(roc), 4) if not np.isnan(roc) else None
    })

per_task_df = pd.DataFrame(per_task_results)
per_task_df.to_csv("reports/per_strain_metrics.csv", index=False)

summary_dict = {
    'device': device,
    'total_compounds': len(df),
    'num_tasks': len(tasks),
    'original_train_size': len(train_dataset),
    'mlsmote_train_size': len(train_dataset_smote),
    'train_roc_auc': round(train_scores['roc_auc'], 4),
    'valid_roc_auc': round(valid_scores['roc_auc'], 4),
    'test_roc_auc': round(test_scores['roc_auc'], 4),
    'test_prc_auc': round(test_scores['prc_auc'], 4),
    'test_accuracy': round(test_scores['accuracy'], 4),
    'final_loss': round(float(loss_history[-1]), 4)
}

with open("reports/model_summary.json", "w") as sf:
    json.dump(summary_dict, sf, indent=2)

print("\nFiltered >4 strains dataset training & evaluation completed successfully!")
