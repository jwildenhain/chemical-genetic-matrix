import pandas as pd
import numpy as np

def clean_predictive_dataset():
    metrics_path = "reports/per_strain_metrics.csv"
    dataset_path = "data/multitask_dataset_filtered_binary.csv"
    output_path = "data/predictive_matrix_cleaned.csv"
    
    print(f"Loading metrics from {metrics_path}...")
    metrics_df = pd.read_csv(metrics_path)
    
    print(f"Loading dataset from {dataset_path}...")
    dataset_df = pd.read_csv(dataset_path)
    
    # We want to keep strains where test_roc_auc > 0.60
    # and we also drop strains where it was NaN (couldn't compute ROC)
    valid_metrics = metrics_df.dropna(subset=['test_roc_auc'])
    predictive_strains = valid_metrics[valid_metrics['test_roc_auc'] > 0.60]['strain_task'].tolist()
    
    print(f"Total strains evaluated: {len(metrics_df)}")
    print(f"Strains with predictive power (ROC-AUC > 0.60): {len(predictive_strains)}")
    
    if len(predictive_strains) == 0:
        print("Warning: No strains met the > 0.60 threshold. Falling back to > 0.55.")
        predictive_strains = valid_metrics[valid_metrics['test_roc_auc'] > 0.55]['strain_task'].tolist()
        print(f"Strains with predictive power (ROC-AUC > 0.55): {len(predictive_strains)}")
    
    keep_cols = ['compound_id', 'smiles'] + predictive_strains
    
    cleaned_df = dataset_df[keep_cols].copy()
    
    print(f"Cleaned dataset shape: {cleaned_df.shape}")
    cleaned_df.to_csv(output_path, index=False)
    print(f"Exported cleaned predictive matrix to {output_path}")

if __name__ == "__main__":
    clean_predictive_dataset()
