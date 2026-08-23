import pandas as pd
import numpy as np
import os

def filter_and_assess():
    raw_path = "data/multitask_dataset_raw.csv"
    if not os.path.exists(raw_path):
        print(f"Error: {raw_path} not found.")
        return
        
    df = pd.read_csv(raw_path)
    print(f"Loaded raw dataset with {len(df)} compounds.")
    
    strain_cols = [c for c in df.columns if c not in ('compound_id', 'smiles')]
    
    # 1. Assess Library/Strain constraints
    print(f"\n=== Strain Analysis (Total Strains: {len(strain_cols)}) ===")
    
    # Active cutoff: Z < -5 is active. Inactive > -3. 
    # But for "frequent flyer" assessment, the user said "active in more than 2/3 of cases"
    # We will define active as Z < -5 for filtering purposes.
    
    actives_mask = df[strain_cols] < -5
    
    # Calculate per-compound activity ratio
    # Number of valid (non-NaN) readings for each compound
    valid_counts = df[strain_cols].notna().sum(axis=1)
    active_counts = actives_mask.sum(axis=1)
    
    hit_ratios = active_counts / valid_counts
    
    # "exclude all molecules that are active in more than 2/3 of cases"
    # "they have to be active at least in 5% of the strains recorded against them"
    mask_valid = (hit_ratios >= 0.05) & (hit_ratios <= (2/3))
    
    df_filtered = df[mask_valid].copy()
    
    print(f"Compounds excluded (Frequent Flyers >66.7%): {(hit_ratios > (2/3)).sum()}")
    print(f"Compounds excluded (Inert <5%): {(hit_ratios < 0.05).sum()}")
    print(f"Compounds remaining: {len(df_filtered)}")
    
    print("\n=== Strain Predictive Power & Bias ===")
    # Assess variance and hit rates per strain on the filtered dataset
    strain_stats = []
    for col in strain_cols:
        col_data = df_filtered[col]
        valid_data = col_data.dropna()
        if len(valid_data) == 0:
            continue
        actives = (valid_data < -5).sum()
        inactives = (valid_data > -3).sum()
        variance = valid_data.var()
        hit_rate = actives / len(valid_data)
        
        strain_stats.append({
            'Strain': col,
            'Total_Screened': len(valid_data),
            'Actives (Z < -5)': actives,
            'Inactives (Z > -3)': inactives,
            'Borderline (-5 <= Z <= -3)': len(valid_data) - actives - inactives,
            'Hit_Rate': hit_rate,
            'Variance': variance
        })
        
    stats_df = pd.DataFrame(strain_stats).sort_values(by='Variance', ascending=False)
    print(stats_df.head(20).to_string(index=False))
    
    print("\n=== Binarizing Dataset for DeepChem ===")
    # Z < -5 -> 1 (Active)
    # Z > -3 -> 0 (Inactive)
    # Else -> NaN
    df_binary = df_filtered.copy()
    for col in strain_cols:
        # We must use float to store NaNs
        conditions = [
            df_binary[col] < -5,
            df_binary[col] > -3
        ]
        choices = [1.0, 0.0]
        df_binary[col] = np.select(conditions, choices, default=np.nan)
        
    # Export filtered dataset
    out_path = "data/multitask_dataset_filtered_binary.csv"
    df_binary.to_csv(out_path, index=False)
    print(f"\nFiltered and binarized dataset saved to {out_path}")

if __name__ == "__main__":
    filter_and_assess()
