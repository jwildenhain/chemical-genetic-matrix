import pandas as pd
import json
import os

def compute_frequencies():
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
    dataset_path = os.path.join(data_dir, "predictive_matrix_cleaned.csv")
    
    print(f"Loading empirical dataset: {dataset_path}")
    df = pd.read_csv(dataset_path)
    
    # Exclude non-strain columns
    exclude = ['compound_id', 'smiles']
    strains = [col for col in df.columns if col not in exclude]
    
    # Calculate baseline frequencies (mean of binary targets)
    frequencies = df[strains].mean().to_dict()
    
    # Save to JSON
    out_path = os.path.join(data_dir, "strain_frequencies.json")
    with open(out_path, 'w') as f:
        json.dump(frequencies, f, indent=4)
        
    print(f"Calculated baseline frequencies for {len(strains)} strains.")
    print(f"Saved to {out_path}")

if __name__ == "__main__":
    compute_frequencies()
