import sys
import os

# Add the parent directory to the Python path so we can import the package
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cgm_predictor.pubchem import get_smiles_from_inchikey
from cgm_predictor.data_mapper import get_task_names, get_orf_mapping, map_tasks_to_orfs
from cgm_predictor.models import CGMPredictor

def run_demo():
    print("=== CGM Predictor programmatic demo ===\n")
    
    # 1. Provide a SMILES string (Example: Mebendazole)
    smiles = "COC(=O)NC1=NC2=C(N1)C=C(C=C2)C(=O)C3=CC=CC=C3"
    print(f"Input SMILES: {smiles}\n")
    
    # 2. Get Metadata
    tasks = get_task_names(dataset_path="../data/predictive_matrix_cleaned.csv")
    mapping = get_orf_mapping(mapping_path="../data/strain_orf_mapping.csv")
    orfs = map_tasks_to_orfs(tasks, mapping)
    print(f"Loaded {len(tasks)} target strain conditions.\n")
    
    # 3. Load the Models
    print("Initializing PyTorch Geometric GCN Model (on CPU to ensure RTX 5080 stability)...")
    
    # We suppress deepchem/tf warnings for clean output
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        predictor_gcn = CGMPredictor(model_type='gcn', n_tasks=len(tasks))
    
    # 4. Run Inference
    print("Running Inference...")
    probs = predictor_gcn.predict(smiles)
    
    # 5. Output Results
    print("\n--- TOP PREDICTION RESULTS (Probability >= 0.90) ---")
    results = sorted(zip(orfs, tasks, probs), key=lambda x: x[2], reverse=True)
    
    print(f"{'ORF (Strain)':<30} | {'Probability':<12}")
    print("-" * 45)
    for orf, task, prob in results:
        if prob >= 0.90:
            display_name = f"{orf} ({task})" if orf != task else task
            print(f"{display_name:<30} | {prob:.2%}")

if __name__ == "__main__":
    run_demo()
