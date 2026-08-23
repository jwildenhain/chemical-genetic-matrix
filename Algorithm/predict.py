import argparse
import sys
from cgm_predictor.pubchem import get_smiles_from_inchikey
from cgm_predictor.data_mapper import get_task_names, get_orf_mapping, map_tasks_to_orfs
from cgm_predictor.models import CGMPredictor

def main():
    parser = argparse.ArgumentParser(description="Predict yeast mutant strain bioactivity using CGM DeepChem/PyG Models.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--smiles", type=str, help="SMILES string of the molecule.")
    group.add_argument("--inchikey", type=str, help="InChIKey of the molecule (will be resolved to SMILES via PubChem API).")
    
    parser.add_argument("--model", type=str, choices=['ecfp', 'gcn', 'hybrid', 'ensemble'], default='ecfp',
                        help="Which model to use: 'ecfp' (DeepChem native), 'gcn' (PyG GCN), 'hybrid' (Early Fusion GCN+ECFP), or 'ensemble' (Late Fusion Average). Default is ecfp.")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Threshold for considering a strain 'Active' (default 0.5). If not binarizing, only prints strains >= threshold. Set to 0.0 to see all.")
    parser.add_argument("--binarize", action="store_true",
                        help="Binarize the output to Active (1) or Inactive (0) based on the threshold. Displays all strains.")
    parser.add_argument("--use-gpu", action="store_true",
                        help="Allow execution on GPU if available (default is CPU).")
    
    args = parser.parse_args()

    smiles = args.smiles
    if args.inchikey:
        print(f"Resolving InChIKey {args.inchikey} to SMILES via PubChem...")
        try:
            smiles = get_smiles_from_inchikey(args.inchikey)
            print(f"Resolved SMILES: {smiles}")
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

    print(f"\nLoading metadata...")
    tasks = get_task_names()
    mapping = get_orf_mapping()
    orfs = map_tasks_to_orfs(tasks, mapping)
    
    print(f"Loading {args.model.upper()} model...")
    # Suppress deepchem TF warnings if possible
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        predictor = CGMPredictor(model_type=args.model, n_tasks=len(tasks), use_gpu=args.use_gpu)
    
    print(f"\nPredicting bioactivity for SMILES: {smiles}")
    try:
        probs = predictor.predict(smiles)
    except Exception as e:
        print(f"Error during prediction: {e}")
        sys.exit(1)
        
    print(f"\n--- PREDICTION RESULTS (Threshold >= {args.threshold}) ---")
    
    results = sorted(zip(orfs, tasks, probs), key=lambda x: x[2], reverse=True)
    
    if args.binarize:
        print(f"{'ORF (Strain)':<30} | {'Status':<12}")
        print("-" * 45)
        for orf, task, prob in results:
            display_name = f"{orf} ({task})" if orf != task else task
            status = "Active" if prob >= args.threshold else "Inactive"
            print(f"{display_name:<30} | {status:<12}")
    else:
        print(f"{'ORF (Strain)':<30} | {'Probability':<12}")
        print("-" * 45)
        count = 0
        for orf, task, prob in results:
            if prob >= args.threshold:
                display_name = f"{orf} ({task})" if orf != task else task
                print(f"{display_name:<30} | {prob:.4f}")
                count += 1
                
        if count == 0:
            print("No strains met the probability threshold.")
            
    print("-" * 45)

if __name__ == "__main__":
    main()
