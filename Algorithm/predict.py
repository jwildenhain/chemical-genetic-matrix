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
                        help="Raw probability threshold (default 0.5). If not binarizing, only prints strains >= threshold. Set to 0.0 to see all.")
    parser.add_argument("--min-confidence", type=float, default=1.0,
                        help="Relative Significance Score (RSS) threshold. A score > 1.0 means it is more likely than random background. (default 1.0)")
    parser.add_argument("--binarize", action="store_true",
                        help="Binarize the output to Active (1) or Inactive (0) based on thresholds.")
    parser.add_argument("--use-gpu", action="store_true",
                        help="Allow execution on GPU if available (default is CPU).")
    parser.add_argument("--show-frequent-flyers", action="store_true",
                        help="If set, includes control strains with a very high baseline hit rate (>80%).")
    
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
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        predictor = CGMPredictor(model_type=args.model, n_tasks=len(tasks), use_gpu=args.use_gpu)
    
    print(f"\nPredicting bioactivity for SMILES: {smiles}")
    try:
        if hasattr(predictor, 'predict_nuanced'):
            nuanced_scores = predictor.predict_nuanced(smiles, tasks)
        else:
            probs = predictor.predict(smiles)
            nuanced_scores = [(p, 0.0) for p in probs]
    except Exception as e:
        print(f"Error during prediction: {e}")
        sys.exit(1)
        
    print(f"\n--- PREDICTION RESULTS ---")
    
    # Bundle results
    results = []
    for i, (orf, task) in enumerate(zip(orfs, tasks)):
        prob, conf = nuanced_scores[i]
        base_freq = predictor.strain_frequencies.get(task, 0.0) if hasattr(predictor, 'strain_frequencies') else 0.0
        
        # Filter frequent flyers
        if not args.show_frequent_flyers and base_freq > 0.8:
            continue
            
        results.append((orf, task, prob, conf, base_freq))
        
    # Sort by confidence score (descending)
    results = sorted(results, key=lambda x: x[3], reverse=True)
    
    if args.binarize:
        print(f"{'ORF (Strain)':<30} | {'Status':<12} | {'Prob':<6} | {'RSS':<6}")
        print("-" * 65)
        for orf, task, prob, conf, base in results:
            display_name = f"{orf} ({task})" if orf != task else task
            status = "Active" if (prob >= args.threshold and conf >= args.min_confidence) else "Inactive"
            print(f"{display_name:<30} | {status:<12} | {prob:.4f} | {conf:.4f}")
    else:
        print(f"{'ORF (Strain)':<30} | {'Probability':<12} | {'Confidence (RSS)':<16} | {'Base Freq':<10}")
        print("-" * 80)
        count = 0
        for orf, task, prob, conf, base in results:
            if prob >= args.threshold and conf >= args.min_confidence:
                display_name = f"{orf} ({task})" if orf != task else task
                print(f"{display_name:<30} | {prob:.4f}       | {conf:.4f}           | {base:.4f}")
                count += 1
                
        if count == 0:
            print("No strains met the probability & confidence thresholds.")
            
    print("-" * 80)

if __name__ == "__main__":
    main()
