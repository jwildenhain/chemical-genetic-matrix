import os
import gzip
import re
import pandas as pd
from rdkit import Chem
from collections import defaultdict
import ast

def extract_dataset():
    sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
    sql_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/db_move2nas/daily/Maybridge/daily_Maybridge_2014-08-28_03h06m_Thursday_1QsrqF6g.sql.gz"
    output_dir = "data"
    os.makedirs(output_dir, exist_ok=True)
    
    target_libs = {"Maybridge", "Maybridge1000", "Bioactive", "Cytotoxic", "Spectrum_ED"}
    
    print("=== 1. Loading SMILES from SDF file ===")
    suppl = Chem.SDMolSupplier(sdf_path)
    cid_to_smiles = {}
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        try:
            smiles = Chem.MolToSmiles(mol)
            if mol.HasProp('id'):
                cid = mol.GetProp('id').strip()
                cid_to_smiles[cid] = smiles
            elif mol.HasProp('code'):
                cid = mol.GetProp('code').strip()
                cid_to_smiles[cid] = smiles
        except Exception:
            continue
    print(f"Loaded {len(cid_to_smiles)} valid compound structures.")
    
    print("\n=== 2. Parsing Z_norm from SQL dump ===")
    
    # Store: data[cid][base_strain] = [(time_snapshot_str, z_score)]
    # e.g. data['LOPAC 00581']['Alg12'] = [('-18h', -2.5), ('-20h', -3.0)]
    extracted_data = defaultdict(lambda: defaultdict(list))
    
    with gzip.open(sql_path, 'rt', encoding='latin1') as f:
        for line in f:
            if line.startswith("INSERT INTO `Z_norm`"):
                idx = line.find("VALUES ")
                if idx != -1:
                    content = line[idx + 7:]
                    items = content.split("),(")
                    for item in items:
                        item = item.lstrip("(").rstrip(");\n")
                        # Format: supplier_obj_id, strain, species, library, plate_row, plate_column, plate_number, value, outlier, z_score, p_value, ...
                        parts = item.split(",")
                        if len(parts) >= 10:
                            cid = parts[0].strip("' ")
                            strain_raw = parts[1].strip("' ")
                            lib = parts[3].strip("' ")
                            z_score_str = parts[9].strip("' ")
                            
                            if lib in target_libs and cid in cid_to_smiles:
                                try:
                                    z_score = float(z_score_str)
                                    # Extract base strain and time snapshot
                                    # Typical formats: 'Alg12-18h', 'Alg12-20h', 'WT', 'WTSpom'
                                    m = re.match(r"^(.*?)(-\d+h)?$", strain_raw)
                                    if m:
                                        base_strain = m.group(1)
                                        snapshot = m.group(2) if m.group(2) else "default"
                                        extracted_data[cid][base_strain].append((snapshot, z_score))
                                except ValueError:
                                    pass

    print(f"Parsed records for {len(extracted_data)} target compounds.")
    
    print("\n=== 3. Resolving Snapshots & Building Matrix ===")
    
    # We will choose -18h if available, otherwise the latest snapshot, or default.
    dataset_rows = []
    all_strains = set()
    
    # To properly resolve snapshots if there are a lot of actives, we need to collect all snapshots for a strain
    # The requirement: "use -18h versus later, however check the overall variability if there a lot of actives take a later snapshot"
    # To keep it simple per compound, we evaluate the best snapshot globally per strain first.
    
    # Let's aggregate all Z-scores per (base_strain, snapshot) to check variance
    strain_snapshot_stats = defaultdict(lambda: defaultdict(list))
    for cid, strains in extracted_data.items():
        for base_strain, snapshots in strains.items():
            for snap, z in snapshots:
                strain_snapshot_stats[base_strain][snap].append(z)
                
    best_snapshot_per_strain = {}
    for base_strain, snapshots in strain_snapshot_stats.items():
        if len(snapshots) == 1:
            best_snapshot_per_strain[base_strain] = list(snapshots.keys())[0]
        else:
            # Check if -18h exists
            if '-18h' in snapshots:
                # evaluate actives
                actives_18h = sum(1 for z in snapshots['-18h'] if z < -5)
                ratio_18h = actives_18h / len(snapshots['-18h'])
                
                # If too many actives (e.g. > 10%), check later snapshots
                if ratio_18h > 0.10:
                    # Find later snapshots (e.g., -20h, -22h)
                    later_snaps = [s for s in snapshots.keys() if s != '-18h' and s != 'default']
                    if later_snaps:
                        later_snaps.sort(reverse=True) # string sort works for '-20h', '-22h'
                        best_snapshot_per_strain[base_strain] = later_snaps[0]
                        print(f"Strain {base_strain}: High active ratio ({ratio_18h:.2%}) in -18h. Selected {later_snaps[0]} instead.")
                    else:
                        best_snapshot_per_strain[base_strain] = '-18h'
                else:
                    best_snapshot_per_strain[base_strain] = '-18h'
            else:
                # pick the one with most readings or first
                best_snapshot_per_strain[base_strain] = sorted(snapshots.keys())[0]

    for cid, strains in extracted_data.items():
        row = {
            'compound_id': cid,
            'smiles': cid_to_smiles[cid]
        }
        for base_strain, snapshots in strains.items():
            best_snap = best_snapshot_per_strain[base_strain]
            # Find the z-score for the best snapshot for this compound
            z_vals = [z for snap, z in snapshots if snap == best_snap]
            if z_vals:
                # If there are multiple replicates for the same snapshot, take mean/median
                z = sum(z_vals) / len(z_vals)
                row[base_strain] = z
                all_strains.add(base_strain)
        
        if len(row) > 2: # Has at least one strain reading
            dataset_rows.append(row)
            
    print(f"Final dataset covers {len(dataset_rows)} compounds and {len(all_strains)} distinct strains.")
    
    df = pd.DataFrame(dataset_rows)
    csv_path = os.path.join(output_dir, "multitask_dataset_raw.csv")
    df.to_csv(csv_path, index=False)
    print(f"Raw dataset exported to {csv_path}")

if __name__ == "__main__":
    extract_dataset()
