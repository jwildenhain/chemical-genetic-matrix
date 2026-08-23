import os
import pandas as pd
from rdkit import Chem

molecules = ['cyclosporin', 'mebendazole', 'fluconazole', 'chrysin', 'ciclosporin', 'flucanozole']
sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
# Check if predicting script exists
sys_path = "/home/jw/Source/CGM_ML_Model"

found_mols = {}
supplier = Chem.SDMolSupplier(sdf_path)
for mol in supplier:
    if mol is None:
        continue
    name = mol.GetProp('_Name') if mol.HasProp('_Name') else mol.GetProp('product_name') if mol.HasProp('product_name') else ''
    
    if name:
        for m in molecules:
            if m.lower() in name.lower():
                smiles = mol.GetProp('parent_smiles') if mol.HasProp('parent_smiles') else Chem.MolToSmiles(mol)
                found_mols[m.lower()] = smiles

for k, v in found_mols.items():
    print(f"{k}: {v}")

