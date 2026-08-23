import gzip
import re
from collections import defaultdict
from rdkit import Chem

sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
sql_path = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

print("Parsing ChemGRID compound ID and strain joins...")

# 1. Map SMILES
cmpd_smiles = {}
supplier = Chem.SDMolSupplier(sdf_path)
for mol in supplier:
    if mol is None:
        continue
    code = mol.GetProp('code') if mol.HasProp('code') else (mol.GetProp('id') if mol.HasProp('id') else None)
    smiles = mol.GetProp('parent_smiles') if mol.HasProp('parent_smiles') else Chem.MolToSmiles(mol)
    if code and smiles:
        clean_code = code.strip().replace(" ", "").upper()
        cmpd_smiles[clean_code] = smiles

print(f"SDF compounds loaded: {len(cmpd_smiles)}")

# 2. Check compound tables in SQL: may_info, lopac_target_list, compounds_SPECMTS3_all
plate_well_cmpd = {}

with open(sql_path, 'r', encoding='latin1') as f:
    for line in f:
        if line.startswith("INSERT INTO `may_info`"):
            # inspect sample row in may_info
            idx = line.find("VALUES ")
            if idx != -1:
                sample = line[idx+7:idx+200]
                print("may_info sample:", sample)
        elif line.startswith("INSERT INTO `lopac_target_list`"):
            idx = line.find("VALUES ")
            if idx != -1:
                sample = line[idx+7:idx+200]
                print("lopac_target_list sample:", sample)
        elif line.startswith("INSERT INTO `compounds_SPECMTS3_all`"):
            idx = line.find("VALUES ")
            if idx != -1:
                sample = line[idx+7:idx+200]
                print("compounds_SPECMTS3_all sample:", sample)
        elif line.startswith("INSERT INTO `info_columns`"):
            idx = line.find("VALUES ")
            if idx != -1:
                sample = line[idx+7:idx+200]
                print("info_columns sample:", sample)

