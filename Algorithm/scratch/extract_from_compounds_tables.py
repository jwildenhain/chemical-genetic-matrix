import os
import re
import numpy as np
import pandas as pd
from collections import defaultdict
from rdkit import Chem

sdf_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/qnap_quake_mv_over/projects/cgm/sdf_files_for_queries/cgm_ChemGRID2011.sdf"
sql_path = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

target_libs = {"Spectrum", "SPECMTS3", "Maybridge1000", "Lopac"}

print("=================================================================")
print("EXTRACTING BIOACTIVITY DATA FROM COMPONENT & SCREEN TABLES")
print("=================================================================")

# Load SDF SMILES
cmpd_smiles = {}
cmpd_libs = {}

supplier = Chem.SDMolSupplier(sdf_path)
for mol in supplier:
    if mol is None:
        continue
    code = mol.GetProp('code') if mol.HasProp('code') else (mol.GetProp('id') if mol.HasProp('id') else None)
    smiles = mol.GetProp('parent_smiles') if mol.HasProp('parent_smiles') else Chem.MolToSmiles(mol)
    db = mol.GetProp('db') if mol.HasProp('db') else ''
    if code and smiles:
        c_code = code.strip()
        cmpd_smiles[c_code] = smiles
        cmpd_libs[c_code] = db

print(f"Loaded {len(cmpd_smiles)} compound SMILES from SDF.")

# Extract SPECMTS3 matrix from compounds_SPECMTS3_all
specmts3_cols = []
specmts3_matrix = {}

with open(sql_path, 'r', encoding='latin1') as f:
    for line in f:
        if line.startswith("CREATE TABLE `compounds_SPECMTS3_all`"):
            # read column names
            for l in f:
                l_str = l.strip()
                if l_str.startswith(");"):
                    break
                m = re.search(r"`([^`]+)`", l_str)
                if m:
                    col_name = m.group(1)
                    if col_name != 'supplier_obj_id':
                        specmts3_cols.append(col_name)
        elif line.startswith("INSERT INTO `compounds_SPECMTS3_all`"):
            idx = line.find("VALUES ")
            if idx != -1:
                content = line[idx+7:]
                items = content.split("),(")
                for item in items:
                    item = item.lstrip("(").rstrip(");\n")
                    parts = [p.strip("' ") for p in item.split(",")]
                    if len(parts) >= 2:
                        cid = parts[0]
                        vals = [float(v) if v != 'NULL' else np.nan for v in parts[1:]]
                        if cid in cmpd_smiles or ("SPE" + cid) in cmpd_smiles:
                            matched_cid = cid if cid in cmpd_smiles else ("SPE" + cid)
                            specmts3_matrix[matched_cid] = vals

print(f"Parsed SPECMTS3 Matrix: {len(specmts3_cols)} strain columns across {len(specmts3_matrix)} compounds!")
print(f"Sample SPECMTS3 Strain Columns: {specmts3_cols[:10]}")

