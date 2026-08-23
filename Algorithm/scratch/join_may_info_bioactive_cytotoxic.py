import os
import re
import gzip
import pandas as pd
import numpy as np
from rdkit import Chem

print("=================================================================")
print("PERFORMING MAY_INFO TABLE JOIN FOR BIOACTIVE & CYTOTOXIC LIBRARIES")
print("=================================================================")

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"
dump_file = "/mnt/zfspool/Backups/Jan/2019_ugreen/db_move2nas/daily/daily_all-databases_2014-08-04_16h52m_Monday.sql.gz"

if not os.path.exists(sql_file) and os.path.exists(dump_file):
    sql_file = dump_file

print(f"Parsing `may_info` table from {os.path.basename(sql_file)}...")

# Parse INSERT statements for may_info
may_info_records = []

def open_sql(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "r", errors="ignore")

in_may_info_data = False
with open_sql(sql_file) as f:
    for line in f:
        if "INSERT INTO `may_info`" in line or "INSERT INTO may_info" in line:
            # Extract tuples inside VALUES (...)
            values_part = line[line.find("VALUES"):].strip()
            # Match tuples
            matches = re.findall(r"\((.*?)\)[,;]", values_part)
            for m in matches:
                # Split fields handling quotes
                fields = [f.strip("'\" ") for f in m.split("','")]
                if len(fields) >= 15:
                    supplier_obj_id = fields[12] if len(fields) > 12 else ''
                    container_id = fields[16] if len(fields) > 16 else ''
                    drug_name = fields[22] if len(fields) > 22 else ''
                    db = fields[28] if len(fields) > 28 else ''
                    therap = fields[34] if len(fields) > 34 else ''
                    chem_src = fields[35] if len(fields) > 35 else ''
                    chem_class = fields[38] if len(fields) > 38 else ''
                    
                    full_text = f"{container_id} {drug_name} {db} {therap} {chem_src} {chem_class}".lower()
                    
                    is_bioactive = 'bioact' in full_text or 'bio' in full_text
                    is_cytotoxic = 'cyto' in full_text or 'toxic' in full_text or 'cytotoxic' in full_text
                    
                    may_info_records.append({
                        'supplier_obj_id': supplier_obj_id,
                        'container_id': container_id,
                        'drug_name': drug_name,
                        'db': db,
                        'therapeutic_application': therap,
                        'chemical_source': chem_src,
                        'chemical_class': chem_class,
                        'is_bioactive': is_bioactive,
                        'is_cytotoxic': is_cytotoxic
                    })

df_may_info = pd.DataFrame(may_info_records)
print(f"Total `may_info` Records Parsed: {len(df_may_info)}")

if not df_may_info.empty:
    print("Bioactive Maybridge Records:", df_may_info['is_bioactive'].sum())
    print("Cytotoxic Maybridge Records:", df_may_info['is_cytotoxic'].sum())
    print("\nSample `may_info` joined entries:\n", df_may_info[['supplier_obj_id', 'container_id', 'drug_name', 'db']].head(10).to_string())

# Load ChemGRID2011 screening data and perform join on supplier_obj_id
cg2011_csv = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/ChemGRID2011_sdf_chemgrid.csv"
if os.path.exists(cg2011_csv):
    print(f"\nJoining `may_info` with ChemGRID2011_sdf_chemgrid.csv...")
    df_cg = pd.read_csv(cg2011_csv, sep='\t', on_bad_lines='skip', low_memory=False)
    
    # Standardize compound identifier keys
    df_cg['code_clean'] = df_cg['code'].astype(str).str.strip()
    df_cg['ID_clean'] = df_cg['ID'].astype(str).str.strip()
    
    # Filter for screened > 4
    df_cg['screened'] = pd.to_numeric(df_cg['screened'], errors='coerce')
    df_cg_gt4 = df_cg[df_cg['screened'] > 4].copy()
    print(f"ChemGRID2011 compounds with screened > 4: {len(df_cg_gt4)}")
    
    if not df_may_info.empty:
        may_ids = set(df_may_info['supplier_obj_id'].astype(str).str.strip())
        joined_gt4 = df_cg_gt4[(df_cg_gt4['code_clean'].isin(may_ids)) | (df_cg_gt4['ID_clean'].isin(may_ids))]
        print(f"Joined Maybridge Bioactive/Cytotoxic Compounds (>4 strains screened): {len(joined_gt4)}")
