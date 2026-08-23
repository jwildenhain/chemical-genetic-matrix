import os
import re
import pandas as pd
import numpy as np

print("=================================================================")
print("PARSING MAY_INFO MAYBRIDGE RECORDS & PERFORMING EXACT TABLE JOIN")
print("=================================================================")

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"
print(f"Reading SQL file: {sql_file}")

may_records = []

with open(sql_file, "r", errors="ignore") as f:
    for line in f:
        if "INSERT INTO `may_info`" in line or "INSERT INTO may_info" in line or "MSH-" in line:
            # Find all tuples starting with ('
            tuples = re.findall(r"\('MSH-[^)]+\)", line)
            for t in tuples:
                # Remove leading (' and trailing )
                content = t[2:-1]
                fields = [f.strip("'\" ") for f in content.split("','")]
                if len(fields) >= 15:
                    msh_id = fields[0]
                    supplier_obj_id = fields[12] if len(fields) > 12 else ''
                    container_id = fields[16] if len(fields) > 16 else ''
                    drug_name = fields[22] if len(fields) > 22 else ''
                    db = fields[28] if len(fields) > 28 else 'Maybridge'
                    
                    may_records.append({
                        'msh_id': msh_id,
                        'supplier_obj_id': supplier_obj_id,
                        'container_id': container_id,
                        'drug_name': drug_name,
                        'db': db
                    })

df_may = pd.DataFrame(may_records)
print(f"Total Parsed Maybridge Compound Records from `may_info`: {len(df_may)}")

if not df_may.empty:
    print("\nSample `may_info` Maybridge Records:")
    print(df_may.head(15).to_string())

# Load Maybridge1000 screened compounds (> 4 strains screened)
may1000_csv = "/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Maybridge1000_sdf_chemgrid.csv"
if os.path.exists(may1000_csv):
    df_may1000 = pd.read_csv(may1000_csv, sep='\t', on_bad_lines='skip', low_memory=False)
    df_may1000['screened'] = pd.to_numeric(df_may1000['screened'], errors='coerce')
    df_may1000_gt4 = df_may1000[df_may1000['screened'] > 4].copy()
    print(f"\nMaybridge1000 compounds with > 4 strains screened: {len(df_may1000_gt4)}")
    
    # Match with supplier_obj_id
    if not df_may.empty:
        may_suppliers = set(df_may['supplier_obj_id'].str.strip())
        df_may1000_gt4['code_clean'] = df_may1000_gt4['code'].astype(str).str.strip()
        df_may1000_gt4['ID_clean'] = df_may1000_gt4['ID'].astype(str).str.strip()
        
        joined = df_may1000_gt4[(df_may1000_gt4['code_clean'].isin(may_suppliers)) | (df_may1000_gt4['ID_clean'].isin(may_suppliers))]
        print(f"Joined Maybridge Bioactive/Cytotoxic compounds (>4 strains screened): {len(joined)}")
