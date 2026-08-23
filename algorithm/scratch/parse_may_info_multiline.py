import os
import re
import pandas as pd

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

print(f"Parsing `may_info` multi-row INSERTs from {os.path.basename(sql_file)}...")

records = []
in_may_info = False
buffer = ""

with open(sql_file, "r", errors="ignore") as f:
    for line in f:
        if "INSERT INTO `may_info`" in line or "INSERT INTO may_info" in line:
            in_may_info = True
            buffer = line
        elif in_may_info:
            buffer += line
            if ";" in line:
                in_may_info = False
                # Extract all tuples (...)
                # Find all 'supplier_obj_id' values (e.g., 'BTB 02074', 'JFD 03144', 'HTS 01153', 'MAYBRIDGE-00023')
                # supplier_obj_id is column index 12 in may_info
                matches = re.findall(r"\('(.*?)'\)", buffer)
                for m in matches:
                    parts = m.split("','")
                    if len(parts) >= 15:
                        obj_id = parts[0]
                        cat_nr = parts[2]
                        supplier_obj_id = parts[12]
                        container_id = parts[16]
                        drug_name = parts[22]
                        db = parts[28] if len(parts) > 28 else 'Maybridge'
                        records.append({
                            'obj_id': obj_id,
                            'cat_nr': cat_nr,
                            'supplier_obj_id': supplier_obj_id,
                            'container_id': container_id,
                            'drug_name': drug_name,
                            'db': db
                        })
                buffer = ""

df_may = pd.DataFrame(records)
print(f"Total Parsed `may_info` Records: {len(df_may)}")
if not df_may.empty:
    print(df_may.head(15).to_string())
