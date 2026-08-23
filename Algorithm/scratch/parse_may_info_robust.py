import os
import pandas as pd

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

print(f"Parsing `may_info` lines from {os.path.basename(sql_file)}...")

may_rows = []
with open(sql_file, "r", errors="ignore") as f:
    for line in f:
        if line.startswith("INSERT INTO `may_info`") or line.startswith("INSERT INTO may_info"):
            # Extract everything after VALUES
            val_str = line[line.find("VALUES") + 6:].strip().rstrip(";")
            # Split rows by '),('
            row_strs = val_str.split("),(")
            for r in row_strs:
                r_clean = r.lstrip("(").rstrip(")")
                # Split CSV tokens
                items = [it.strip("'\" ") for it in r_clean.split("', '")]
                if len(items) >= 13:
                    may_rows.append({
                        'obj_id': items[0],
                        'cat_nr': items[2] if len(items) > 2 else '',
                        'supplier_obj_id': items[12] if len(items) > 12 else '',
                        'container_id': items[16] if len(items) > 16 else '',
                        'drug_name': items[22] if len(items) > 22 else '',
                        'db': items[28] if len(items) > 28 else ''
                    })

df_may = pd.DataFrame(may_rows)
print(f"Total Parsed `may_info` Rows: {len(df_may)}")
print(df_may['db'].value_counts().to_string())
print("\nSample Maybridge Entries:\n", df_may.head(20).to_string())
