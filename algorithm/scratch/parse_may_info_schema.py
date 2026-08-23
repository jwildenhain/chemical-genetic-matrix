import os
import re
import pandas as pd

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"

print(f"Reading schema for may_info in {os.path.basename(sql_file)}...")

with open(sql_file, 'r', errors='ignore') as f:
    create_stmt = []
    in_may_info = False
    for line in f:
        if "CREATE TABLE `may_info`" in line or "CREATE TABLE may_info" in line:
            in_may_info = True
            create_stmt.append(line)
        elif in_may_info:
            create_stmt.append(line)
            if ";" in line:
                break

print("`may_info` CREATE TABLE Schema:\n")
print("".join(create_stmt))
