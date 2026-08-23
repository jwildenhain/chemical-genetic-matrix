import os
import re
import gzip
import pandas as pd

print("=================================================================")
print("INSPECTING MAY_INFO & MAYBRIDGE TABLE JOINS IN MYSQL DUMP")
print("=================================================================")

sql_file = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"
dump_file = "/mnt/zfspool/Backups/Jan/2019_ugreen/db_move2nas/daily/daily_all-databases_2014-08-04_16h52m_Monday.sql.gz"

if not os.path.exists(sql_file) and os.path.exists(dump_file):
    sql_file = dump_file

print(f"Reading SQL Dump: {sql_file}")

# Open SQL file (gzip or text)
def open_sql(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "r", errors="ignore")

may_info_lines = []
may_keys_lines = []
may_db_key_lines = []

current_table = None

with open_sql(sql_file) as f:
    for line in f:
        if "CREATE TABLE" in line:
            tokens = line.split()
            if len(tokens) >= 3:
                current_table = tokens[2].strip("`(")
        elif "INSERT INTO" in line and current_table:
            if "may_info" in current_table.lower():
                may_info_lines.append(line)
            elif "may_keys" in current_table.lower():
                may_keys_lines.append(line)
            elif "may_db_key" in current_table.lower():
                may_db_key_lines.append(line)

print(f"\nExtracted table INSERT statements:")
print(f"  - may_info INSERTs:   {len(may_info_lines)}")
print(f"  - may_keys INSERTs:   {len(may_keys_lines)}")
print(f"  - may_db_key INSERTs: {len(may_db_key_lines)}")

# Inspect may_info content
print("\nSample lines from may_info:")
for l in may_info_lines[:10]:
    print(" ", l[:200])

print("\nSample lines from may_keys:")
for l in may_keys_lines[:10]:
    print(" ", l[:200])
