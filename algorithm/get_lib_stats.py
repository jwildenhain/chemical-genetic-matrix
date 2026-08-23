import gzip
from collections import defaultdict
import time

def get_stats():
    sql_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/db_move2nas/daily/Maybridge/daily_Maybridge_2014-08-28_03h06m_Thursday_1QsrqF6g.sql.gz"
    
    # Store sets of strains and compounds for each library
    lib_strains = defaultdict(set)
    lib_compounds = defaultdict(set)
    
    print("Parsing Z_norm table...")
    start_time = time.time()
    
    with gzip.open(sql_path, 'rt', encoding='latin1') as f:
        for line in f:
            if line.startswith("INSERT INTO `Z_norm`"):
                idx = line.find("VALUES ")
                if idx != -1:
                    content = line[idx + 7:]
                    items = content.split("),(")
                    for item in items:
                        item = item.lstrip("(").rstrip(");\n")
                        parts = item.split(",")
                        if len(parts) >= 10:
                            cid = parts[0].strip("' ")
                            strain = parts[1].strip("' ")
                            lib = parts[3].strip("' ")
                            
                            lib_strains[lib].add(strain)
                            lib_compounds[lib].add(cid)

    print(f"Parsing took {time.time() - start_time:.2f} seconds.\n")
    
    print(f"{'Library Name':<25} | {'Unique Strains':<15} | {'Unique Compounds':<18}")
    print("-" * 65)
    for lib in sorted(lib_strains.keys()):
        num_strains = len(lib_strains[lib])
        num_comps = len(lib_compounds[lib])
        print(f"{lib:<25} | {num_strains:<15} | {num_comps:<18}")

if __name__ == "__main__":
    get_stats()
