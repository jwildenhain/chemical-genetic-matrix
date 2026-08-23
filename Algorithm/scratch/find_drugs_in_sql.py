import re

def parse_sql():
    dump_path = "/mnt/zfspool/Backups/DS712Plus/DS712PlusHomeFolder/ExamineOldDataBackups/sbsr-2059 - ok/ChemGrid.latest.sql"
    drugs = ['cyclosporin', 'mebendazole', 'fluconazole', 'chrysin']
    
    with open(dump_path, 'r', encoding='latin1') as f:
        for line in f:
            if 'INSERT INTO' not in line:
                continue
            
            # The line looks like: INSERT INTO `table_name` VALUES ...
            table_name = "UNKNOWN_TABLE"
            if "INSERT INTO" in line:
                parts = line.split(" VALUES ")
                table_name = parts[0].replace("INSERT INTO", "").strip(" `")
            
            # Split the massive line into individual tuples
            tuples = line.split("),(")
            for t in tuples:
                tl = t.lower()
                for d in drugs:
                    if d in tl:
                        print(f"--- FOUND {d.upper()} IN TABLE: {table_name} ---")
                        clean_t = t.strip("(); \n\r")
                        fields = clean_t.split("','")
                        print("FIELDS:", len(fields))
                        for i, field in enumerate(fields):
                            field = field.strip("'")
                            if len(field) > 0 and len(field) < 100:
                                print(f"  [{i}]: {field}")
                        print("Raw snippet:", clean_t[:200])

if __name__ == "__main__":
    parse_sql()
