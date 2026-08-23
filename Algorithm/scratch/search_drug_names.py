import mysql.connector
import sys

def search_drugs():
    try:
        conn = mysql.connector.connect(
            host="localhost",
            user="root",
            password="password",
        )
        cursor = conn.cursor()
    except Exception as e:
        print(f"Error connecting: {e}")
        return

    cursor.execute("SHOW DATABASES;")
    databases = [d[0] for d in cursor.fetchall() if d[0] not in ('information_schema', 'mysql', 'performance_schema', 'sys')]

    drugs_to_find = ['%cyclosporin%', '%mebendazole%', '%fluconazole%', '%chrysin%']
    
    for db in databases:
        try:
            cursor.execute(f"USE {db};")
            cursor.execute("SHOW TABLES;")
            tables = [t[0] for t in cursor.fetchall()]
            
            for table in tables:
                cursor.execute(f"DESCRIBE {table};")
                columns = [c[0].lower() for c in cursor.fetchall()]
                
                # Look for name-like columns
                name_cols = [c for c in columns if 'name' in c or 'synonym' in c or 'common' in c or c == 'id']
                if not name_cols: continue
                
                for drug in drugs_to_find:
                    for col in name_cols:
                        try:
                            # Also pull smiles and compound_id if they exist
                            sel_cols = [col]
                            if 'smiles' in columns: sel_cols.append('smiles')
                            if 'compound_id' in columns: sel_cols.append('compound_id')
                            if 'id' in columns and 'id' not in sel_cols: sel_cols.append('id')
                            
                            q = f"SELECT {', '.join(sel_cols)} FROM {table} WHERE {col} LIKE '{drug}' LIMIT 5;"
                            cursor.execute(q)
                            results = cursor.fetchall()
                            for r in results:
                                print(f"FOUND {drug.replace('%', '')} in {db}.{table}: {r}")
                        except Exception:
                            pass
        except Exception as e:
            pass

    cursor.close()
    conn.close()

if __name__ == "__main__":
    search_drugs()
