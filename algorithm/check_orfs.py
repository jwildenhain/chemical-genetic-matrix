import pandas as pd
import gzip
import re

def check_orfs():
    df = pd.read_csv("data/multitask_dataset_filtered_binary.csv")
    strains = [c for c in df.columns if c not in ('compound_id', 'smiles')]
    
    print(f"Total strains in our dataset: {len(strains)}")
    
    sql_path = "/mnt/zfspool/Backups/Jan/2019_ugreen/db_move2nas/daily/Maybridge/daily_Maybridge_2014-08-28_03h06m_Thursday_1QsrqF6g.sql.gz"
    
    # Extract yeast_gene_info mapping
    sym_to_orf = {}
    orf_to_sym = {}
    
    print("Parsing yeast_gene_info...")
    with gzip.open(sql_path, 'rt', encoding='latin1') as f:
        in_insert = False
        for line in f:
            if line.startswith("INSERT INTO `yeast_gene_info`"):
                idx = line.find("VALUES ")
                content = line[idx+7:]
                items = content.split("),(")
                for item in items:
                    item = item.lstrip("(").rstrip(");\n")
                    # Fields: sym, synonyms, description, function, info, orf, sgd
                    import csv
                    from io import StringIO
                    reader = csv.reader(StringIO(item), quotechar="'", escapechar='\\')
                    try:
                        parsed = next(reader)
                        if len(parsed) >= 6:
                            sym = parsed[0].strip()
                            orf = parsed[5].strip()
                            sym_to_orf[sym.lower()] = orf
                            orf_to_sym[orf.lower()] = sym
                    except Exception:
                        pass

    print(f"Loaded {len(sym_to_orf)} gene info records.")
    
    matched = []
    unmatched = []
    
    # Sort symbols by length descending to match longest possible prefix first
    sorted_syms = sorted(sym_to_orf.keys(), key=len, reverse=True)
    
    for s in strains:
        s_clean = s.lower().replace("-", "").replace("_", "")
        if s.lower() in sym_to_orf:
            matched.append((s, sym_to_orf[s.lower()]))
        elif s_clean in sym_to_orf:
            matched.append((s, sym_to_orf[s_clean]))
        else:
            found = False
            for sym in sorted_syms:
                if s.lower().startswith(sym):
                    matched.append((s, sym_to_orf[sym]))
                    found = True
                    break
            if not found:
                unmatched.append(s)

    print(f"\nStrains mapped to an ORF: {len(matched)}")
    print(f"Strains NOT mapped to an ORF: {len(unmatched)}")
    
    # Save the mapping to CSV
    mapping_df = pd.DataFrame(matched, columns=["Strain", "ORF"])
    mapping_df.to_csv("data/strain_orf_mapping.csv", index=False)
    print("Saved mapping to data/strain_orf_mapping.csv")
        
    print("\nExamples of Unmapped Strains:")
    for s in unmatched[:20]:
        print(f"  {s}")

if __name__ == "__main__":
    check_orfs()
