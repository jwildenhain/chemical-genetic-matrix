import os
import pandas as pd

files = {
    'Maybridge1000': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Maybridge1000_sdf_chemgrid.csv',
    'Lopac': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Lopac_sdf_chemgrid.csv',
    'SPECMTS3': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/SPECTMTS3_sdf_chemgrid.csv',
    'Spectrum128': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/Spectrum128_sdf_chemgrid.csv',
    'ChemGRID2011': '/mnt/zfspool/wdc1/jw/public_html/chemgrid/data/ChemGRID2011_sdf_chemgrid.csv'
}

for name, filepath in files.items():
    if os.path.exists(filepath):
        df = pd.read_csv(filepath, sep='\t', nrows=2)
        if df.shape[1] == 1:
            df = pd.read_csv(filepath, nrows=2)
        cols = list(df.columns)
        
        non_strain = {'sdf_structure', 'ID', 'mol_formular', 'mol_mass', 'product_name', 'db', 'c_log_p', 
                      'PM_tpsa', 'parent_smiles', 'code', 'cid', 'xsd', 'xmean', 'xmedian', 'xmin', 'xmax', 
                      'hg', 'active', 'resistent', 'nonactive', 'screened', 'ratio'}
        strains = [c for c in cols if c not in non_strain and not c.startswith('BCUT') and not c.startswith('PPSA') and not c.startswith('PNSA')]
        
        print(f"\nFile: {name:15s} | Total Columns: {len(cols)}")
        print(f"  Strain/Screen Columns Count: {len(strains)}")
        print(f"  Sample Columns: {strains[:15]}")
