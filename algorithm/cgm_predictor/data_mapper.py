import pandas as pd
import os

def get_task_names(dataset_path: str = "data/predictive_matrix_cleaned.csv") -> list:
    """
    Reads the header of the dataset to get the exact list of 105 task names
    in the exact order they were used for training.
    """
    if not os.path.exists(dataset_path):
        raise FileNotFoundError(f"Dataset {dataset_path} not found. Ensure you are running from the project root.")
    
    # Read only the first row (header)
    df = pd.read_csv(dataset_path, nrows=0)
    tasks = [col for col in df.columns if col not in ['compound_id', 'smiles']]
    return tasks

def get_orf_mapping(mapping_path: str = "data/strain_orf_mapping.csv") -> dict:
    """
    Reads the mapping file and returns a dictionary mapping Strain -> ORF.
    """
    if not os.path.exists(mapping_path):
        return {} # Return empty mapping if file is missing
        
    df = pd.read_csv(mapping_path)
    mapping = dict(zip(df['Strain'], df['ORF']))
    return mapping

def map_tasks_to_orfs(tasks: list, mapping: dict) -> list:
    """
    Given a list of tasks (strains), returns a list of corresponding ORF names.
    Falls back to the original strain name if no ORF mapping is found.
    """
    return [mapping.get(task, task) for task in tasks]
