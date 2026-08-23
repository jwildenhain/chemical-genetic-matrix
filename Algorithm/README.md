# CGM ML Model: Yeast Bioactivity Predictor

This repository contains machine learning models designed to predict chemical-genetic interactions (bioactivity) in *Saccharomyces cerevisiae* (yeast) mutant strains based on compound structures.

By training on historical chemical-genomics screening data, the model can take a novel molecule (represented as a SMILES string or InChIKey) and predict its probability of exhibiting bioactivity across 105 distinct yeast deletion strains (ORFs).

## Features

- **Multiple Model Architectures**:
  - `ecfp`: A Dense Neural Network trained on Extended-Connectivity Fingerprints (ECFPs). Highly robust baseline for small datasets.
  - `gcn`: A Graph Convolutional Network (PyTorch Geometric) that learns directly from 2D molecular topologies.
  - `hybrid`: An Early Fusion model that concatenates GCN node features and ECFP embeddings.
  - `ensemble`: A Late Fusion approach that averages the probabilities of the ECFP and GCN models for maximum predictive power.
- **InChIKey Resolution**: Automatically queries PubChem PUG REST API if provided with an InChIKey rather than SMILES.
- **Strain Mapping**: Outputs predictions mapped to human-readable yeast strain ORF names.

## Environment Setup

We recommend using Anaconda or Miniconda. A complete conda `environment.yml` is provided, as well as a `requirements.txt` for pip users.

```bash
# Using Conda
conda env create -f environment.yml
conda activate deepchem_modern

# Using Pip
pip install -r requirements.txt
```

## Usage

Use the `predict.py` command-line interface to predict bioactivities. 

### Basic Prediction
Provide a SMILES string:
```bash
python predict.py --smiles "O=C1NC(=O)C(=C1c2ccccc2)c3ccccc3" --model ecfp
```

### Predict via InChIKey
If you only have an InChIKey, the script will automatically resolve it to SMILES:
```bash
python predict.py --inchikey "QGUGTBBHSHDMMO-UHFFFAOYSA-N" --model ensemble
```

### Adjust Threshold
By default, the script only prints strains where the probability is >= 0.5 (50%). You can adjust this:
```bash
# See all 105 strain predictions
python predict.py --smiles "O=C1NC(=O)C(=C1c2ccccc2)c3ccccc3" --threshold 0.0
```

## Repository Structure

- `predict.py`: Main CLI tool for running predictions.
- `cgm_predictor/`: Python package containing the model classes (`models.py`), PubChem API integration (`pubchem.py`), and strain mapping logic (`data_mapper.py`).
- `data/`: Cleaned dataset matrix used for training and mapping files.
- `models/`: Serialized model weights (DeepChem formats and PyTorch `.pt` files).
- `src/`: Training, evaluation, and data extraction scripts.
  - `train_pyg_model.py`: Script used to train the native GCN model.
  - `train_hybrid_model.py`: Script used to train the Early Fusion model.
  - `evaluate_ensemble.py`: Validation script for testing late fusion performance.
- `archive/`: Old iterations of models and standalone scripts kept for historical reference.

## Performance Note

Due to the limited dataset size (~2,550 unique molecules), **the DeepChem ECFP (`--model ecfp`) or the Ensemble (`--model ensemble`) architectures are recommended** for the best generalization to novel structures (Test ROC-AUC ~ 0.74 - 0.77). Early fusion networks (`--model hybrid`) have a tendency to overfit this specific dataset.
