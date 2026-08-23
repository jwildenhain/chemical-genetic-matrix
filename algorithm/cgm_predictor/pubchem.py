import requests

def get_smiles_from_inchikey(inchikey: str) -> str:
    """
    Fetches the Canonical SMILES string for a given InChIKey using the PubChem PUG REST API.
    Raises an exception if the API request fails or the InChIKey is not found.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/CanonicalSMILES/JSON"
    response = requests.get(url, timeout=10)
    
    if response.status_code != 200:
        raise ValueError(f"Failed to fetch SMILES for InChIKey {inchikey}. HTTP Status: {response.status_code}")
        
    data = response.json()
    try:
        smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
        return smiles
    except (KeyError, IndexError):
        raise ValueError(f"Could not parse SMILES from PubChem response for InChIKey {inchikey}.")
