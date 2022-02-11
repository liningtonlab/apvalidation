from rdkit import Chem
from rdkit.Chem import inchi


def to_inchikey(smiles):
    """
    convert a smiles string to an inchi key
    :param smiles: smiles string
    :return: inchi key
    """

    mol = Chem.MolFromSmiles(smiles)
    inchi_key = inchi.MolToInchiKey(mol)
    return inchi_key
