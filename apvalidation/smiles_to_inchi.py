from rdkit import Chem
from rdkit.Chem import inchi


def smiles_to_inchi(smiles):
    """
    convert a smiles string to an inchi key
    :param smiles: smiles string
    :return: inchi key
    """

    mol = Chem.MolFromSmiles(smiles)
    inchi_key = inchi.MolToInchiKey(mol)
    return inchi_key
