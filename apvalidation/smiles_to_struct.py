from rdkit import Chem
from rdkit.Chem import Draw


def smiles_to_struct(smiles, output_file='output.png'):
    """
    convert a smiles string into a chemical structure. Save this structure png to a file.
    :param smiles: smiles of the input compound
    :param output_file: filename of the structure png out
    :return: none
    """

    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, f'./{output_file}')
