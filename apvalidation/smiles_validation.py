from cmath import asin
from rdkit import Chem
from rdkit.Chem import Draw
from .smiles_to_inchi import smiles_to_inchi


def validate_smiles(smiles, img_path, asInchi=False):
    """
    Check to see if the smiles string is valid by attempting to convert it to a mol object.
    If the smiles is valid then create a structure image and return the filepath to the image along
    with the submitted smiles string. If the smiles is not valid then return the filepath to an error image
    and an empty string.

    :param smiles: a smiles string
    :return: a string and a filepath
    """

    smiles = "".join(smiles.split())
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and asInchi == False:
        generated_img_path = f"{img_path}\{smiles}.png"
        Draw.MolToFile(mol, generated_img_path, size=(200, 200), fitImage=True)
        return smiles, generated_img_path
    elif mol is not None and asInchi == True:
        inchi = smiles_to_inchi(smiles)
        generated_img_path = f"{img_path}\{inchi}.png"
        Draw.MolToFile(mol, generated_img_path, size=(200, 200), fitImage=True)
        return smiles, generated_img_path
    else:
        generated_img_path = "{img_path}/Invalid Smiles.png"
        return "", generated_img_path
