from rdkit import Chem
from rdkit.Chem import Draw
from apvalidation.smiles_to_inchikey import to_inchikey


def validate_struct(smiles, img_path, asInchiKey=False):
    """
    Check to see if the smiles string is valid by attempting to convert it to a mol object.
    If the smiles is valid then create a structure image and return the filepath to the image along
    with the submitted smiles string. If the smiles is not valid then return the filepath to an error image
    and an empty string.

    :param smiles: a smiles string to validate
    :param img_path: the file path which the output image is sent to
    :param asInchiKey: save the image name as the InchiKey rather than the Smiles string
    :return: a string and a filepath
    """

    smiles = "".join(smiles.split())
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and asInchiKey == False:
        generated_img_path = f"{img_path}/{smiles}.png"
        Draw.MolToFile(mol, generated_img_path, size=(200, 200), fitImage=True)
        return smiles, generated_img_path
    elif mol is not None and asInchiKey == True:
        inchikey = to_inchikey(smiles)
        generated_img_path = f"{img_path}/{inchikey}.png"
        Draw.MolToFile(mol, generated_img_path, size=(200, 200), fitImage=True)
        return smiles, generated_img_path
    else:
        generated_img_path = "{img_path}/Invalid Smiles.png"
        return "", generated_img_path
