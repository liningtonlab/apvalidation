from rdkit import Chem
from rdkit.Chem import Draw
from apvalidation.smiles_to_inchikey import to_inchikey
from apvalidation.chem import convert
import io
from PIL import Image

def validate_struct(smiles):
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
    if mol is not None:
        svg_text = convert(mol)
        img = Draw.MolToImage(mol)
        output = io.BytesIO()
        img.save(output, format="PNG")
        img_text = output.getvalue()
        return img_text ,output, output.read()
    else:
        # return the svg text of an error image

        return ""

# def convert_structure(inp: str, fmt: Format = Format.sdf, get3d: bool = False):
#     try:
#         out = chem.convert(structure=inp, fmt=fmt, get3d=get3d)
#     except:
#         raise APIException(400, detail="Structure could not be converted")
#     if fmt == Format.sdf:
#         return StreamingResponse(
#             io.BytesIO(out.encode()), media_type="chemical/x-mdl-sdfile"
#         )
#     if fmt == Format.svg:
#         return StreamingResponse(io.BytesIO(out.encode()), media_type="image/svg+xml")
#     return