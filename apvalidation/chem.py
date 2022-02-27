from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDepictor import Compute2DCoords
from rdkit.Chem.Draw import rdMolDraw2D
import io
import base64

def convert(structure: str) -> str:
    """Convenience function for conversion"""
    m_canon = rdkit_atom_order(structure)
    m_canon.SetProp("_Name", AllChem.CalcMolFormula(m_canon))
    structure_svg = mol_to_svg(m_canon)
    structure_svg = io.BytesIO(structure_svg.encode())
    string_structure_svg = base64.b64encode(structure_svg.getvalue()).decode('utf8')
    return string_structure_svg


def mol_to_svg(m: Chem.Mol, size=(300, 300)) -> str:
    """Take SMILES and return SVG string"""
    d = rdMolDraw2D.MolDraw2DSVG(*size)
    d.drawOptions().addStereoAnnotation = True
    d.DrawMolecule(m)
    d.FinishDrawing()
    return d.GetDrawingText()


def rdkit_atom_order(m, add_hs=True):
    """Canonicalize using RDKit SMILES export
    Args:
        m (rdkit.Chem.Mol): Mol object for RDKit
    Returns:
        rdkit.Chem.Mol: New canonicalized RDKit mol
    """
    m_renum = Chem.MolFromSmiles(Chem.MolToSmiles(m))
    if add_hs:
        m_canon = Chem.AddHs(m_renum)
    else:
        m_canon = m_renum
    add_atom_indices(m_canon)
    return m_canon


def smi_to_mol(smi):
    return Chem.MolFromSmiles(smi)


def add_atom_indices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i + 1)
        
if __name__ == "__main__":
    convert("O=Cc1ccc(O)c(OC)c1")