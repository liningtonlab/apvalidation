from rdkit import Chem



def to_inchikey(smiles):
    """
    convert a smiles string to an inchi key
    :param smiles: smiles string
    :return: inchi key
    """
    if smiles == "":
        return smiles

    mol = Chem.MolFromSmiles(smiles)
    inchi_key = Chem.rdinchi.MolToInchiKey(mol)

    return inchi_key
