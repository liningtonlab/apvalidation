from apvalidation.smiles_validation import validate_struct

print(type(validate_struct(smiles="CC(=O)NCCC1=CNc2c1cc(OC)cc2")))