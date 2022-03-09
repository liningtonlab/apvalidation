from apvalidation.smiles_validation import validate_struct
from apvalidation import smiles_to_inchikey
from apvalidation.extract import Varian, Bruker, Jcampdx
import os


for filename in os.listdir("./test_files/Lobosamide C"):
    print(filename)
    varian_dict = Varian.read(f"test_files/Lobosamide C/{filename}/procpar")
    output = Varian.find_params(varian_dict)
    print(output)


