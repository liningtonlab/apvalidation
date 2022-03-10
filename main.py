from apvalidation.smiles_validation import validate_struct
from apvalidation import smiles_to_inchikey
from apvalidation.extract import Varian, Bruker, Jcampdx
import os
import nmrglue as ng
import numpy as np


# for filename in os.listdir("./test_files/Lobosamide C"):
#     print(filename)
#     varian_dict = Varian.read(f"test_files/Lobosamide C/{filename}/procpar")
#     output = Varian.find_params(varian_dict)
#     print(output)

dict, data = ng.varian.read_fid("./test_files/Lobosamide C/Lobos9512D_A_2_gHMBCAD_DMSO_600MHz_5mmShigemi.fid/fid 2")
data_1D = np.array(data)
print(data_1D)
print(data_1D.shape)

dict, data = ng.varian.read_fid("./test_files/Lobosamide C/Lobos9512D_A_2_gHMBCAD_DMSO_600MHz_5mmShigemi.fid/fid")
data_2D = np.array(data)
print(data_2D)
print(data_2D.shape)

check = data_1D==data_2D
print(False in check)
