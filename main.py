from ast import Assert
from apvalidation.smiles_validation import validate_struct
from apvalidation import smiles_to_inchikey
from apvalidation.extract import Varian, Bruker, Jcampdx
import os
import nmrglue as ng
import numpy as np



"""
VARIAN TESTS
"""
for filename in os.listdir("./test_files/Lobosamide C"):
    print(filename)
    varian_dict = Varian.read([f"test_files/Lobosamide C/{filename}/procpar"])
    output = Varian.find_params(varian_dict)
    print(output)
print("\n2nd folder ---------------------------------------------------------------------------------------\n")
for filename in os.listdir("./test_files/Borrelidin analog"):
    print(filename)
    varian_dict = Varian.read([f"test_files/Borrelidin analog/{filename}/procpar"])
    output = Varian.find_params(varian_dict)
    print(output)



# dict, data = ng.varian.read_fid("./test_files/Lobosamide C/Lobos9512D_A_2_gHMBCAD_DMSO_600MHz_5mmShigemi.fid/fid 2")
# data_1D = np.array(data)
# print(data_1D)
# print(data_1D.shape)

# dict, data = ng.varian.read_fid("./test_files/Lobosamide C/Lobos9512D_A_2_gHMBCAD_DMSO_600MHz_5mmShigemi.fid/fid")
# data_2D = np.array(data)
# print(data_2D)
# print(data_2D.shape)

# check = data_1D==data_2D
# print(False in check)

"""
BRUKER TESTS
"""
# for filename in os.listdir("./test_files/Granaticin_D"):
#     print(filename)
#     if os.path.exists(f"test_files/Granaticin_D/{filename}/acqu2"):
#         file_list = [f"test_files/Granaticin_D/{filename}/acqu", f"test_files/Granaticin_D/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Granaticin_D/{filename}/acqu"]
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)

# print("\n2nd folder ---------------------------------------------------------------------------------------\n")

# for filename in os.listdir("./test_files/Demethoxy-cornuside"):
#     print(filename)
#     if os.path.exists(f"test_files/Demethoxy-cornuside/{filename}/acqu2"):
#         file_list = [f"test_files/Demethoxy-cornuside/{filename}/acqu", f"test_files/Demethoxy-cornuside/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Demethoxy-cornuside/{filename}/acqu"]
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)
