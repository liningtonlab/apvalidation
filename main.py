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

# print("\nOUTPUT FOR Lobosamide C---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/Lobosamide C"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/Lobosamide C/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
print("\nOUTPUT FOR Borrelidin analog---------------------------------------------------------------------------------------\n")
for filename in os.listdir("./test_files/Borrelidin analog"):
    print(filename)
    varian_dict = Varian.read([f"test_files/Borrelidin analog/{filename}/procpar"])
    output = Varian.find_params(varian_dict)
    print(f"{output}\n")
# print("\nOUTPUT FOR 1353_Day3_A1_B---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/1353_Day3_A1_B"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/1353_Day3_A1_B/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR 1353_SYP Original Fill Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/1353_SYP Original Full Data Set"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/1353_SYP Original Full Data Set/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")



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
print("\nOUTPUT FOR Granaticin_D---------------------------------------------------------------------------------------\n")
for filename in os.listdir("./test_files/Granaticin_D"):
    print(filename)
    if os.path.exists(f"test_files/Granaticin_D/{filename}/acqu2"):
        file_list = [f"test_files/Granaticin_D/{filename}/acqu", f"test_files/Granaticin_D/{filename}/acqu2"]
    else:
        file_list = [f"test_files/Granaticin_D/{filename}/acqu"]
    dict_list = Bruker.read(file_list)
    params = Bruker.find_params(dict_list)
    print(params)

# print("\nOUTPUT FOR Demethoxy-cornuside---------------------------------------------------------------------------------------\n")

# for filename in os.listdir("./test_files/Demethoxy-cornuside"):
#     print(filename)
#     if os.path.exists(f"test_files/Demethoxy-cornuside/{filename}/acqu2"):
#         file_list = [f"test_files/Demethoxy-cornuside/{filename}/acqu", f"test_files/Demethoxy-cornuside/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Demethoxy-cornuside/{filename}/acqu"]
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)

# print("\nOUTPUT FOR SL_RLUS-2152D-1-1---------------------------------------------------------------------------------------\n")

# for filename in os.listdir("./test_files/SL_RLUS-2152D-1-1"):
#     print(filename)
#     if os.path.exists(f"test_files/SL_RLUS-2152D-1-1/{filename}/acqu2"):
#         file_list = [f"test_files/SL_RLUS-2152D-1-1/{filename}/acqu", f"test_files/SL_RLUS-2152D-1-1/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/SL_RLUS-2152D-1-1/{filename}/acqu"]
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)


"""
JEOL DATA
"""
# print("\nOUTPUT FOR RGL1617G1B proton(Jcampdx).jdx---------------------------------------------------------------------------------------\n")
# dict_list = Jcampdx.read(['test_files/RGL1617G1B proton(Jcampdx).jdx'])
# params = Jcampdx.find_params(dict_list)
# print(params)