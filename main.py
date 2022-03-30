from ast import Assert
# from apvalidation.smiles_validation import validate_struct
# from apvalidation import smiles_to_inchikey
from apvalidation.extract import Varian, Bruker, JEOL, Jcampdx_Handler
from apvalidation.smiles_validation import validate_struct
from rdkit import Chem
import os
import pandas as pd
import nmrglue as ng
import numpy as np
import json



print(validate_struct("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5", "./" ))


"""
VARIAN TESTS
"""

# print("\nOUTPUT FOR Lobosamide C---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/Lobosamide C"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/Lobosamide C/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR Borrelidin analog---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/Borrelidin analog"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/Borrelidin analog/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
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

# print("\nOUTPUT FOR Aspochalasin I Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/Aspochalasin I"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/Aspochalasin I/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR echinulin Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/echinulin"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/echinulin/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR emerimicin V Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/emerimicin V"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/emerimicin V/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR MB0593C-BAC-AA Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/MB0593C-BAC-AA"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/MB0593C-BAC-AA/{filename}/procpar"])
#     output = Varian.find_params(varian_dict)
#     print(f"{output}\n")
# print("\nOUTPUT FOR MB0593E-BAC-DD Data Set---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/MB0593E-BAC-DD"):
#     print(filename)
#     varian_dict = Varian.read([f"test_files/MB0593E-BAC-DD/{filename}/procpar"])
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
# print("\nOUTPUT FOR Granaticin_D---------------------------------------------------------------------------------------\n")
# for filename in os.listdir("./test_files/Granaticin_D"):
#     print(filename)
#     if os.path.exists(f"test_files/Granaticin_D/{filename}/acqu2"):
#         file_list = [f"test_files/Granaticin_D/{filename}/acqu", f"test_files/Granaticin_D/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Granaticin_D/{filename}/acqu"]
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)

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

# print("\nOUTPUT FOR Tetronasin---------------------------------------------------------------------------------------\n")

# for filename in os.listdir("./test_files/Tetronasin"):
#     print(filename)
#     if os.path.exists(f"test_files/Tetronasin/{filename}/acqu2"):
#         file_list = [f"test_files/Tetronasin/{filename}/acqu", f"test_files/Tetronasin/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Tetronasin/{filename}/acqu"]
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


"""
JDX testing
"""
print("\nOUTPUT Lobos9512D_A_proton_Pyridine_600MHz_5mmregtube.jdx ---------------------------------------------------------------------------------------\n")
# varian_dict = Varian.read(["test_files/Lobos9512D_A_proton_Pyridine_600MHz_5mmregtube.jdx"])
# jcamp_dict = Jcampdx_Handler.read(["test_files/Lobos9512D_A_proton_Pyridine_600MHz_5mmregtube.jdx"])
# find_params_output = Jcampdx_Handler.find_params(jcamp_dict)
# print(find_params_output)

# jcamp_dict = Jcampdx_Handler.read(["test_files/Granaticin_C_Proton.jdx"])
# find_params_output = Jcampdx_Handler.find_params(jcamp_dict)
# print(find_params_output)

print("FIRST FILE ------\n\n")
jcamp_dict = Jcampdx_Handler.read(["test_files/Bruker_HMBC.jdx"])
output = Jcampdx_Handler.find_params(jcamp_dict)


print("SECOND FILE ------\n\n")
jcamp_dict2 = Jcampdx_Handler.read(["test_files/bruker_j.jdx"])
with open("./test_1d.json", "w") as f:
    json.dump(jcamp_dict2[0][0], f)
output2 = Jcampdx_Handler.find_params(jcamp_dict2)


print("FINAL OUTPUTS!!!")
print(output)
print(output2)
# print(jcamp_dict[0][0]["_datatype_LINK"][0]['_comments'])


# with open("./dump.json", "w") as f:
#     json.dump(jcamp_dict[0][0]["_datatype_LINK"][0], f)



# find_params_output = Jcampdx_Handler.find_params(jcamp_dict)
# print("This is the output for the Bruker file.")
# print(find_params_output)
# jcamp_string = '\n'.join(list_of_data).split('\n')
# print(jcamp_string)
# with open("/workspaces/apvalidation/test_files/fake_procpar.txt", "w") as f:
#     f.write(jcamp_string)
#     f.close()
# varian_dict = Varian.read(["/workspaces/apvalidation/test_files/fake_procpar.txt"])
# print(varian_dict)


# print(dict(zip(list_of_data[::2], list_of_data[1::2])))
# output = Jcampdx_Handler.find_params(jcamp_dict)
# print(output)
# print(jcamp_dict)
# print(jcamp_dict[0][0]["_datatype_LINK"][0]["$ORIGINALFORMAT"])
# afile = open(r'/workspaces/apvalidation/json_dump', 'w', encoding='utf8')
# json.dump(jcamp_dict, afile)

# print(Chem.MolFromSmiles("C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C"))
# validate_struct(smiles="C[C@H]1C[C@H](C[C@@H]([C@H](/C(=C\C=C\C[C@H](OC(=O)C[C@@H]([C@H](C1)C)O)[C@@H]2CCC[C@H]2C(=O)O)/C#N)O)C)C", img_path=".", asInchiKey=True)





