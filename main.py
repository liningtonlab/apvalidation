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
import re
import nmrglue as ng



# print(validate_struct("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5", "./" ))


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
# print("REFERENCE FILES OUTPUT--------------------")
# for filename in os.listdir("test_files/Granaticin_C"):
#     print(filename)
#     if filename == '.DS_Store':
#         continue
#     if os.path.exists(f"test_files/Granaticin_C/{filename}/acqu2"):
#         file_list = [f"test_files/Granaticin_C/{filename}/acqu", f"test_files/Granaticin_C/{filename}/acqu2"]
#     else:
#         file_list = [f"test_files/Granaticin_C/{filename}/acqu"]
#     print(file_list)
#     dict_list = Bruker.read(file_list)
#     params = Bruker.find_params(dict_list)
#     print(params)



# print("JDX FILES OUTPUT--------------------------")
# print("VARIAN 2D:")
# jcamp_dict_varian2d = Jcampdx_Handler.read(["test_files/Lobos9512D_A_proton_Pyridine_600MHz_5mmregtube.jdx"])

# print(Jcampdx_Handler.find_params(jcamp_dict_varian2d))

# print("BRUKER #################################################################################################################################################")
# print("BRUKER 1D--------------------------")
jcamp_dict_bruker1d = Jcampdx_Handler.read(["test_files/bruker_j.jdx"])
print(f"Length of the dict.keys() is: {jcamp_dict_bruker1d[0].keys()} ")
print(Jcampdx_Handler.find_params(jcamp_dict_bruker1d))

print("BRUKER 2D--------------------------")
jcamp_dict_bruker2d = Jcampdx_Handler.read(["test_files/Bruker_HMBC.jdx"])
print(Jcampdx_Handler.find_params(jcamp_dict_bruker2d))

print("AZAMERONE NITRATE ---------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/azamerone_nitrate"):
    file = f"test_files/MNOVA_jdx/Bruker/azamerone_nitrate/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("DEMETHOXY CORNUSIDE ---------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/Demethoxy cornuside"):
    file = f"test_files/MNOVA_jdx/Bruker/Demethoxy cornuside/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("DEMETHOXY CORNUSIDE ---------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/Demethoxy cornuside"):
    file = f"test_files/MNOVA_jdx/Bruker/Demethoxy cornuside/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("GRANATICIN C ----------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/Granaticin_C"):
    file = f"test_files/MNOVA_jdx/Bruker/Granaticin_C/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("NMR RAW DATA FOR VOATRIAFRICANINE A -----------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/NMR raw data for voatriafricanine A"):
    file = f"test_files/MNOVA_jdx/Bruker/NMR raw data for voatriafricanine A/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("TETRONASIN -----------------------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Bruker/Tetronasin"):
    file = f"test_files/MNOVA_jdx/Bruker/Tetronasin/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

# print("VARIAN #############################################################################################################################################")
print("MB0593C-BAC-AA ------------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Varian/MB0593C-BAC-AA"):
    file = f"test_files/MNOVA_jdx/Varian/MB0593C-BAC-AA/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("EMERIMICIN V --------------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Varian/emerimicin V"):
    file = f"test_files/MNOVA_jdx/Varian/emerimicin V/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)

print("ECHINULIN -----------------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Varian/echinulin"):
    file = f"test_files/MNOVA_jdx/Varian/echinulin/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)


print("ASPOCHALASIN 1 -----------------------------------")
for filename in os.listdir("test_files/MNOVA_jdx/Varian/Aspochalasin I"):
    file = f"test_files/MNOVA_jdx/Varian/Aspochalasin I/{filename}"
    dict_list = Jcampdx_Handler.read([file])
    params = Jcampdx_Handler.find_params(dict_list)
    print(params)


print("JEOL ####################################################################################################################################################")
# for filename in os.listdir("test_files/Raw/JEOL"):
#     file = f"test_files/Raw/JEOL/{filename}"
#     dict_list = JEOL.read([file])
#     params = JEOL.find_params(dict_list)
#     print(params)
# with open("./test_varian_jdx.json", "w") as f:
#     comment_list = jcamp_dict_varian2d[0][0]["_datatype_LINK"][0]["_comments"]
#     for item in comment_list:
#         f.write(item+"\n")






# with open("./test_bruker_jdx.txt", "w") as f:
#     comment_list = jcamp_dict2[0][0]["_datatype_LINK"][0]["_comments"]
#     for item in comment_list:
#         f.write(item+"\n")


print("BRUKER COMBINED FILES ##############################################################################################################################")

read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Bruker/GranC.jdx"])
params = Jcampdx_Handler.find_params(read_output)
print(params)


print("VARIAN COMBINED FILES ##############################################################################################################################")

# read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Combined jdx Varian/combinedvarian.jdx"])
# params = Jcampdx_Handler.find_params(read_output)
# print(params)

read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Combined jdx Varian/combinedvarian2.jdx"])
params = Jcampdx_Handler.find_params(read_output)
print(params)
        
        
        

# read_fake_procpar(jcamp_dict)

    
# find_params_output = Jcampdx_Handler.find_params(jcamp_dict)
# print(find_params_output)

# jcamp_dict = Jcampdx_Handler.read(["test_files/Granaticin_C_Proton.jdx"])
# find_params_output = Jcampdx_Handler.find_params(jcamp_dict)
# print(find_params_output)

# print("FIRST BRUKER FILE ------\n\n")
# jcamp_dict = Jcampdx_Handler.read(["test_files/Bruker_HMBC.jdx"])
# output = Jcampdx_Handler.find_params(jcamp_dict)
# print(output)


# print("SECOND BRUKER FILE ------\n\n")
# jcamp_dict2 = Jcampdx_Handler.read(["test_files/bruker_j.jdx"])
# with open("./test_1d.json", "w") as f:
#     json.dump(jcamp_dict2[0][0], f)
# output2 = Jcampdx_Handler.find_params(jcamp_dict2)
# print(output2)


"""
NMRML Testing
"""
# nmrml_dict = nmrML.read(["/workspaces/apvalidation/test_files/nmrml/FAM013_AHTM.PROTON_04.nmrML"])

# with open("/workspaces/apvalidation/nmrml_dump.json", "w") as f:
#     json.dump(nmrml_dict,f)


# nmrml_dict2 = nmrML.read(["/workspaces/apvalidation/test_files/nmrml/HMDB00005.nmrML"])

# with open("/workspaces/apvalidation/nmrml_dump_2d.json", "w") as f:
#     json.dump(nmrml_dict2,f)



# print(nmrML.read([]))




print("CALIFORNIA DATA ###############################################################################################")
# for filename in os.listdir("test_files/North Carolina Data/1-mitragynine"):
#     file = f"test_files/North Carolina Data/1-mitragynine/{filename}"
#     read_output = Jcampdx_Handler.read([file])
#     params = Jcampdx_Handler.find_params(read_output)
#     print(params)

# for filename in os.listdir("test_files/Raw/JEOL"):
#     file = f"test_files/Raw/JEOL/{filename}"
#     dict_list = JEOL.read([file])
#     params = JEOL.find_params(dict_list)
#     print(params)