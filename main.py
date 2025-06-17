import os
import nmrglue as ng
import apvalidation
from apvalidation.extract.extract_jcampdx import Jcampdx
from apvalidation.extract_core import extract_jdx
from apvalidation.mnova_jdx_reader import separate_mnova_jdx
from apvalidation.peak_validator import Validate as Peak_Validate
from apvalidation.file_validation import find_path_and_extract
from apvalidation.smiles_to_inchikey import to_inchikey

# input_file = "test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx"
# save_location = "test_files/testing_save_singles/test_folder_6"
# saved_path = separate_mnova_jdx(input_file, save_location)
# print(saved_path)

# for filename in os.listdir(saved_path):
#     read_output = Jcampdx.read([os.path.join(saved_path, filename)])
#     find_param_output = Jcampdx.find_params(read_output)
#     print(f"find_paramoutput = {find_param_output}")


# result = Peak_Validate.validate(
#     H_text_block="13.0,dhg (13.5-13.9), 12.4sdf - 14.6, 10.0",
#     C_text_block="20, 20, 20",
#     smiles="CCCC=CCCC=CCC\CCC/CCC\CCCC",
#     solvent="D2O",
#     h_frequency=400,
#     h_temperature=300,
#     c_frequency=300,
#     c_temperature=300,
#     reference="DMSO",
# )



# result = Peak_Validate.validate(
#     H_text_block="",
#     C_text_block="10.1, 20, 20",
#     smiles="CCCC=CCCC=CCC\CCC/CCC\CCCC",
#     solvent="D2O",
#     h_frequency=None,
#     h_temperature=None,
#     c_frequency=300,
#     c_temperature=300,
#     reference="DMSO",
# )
# print("peak result")
# print(result)

# from apvalidation.extract_varian import Varian

# input_file = "./procpar"
# Varian.remove_personal_info(input_file)

# metadata = find_path_and_extract("./apvalidation/test/test_bruker_jdx.zip", is_second_time = False)
# print(metadata)


metadata = find_path_and_extract("./apvalidation/test/test_roger_jdx.zip", is_second_time = False)
print(metadata)
# print("---------------------------------------------------------")
# metadata = find_path_and_extract("./apvalidation/test/error_jdx.zip", is_second_time = False)
# print(metadata)


# regular_jdx = ng.jcampdx.read(filename="/workspaces/apvalidation/apvalidation/test/JEOL/ I1_85_02_ PULSE ACQUISITION_exp_1.jdx")
# print("regular_jdx is")
# print(regular_jdx)

# print("\n--------------------\n")

# regular_jdx = ng.jcampdx.read(filename="/workspaces/apvalidation/apvalidation/test/JEOL_DX/1H NMR to check_ zg30_exp_1.dx")
# print("DX is")
# print(regular_jdx)


# print("\n\n\n-------JEOL.zip-------")
# metadata = find_path_and_extract("./apvalidation/test/JEOL.zip", is_second_time = False)
# print(metadata)
# extract_jdx(
#     "./apvalidation/test/",
#     "JEOL.zip",
#     "jeol_output",
#     "./apvalidation/test/output/",
# )



# print("\n\n\n-------test_bruker_jdx.zip-------")
# metadata = find_path_and_extract("./apvalidation/test/test_bruker_jdx.zip", is_second_time = False)
# print(metadata)


# print("\n\n\n-------RGL1617G1B.zip-------")
# metadata = find_path_and_extract("./apvalidation/test/RGL1617G1B.zip", is_second_time = False)
# print(metadata)


# print("\n\n\n-------JEOL_DX.zip-------")
# metadata = find_path_and_extract("./apvalidation/test/JEOL_DX.zip", is_second_time = False)
# print(metadata)
# extract_jdx(
#     "./apvalidation/test/",
#     "JEOL_DX.zip",
#     "jeol_output",
#     "./apvalidation/test/output/",
# )

# metadata = find_path_and_extract("./apvalidation/test/JEOL.zip")
# print(metadata)

# print("----")

# metadata = find_path_and_extract("./apvalidation/test/salarin_C_failed_exp_type.zip")
# print(metadata)

# metadata = find_path_and_extract("./apvalidation/test/endolide_E_NMR_RAW_HMBC_ONLY.zip")
# print(metadata)

# test_dir_path = "./apvalidation/test"
# for filename in os.listdir(test_dir_path):
#     file_path = os.path.join(test_dir_path, filename)

#     # Check if the file ends with ".zip"
#     if filename.endswith(".zip") and os.path.isfile(file_path):
#         # Execute the code on the ZIP file
#         print("--------------------------------")
#         print(f"processing file {file_path}")
#         try:
#             metadata = find_path_and_extract(file_path, is_second_time = False)
#             print("VVVVVVVVVVVV sucessfully processed VVVVVVVVVVV")
#             print(metadata)
#         except Exception as e:
#             print("FFFFFFFFFFFFF failed to process FFFFFFFFFFFFF")
#             print(e)

# metadata = find_path_and_extract("./apvalidation/test/original_data_LELBFTMXCIIKKX-QVRQZEMUSA-N.zip")
# print(metadata)



# metadata = find_path_and_extract("./apvalidation/test/test_dept.zip")
# print(metadata)


print("test to_inchikey(smiles)")
to_inchikey("CC=CCC=CC/CCC=CC\CC=CC/CC=CC\CC=CCCCC/CC=CCCC\CC=CC/CCC")