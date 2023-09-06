import os
from apvalidation.extract_classes import Jcampdx
from apvalidation.mnova_jdx_reader import separate_mnova_jdx
from apvalidation.peak_validator import Validate
from apvalidation.file_validation import find_path_and_extract

# input_file = "test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx"
# save_location = "test_files/testing_save_singles/test_folder_6"
# saved_path = separate_mnova_jdx(input_file, save_location)
# print(saved_path)

# for filename in os.listdir(saved_path):
#     read_output = Jcampdx.read([os.path.join(saved_path, filename)])
#     find_param_output = Jcampdx.find_params(read_output)
#     print(f"find_paramoutput = {find_param_output}")


# result = Validate.validate("13.0, (13.5-13.9), 12.4 - 14.6, 10.0", "(11.5-15.2), 64.7", "CC", "D2O")
# print("result")
# print(result)

# from apvalidation.extract import Varian

# input_file = "./procpar"
# Varian.remove_personal_info(input_file)

metadata = find_path_and_extract("./apvalidation/test/test_bruker_jdx.zip", is_second_time = False)
print(metadata)

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
            


