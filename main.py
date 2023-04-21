import os
from apvalidation.extract import Jcampdx_Handler
from apvalidation.mnova_jdx_reader import separate_mnova_jdx
from apvalidation.peak_validator import Validate
from apvalidation.file_validation import find_path_and_extract

# input_file = "test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx"
# save_location = "test_files/testing_save_singles/test_folder_6"
# saved_path = separate_mnova_jdx(input_file, save_location)
# print(saved_path)

# for filename in os.listdir(saved_path):
#     read_output = Jcampdx_Handler.read([os.path.join(saved_path, filename)])
#     find_param_output = Jcampdx_Handler.find_params(read_output)
#     print(f"find_paramoutput = {find_param_output}")


# result = Validate.validate("13.0, (13.5-13.9), 12.4 - 14.6, 10.0", "(11.5-15.2), 64.7", "CC", "D2O")
# print("result")
# print(result)

# from apvalidation.extract import Varian

# input_file = "./procpar"
# Varian.remove_personal_info(input_file)

# find_path_and_extract("./apvalidation/test/jdx_split/TEST_JEOL.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/jdx_split/TEST ROGER JDX proton_only.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/jdx_split/TEST ROGER JDX FULL.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/jdx_split/TEST ROGER JDX NO COSY.zip", is_second_time = False)

# find_path_and_extract("./apvalidation/test/jdx_split/SSW03M2 D3-3 2.zip", is_second_time = False)

# find_path_and_extract("./apvalidation/test/ycld34.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/Borrelidin_B_NMR_RAW.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/Borrelidin_B.mnova", is_second_time = False)
find_path_and_extract("./apvalidation/test/compound_zip.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/Granaticin_C.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/Borrelidin_B_NMR_RAW.zip", is_second_time = False)
# find_path_and_extract("./apvalidation/test/compound_KZFMHNJGNPWLAV-DEAAMUPLSA-N.zip", is_second_time = False)


# find_path_and_extract("./apvalidation/test/jdx_split/TEST ROGER JDX PROTON_NOSEY.zip", is_second_time = False)

# find_path_and_extract("./apvalidation/test/original_data_SIIRGYXKTQZYMV-DGPRMSGASA-N.zip", is_second_time = False)
