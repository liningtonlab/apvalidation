import os
from apvalidation.extract import Jcampdx_Handler
from apvalidation.mnova_jdx_reader import separate_mnova_jdx
from apvalidation.peak_validator import Validate

# input_file = "test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx"
# save_location = "test_files/testing_save_singles/test_folder_6"
# saved_path = separate_mnova_jdx(input_file, save_location)
# print(saved_path)

# for filename in os.listdir(saved_path):
#     read_output = Jcampdx_Handler.read([os.path.join(saved_path, filename)])
#     find_param_output = Jcampdx_Handler.find_params(read_output)
#     print(f"find_paramoutput = {find_param_output}")


result = Validate.validate("13.0, (13.5-13.9), 12.4 - 14.6, 10.0", "(11.5-15.2), 64.7", "CC", "D2O")
print("result")
print(result)
