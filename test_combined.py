from apvalidation.extract import Varian, Bruker, JEOL, Jcampdx_Handler
from apvalidation.smiles_validation import validate_struct
import os
import json

# read in the combined 
# with open("combined_read_output.txt", "w") as f:
#     read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Combined Test Bruker/combinedGranaticinD.jdx"])
#     print(read_output)
#     f.write(str(read_output[0]))

"""
TESTING THE COMBINED FILE WITH NMR GLUE
"""
read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Combined Test Bruker/combinedGranaticinD.jdx"])
params = Jcampdx_Handler.find_params(read_output)
print(params)
# for item in read_output[0][0]['_datatype_LINK']:

#     print(item.keys())




# with open('test_files/MNOVA_jdx/Combined Test Bruker/combinedGranaticinD.jdx', 'r') as f:
#     line_list = f.readlines()
#     section_list = []
#     section_lines = []
#     for line in line_list:
#         if line.startswith("##TITLE") or line.startswith("$$ ##TITLE"):
#             section_list.append(section_lines)
#             section_lines = []
#             section_lines.append(line)
#         else:
#             section_lines.append(line)

#     with open("combined_minus_raw.txt", "w") as write_f:
        
#         for section in section_list[1:]:
#             if section[-1].startswith("##END"):
#                 print("SECTION SKIPPED")
#                 continue
#             write_f.write("SECTION BREAK -----------------------------------------------------------------------------------------------------------------\n")
#             for line in section:
#                 write_f.write(line)
        