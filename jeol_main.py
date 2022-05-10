import nmrglue as ng
from apvalidation.extract import Jcampdx_Handler


file = f"/workspaces/apvalidation/apvalidation/test_files/mnova_combined_jdx.jdx"
dict_list = Jcampdx_Handler.read([file])
params = Jcampdx_Handler.find_params(dict_list)
print(params)