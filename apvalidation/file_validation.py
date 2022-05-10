from apvalidation import extract as extractor
from apvalidation.simple_file_finder import MetaFinder
from apvalidation.extract_core import extract_core_file
from apvalidation.patoolutil import is_zip, repack_to_zip

# Local Test Import
# import extract as extractor
# from simple_file_finder import MetaFinder
# from extract_core import extract_core_file
# from patoolutil import is_zip, repack_to_zip

import sys
import os
import zipfile
import tempfile
import json
from pathlib import Path
import re

def find_path_and_extract(submitted_zip_file: str) -> json:
    """
    This function integrates file_finder and paramExtractor.
    :param submitted_zip_file: user submitted zip file path
    :return: Experiment parameters
    """
       
    
    if not is_zip(submitted_zip_file):
        submitted_zip_file = repack_to_zip(submitted_zip_file)
    
    meta = MetaFinder(submitted_zip_file)
    assert meta.error_message == [], meta.error_message  
    
    meta_file = meta.meta_info
    vendor_type = meta_file["vendor_name"]
    file_root = meta_file["meta_file"]
    existing_folder_names = []

    with zipfile.ZipFile(submitted_zip_file, 'r') as zipObj:
        # Extract all the contents of zip file in current directory
        res_dict = []
        # for params_path, vendor in file_root, vendor_type:
        for i, path_list in enumerate(file_root):
            unzipped_path_name = []
            for path in path_list:
                core_file_read = zipObj.read(path)
                tf = create_temporary_file(core_file_read)
                unzipped_path_name.append(tf.name)
                

            # Get param according to the vendor name
            if vendor_type[i] == "Varian":
                param_dict = extractor.Varian.read(unzipped_path_name)
                params = extractor.Varian.find_params(param_dict)
            elif vendor_type[i] == "Bruker":
                param_dict = extractor.Bruker.read(unzipped_path_name)
                params = extractor.Bruker.find_params(param_dict)
            elif vendor_type[i] == "Jcampdx":
                param_dict = extractor.Jcampdx_Handler.read(unzipped_path_name)
                params = extractor.Jcampdx_Handler.find_params(param_dict)

            # file_root_without_file_name = str(Path(path).parent)
            # res_dict[file_root_without_file_name] = params
            # res_dict[file_root_without_file_name]["vendor"] = vendor_type[i]
            
            file_root_without_file_name = str(Path(path).parent)
            if file_root_without_file_name == ".":
                file_root_without_file_name = "/"
            if type(params) == list:
                for param in params:
                    # res_dict[file_root_without_file_name] = param
                    param["original_data_path"] = file_root_without_file_name
                json_params = json.dumps(params, indent=4)
                print(json_params)
                return json_params
            
            # res_dict[file_root_without_file_name] = params
            # res_dict[file_root_without_file_name]["vendor"] = vendor_type[i]
            params["original_data_path"] = file_root_without_file_name
            params["vendor"] = vendor_type[i]
            res_dict.append(params)
            
            # # Select core files and extract under name_format directory
            # # Directory name format : <nuc_1>_<nuc_2>_<experiment_type>
            # two_d_name = res_dict[file_root_without_file_name]["nuc_2"] + "_" if res_dict[file_root_without_file_name]["nuc_2"] else ""
            # #NULL value is replaced by an empty string
            # if res_dict[file_root_without_file_name]["nuc_1"]:
            #     one_d_name = res_dict[file_root_without_file_name]["nuc_1"] + "_"
            # else:
            #     one_d_name = ""
                
            # folder_name = one_d_name + two_d_name + res_dict[file_root_without_file_name]["experiment_type"]
            # repeat_exp_num = existing_folder_names.count(folder_name)
            # existing_folder_names.append(folder_name)
            # if repeat_exp_num >= 1:
            #     folder_name = folder_name + " ({})".format(repeat_exp_num)
                        
            # parent_dir = os.getcwd()
            # indiv_exp_path = str(re.search("^(.+)/([^/]+)$", file_root[i][0])[1])
            # print("path")
            # print(indiv_exp_path)
            # extract_core_file(submitted_zip_file, indiv_exp_path, vendor_type[i], folder_name, parent_dir)

            os.unlink(tf.name) # Delete temporary file

        json_params = json.dumps(res_dict, indent=4)

        print(json_params)
        return json_params


def create_temporary_file(core_file_read):
    tf = tempfile.NamedTemporaryFile(delete=False)  # Create a temporary file that has path
    tf.write(core_file_read)  # Paste parameter file from zip file to temporary file
    assert os.path.isfile(tf.name)
    tf.close()
    return tf


if __name__ == '__main__':
    find_path_and_extract(sys.argv[1])
