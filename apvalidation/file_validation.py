from apvalidation import extract as extractor
from apvalidation.simple_file_finder import MetaFinder
from apvalidation.extract_core import extract_core_file

# import extract as extractor
# from simple_file_finder import MetaFinder
# from extract_core import extract_core_file

# Django import version
# from .paramExtract.packages import extract as extractor
# from .file_finder.simple_file_finder import MetaFinder
# from .name_format.extract_core import extract_core_file

import sys
import os
import zipfile
import tempfile
import json


def find_path_and_extract(submitted_zip_file: str) -> json:
    """
    This function integrates file_finder and paramExtractor.
    :param submitted_zip_file: user submitted zip file path
    :return: Experiment parameters
    """

    meta_file = MetaFinder(submitted_zip_file).meta_info
    vendor_type = meta_file["vendor_name"]
    file_root = meta_file["meta_file"]

    with zipfile.ZipFile(submitted_zip_file, 'r') as zipObj:
        # Extract all the contents of zip file in current directory
        res_dict = {}
        # for params_path, vendor in file_root, vendor_type:
        for i in range(len(file_root)):
            core_file_read = zipObj.read(file_root[i])
            tf = create_temporary_file(core_file_read)

            # Below Python 3.10 doesn't support match case
            # match vendor_type[i]:
            #     case "Varian":
            #         param_dict = paramExtract.extract.Varian.read(tf.name)
            #         params = paramExtract.extract.Varian.find_params(param_dict)
            #     case "Bruker":
            #         param_dict = paramExtract.extract.Bruker.read(tf.name)
            #         params = paramExtract.extract.Bruker.find_params(param_dict)
            #     case "JEOL":
            #         param_dict = paramExtract.extract.Jcampdx.read(tf.name)
            #         params = paramExtract.extract.Jcampdx.find_params(param_dict)

            # Get param according to the vendor name
            if vendor_type[i] == "Varian":
                param_dict = extractor.Varian.read(tf.name)
                params = extractor.Varian.find_params(param_dict)
            elif vendor_type[i] == "Bruker":
                param_dict = extractor.Bruker.read(tf.name)
                params = extractor.Bruker.find_params(param_dict)
            elif vendor_type[i] == "JEOL":
                param_dict = extractor.Jcampdx.read(tf.name)
                params = extractor.Jcampdx.find_params(param_dict)

            res_dict[file_root[i]] = params
            res_dict[file_root[i]]["vendor"] = vendor_type[i]
            # Select core files and extract under name_format directory
            # Directory name format : <nuc_1>_<nuc_2>_<experiment_type>
            two_d_name = res_dict[file_root[i]]["nuc_2"] + "_" if res_dict[file_root[i]]["nuc_2"] else ""
            #NULL value is replaced by an empty string
            if res_dict[file_root[i]]["nuc_1"]:
                one_d_name = res_dict[file_root[i]]["nuc_1"] + "_"
            else:
                one_d_name = ""
            folder_name = one_d_name + two_d_name + res_dict[file_root[i]]["experiment_type"]
            parent_dir = os.getcwd()
            # extract_core_file(submitted_zip_file, file_root[i], vendor_type[i], folder_name, parent_dir)

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
