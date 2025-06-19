# from apvalidation.extract import Varian, Bruker, JEOL, Jcampdx
from apvalidation.extract.extract_bruker import Bruker as bruker_extractor
from apvalidation.extract.extract_varian import Varian as varian_extractor
from apvalidation.extract.extract_jeol import JEOL as jeol_extractor
from apvalidation.extract.extract_jcampdx import Jcampdx as jcampdx_extractor
from apvalidation.simple_file_finder import MetaFinder
from apvalidation.extract_core import extract_core_file
from apvalidation.patoolutil import is_zip, repack_to_zip, is_compressed_but_not_zip
from apvalidation.mnova_jdx_reader import separate_mnova_jdx

# Local Test Import
# from apvalidation.extract import Varian, Bruker, Jcampdx
# from simple_file_finder import MetaFinder
# from extract_core import extract_core_file, extract_jdx
# from patoolutil import is_zip, repack_to_zip
# from mnova_jdx_reader import separate_mnova_jdx

import sys
import os
import zipfile
import tempfile
import json
from pathlib import Path
import re

zip_file_extention = [
    ".zip",
    ".rar",
    ".7z",
    ".tar",
    ".ace",
    ".adf",
    ".alz",
    ".ape",
    ".a",
    ".arc",
    ".arj",
    ".bz2",
    ".cab",
    ".Z",
    ".cpio",
    ".deb",
    ".dms",
    ".flac",
    ".gz",
    ".iso",
    ".lrz",
    ".lha",
    ".lzh",
    ".lz",
    ".lzma",
    ".lzo",
    ".rpm",
    ".rz",
    ".shn",
    ".xz",
    ".jar",
    ".zoo",
    ".zpaq",
]


def find_path_and_extract(
    submitted_zip_file: str,
    is_second_time: bool = False
) -> json:
    """
    This function integrates file_finder and paramExtractor.
    :param submitted_zip_file: user submitted zip file path
    :return: Experiment parameters
    """

    if not is_zip(submitted_zip_file):
        if is_compressed_but_not_zip(
            submitted_zip_file, zip_file_extention
        ) == True:
            old_zip_file = submitted_zip_file
            submitted_zip_file = repack_to_zip(submitted_zip_file)
            os.unlink(old_zip_file)
        else:
            raise ValueError(json.dumps({
                "Compression Error": {
                    submitted_zip_file: "File is not in a compressed (.zip, .rar, .7z, etc.) format."
                    + " Please compress your NMR data before uploading."
                }
            }))

    print("submitted_zip_file is", submitted_zip_file)

    meta = MetaFinder(
        submitted_zip_file,
        zip_file_extention
    )
    
    assert meta.error_message == {}, json.dumps(meta.error_message)

    meta_file = meta.meta_info
    print("meta_file is", meta_file)
    vendor_type = meta_file["vendor_name"]
    filetype = meta_file["filetype"]
    file_root = meta_file["file_root"]
    existing_folder_names = []
    jcamp = False

    with zipfile.ZipFile(submitted_zip_file, "r") as zipObj:
        # Extract all the contents of zip file in current directory
        res_dict = []
        for i, vendor in enumerate(vendor_type):
            if vendor == "Jcampdx" and len(file_root[i]) > 1:
                vendor_type.pop(i)
                filetype.pop(i)
                tmp = file_root.pop(i)
                for j in range(len(tmp)):
                    vendor_type.append("Jcampdx")
                    filetype.append("Jcampdx")
                    file_root.append([tmp[j]])
        for i, path_list in enumerate(file_root):
            unzipped_path_name = []
            
            for path in path_list:
                core_file_read = zipObj.read(path)
                tf = create_temporary_file(core_file_read)
                unzipped_path_name.append(tf.name)
            
            # Get param according to the vendor name
            if vendor_type[i] == "Varian":
                param_dict = varian_extractor.read(unzipped_path_name)
                params = varian_extractor.find_params(param_dict)
                add_path_vendor(path, params, vendor_type[i], filetype[i], res_dict)
            elif vendor_type[i] == "Bruker":
                param_dict = bruker_extractor.read(unzipped_path_name)
                params = bruker_extractor.find_params(param_dict)
                add_path_vendor(path, params, vendor_type[i], filetype[i], res_dict)

            # file_root_without_file_name = str(Path(path).parent)

            if vendor_type[i] == "Jcampdx":
                jcamp = True
                jcamp_file_extension = path_list[0].split(".")[-1]
                
                # spilit_file_dir = f"{str(Path(submitted_zip_file).parent)}/jdx_spilt"
                loc = os.path.splitext(submitted_zip_file)[0]
                
                if not is_second_time:
                    loc = separate_mnova_jdx(unzipped_path_name[0], loc, jcamp_file_extension)

            os.unlink(tf.name)  # Delete temporary file
        
        if jcamp:
            res_dict = extract_jcamp(loc)
        
        json_params = json.dumps(res_dict, indent=4)

        # for i, vendor in enumerate(vendor_type):
        #     # Select core files and extract under name_format directory
        #     # Directory name format : <nuc_1>_<nuc_2>_<experiment_type>
        #     two_d_name = res_dict[i]["nuc_2"] + "_" if res_dict[i]["nuc_2"] else ""
        #     #NULL value is replaced by an empty string
        #     if res_dict[i]["nuc_1"]:
        #         one_d_name = res_dict[i]["nuc_1"] + "_"
        #     else:
        #         one_d_name = ""

        #     folder_name = one_d_name + two_d_name + res_dict[i]["experiment_type"]
        #     repeat_exp_num = existing_folder_names.count(folder_name)
        #     existing_folder_names.append(folder_name)
        #     if repeat_exp_num >= 1:
        #         folder_name = folder_name + " ({})".format(repeat_exp_num)

        #     parent_dir = os.getcwd()
        #     # indiv_exp_path = str(re.search("^(.+)/([^/]+)$", file_root[i][0])[1])
        #     indiv_exp_path = f"{str(Path(file_root[i][0]).parent)}/"
        #     param_file = res_dict[i]['original_data_path']
        #     # extract_core_file(submitted_zip_file, indiv_exp_path, vendor_type[i], folder_name, parent_dir, param_file)
        #     extract_jdx(loc,param_file,folder_name,parent_dir)

        return json_params




def extract_jcamp(loc):
    print("loc is", loc)
    res_dict = []
    for path in os.listdir(loc):
        if Path(path).suffix == ".jdx" or Path(path).suffix == ".dx":
            full_path = os.path.join(loc, path)
            print("full_path is", full_path)
            param_dict, json_nmr_data_dict = jcampdx_extractor.read(full_path)
            print("json_nmr_data_dict.keys() is")
            print(json_nmr_data_dict.keys())
            manuf = jcampdx_extractor.find_manuf(param_dict=param_dict, json_nmr_data_dict=json_nmr_data_dict)
            print("manuf is", manuf)
            found_params = jcampdx_extractor.find_params(
                param_dict,
                json_nmr_data_dict=json_nmr_data_dict,
                manuf=manuf
            )
            params = found_params[0]
            add_path_vendor(path, params, manuf, "Jcampdx", res_dict)
    
    print("res_dict is")
    print(res_dict)
    return res_dict


def create_temporary_file(core_file_read):
    tf = tempfile.NamedTemporaryFile(
        delete=False
    )  # Create a temporary file that has path
    tf.write(core_file_read)  # Paste parameter file from zip file to temporary file
    assert os.path.isfile(tf.name)
    tf.close()
    return tf


def add_path_vendor(path, params, vendor_type, filetype, res_dict):    
    file_root_without_file_name = str(path)
    if file_root_without_file_name == ".":
        file_root_without_file_name = "/"
    params["original_data_path"] = file_root_without_file_name.strip()
    params["vendor"] = vendor_type
    params["filetype"] = filetype
    res_dict.append(params)


if __name__ == "__main__":
    find_path_and_extract(sys.argv[1])
    # find_path_and_extract("/Users/jonghyeokkim/Downloads/NMR/12-speciofoline.zip", True)
    # extract_jcamp(sys.argv[1])
