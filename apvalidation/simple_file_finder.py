import os
import sys
from zipfile import ZipFile
import re

class MetaFinder:

    """
    When a user submits a zip file of MNR data, script unzips it in memory and find path of parameter files.
    """

    meta_name_by_vendor = {
        ".jdf": "JEOL", ".jdx": "JEOL", "acqu": "Bruker", "procpar": "Varian"
    }

    def __init__(self, input_zip: str):
        self.meta_info = self.find_meta(input_zip)

    def find_meta(self, input_zip: str) -> dict:
        input_zip = os.path.normpath(input_zip)
        submitted_zip = ZipFile(input_zip)
        all_path_list = submitted_zip.namelist()
        meta_file_name_list = list(MetaFinder.meta_name_by_vendor.keys())

        # Search for a meta data file names
        vendor_list = []
        param_path_list = []
        for name in meta_file_name_list:
            lst = self.param_file_finder(all_path_list, name)
            vendor_list += lst[0]
            param_path_list += lst[1]
        meta_info = {"vendor_name": vendor_list, "meta_file": param_path_list}

        # If meta data file is not found, raise an assertion
        assert meta_info["meta_file"], "Only Varian, JEOL, Bruker files are accepted"

        # Based on found meta data, go through file validation
        # for vendor in meta_info["vendor_name"]:
        for i in range(len(meta_info["vendor_name"])):
            parent_dir = re.search("^(.+)/([^/]+)$", meta_info["meta_file"][i])
            target_exp = parent_dir[1] if parent_dir is not None else ""

            # match meta_info["vendor_name"][i]:
            #     case "Varian":
            #         self.__varian_validation(all_path_list, target_exp)
            #     case "Bruker":
            #         self.__bruker_validation(all_path_list, target_exp)
            #     case "JEOL":
            #         self.__jeol_validation(all_path_list, target_exp)

            if meta_info["vendor_name"][i] == "Varian":
                self.__varian_validation(all_path_list, target_exp)
            elif meta_info["vendor_name"][i] == "Bruker":
                self.__bruker_validation(all_path_list, target_exp)
            elif meta_info["vendor_name"][i] == "JEOL":
                self.__jeol_validation(all_path_list, target_exp)

        return meta_info

    @staticmethod
    def param_file_finder(path_list: str, keyword: str) -> list:
        core_path_list = []
        vendor_list = []
        for path in path_list:
            if path.endswith(keyword):
                vendor_list.append(MetaFinder.meta_name_by_vendor[keyword])
                core_path_list.append(path)
        return [vendor_list, core_path_list]

    @staticmethod
    def key_file_finder(path_list: str, keyword: str, start_with: str) -> list:
        key_path_list = []
        for path in path_list:
            if path.endswith(keyword) and path.startswith(start_with):
                key_path_list.append(path)
        return key_path_list

    def __varian_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        assert fid_path, "fid file is missing"

    def __bruker_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        ser_path = self.key_file_finder(all_path_list, "ser", individual_folder_path)
        assert fid_path+ser_path, "fid file is missing"

    def __jeol_validation(self, all_path_list: str, individual_folder_path: str):
        jdx_path = self.key_file_finder(all_path_list, "jdx", individual_folder_path)
        assert jdx_path, ".jdf is not supported. Convert to .jdx file"


if __name__ == '__main__':
    res = MetaFinder(sys.argv[1]).meta_info
    print(res)