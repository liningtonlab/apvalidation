import os
import sys
from zipfile import ZipFile
import re

class MetaFinder:

    """
    When a user submits a zip file of MNR data, script unzips it in memory and find path of parameter files.
    """

    meta_name_by_vendor = {
        ".jdf": "Jcampdx", ".jdx": "Jcampdx", "acqu": "Bruker", "procpar": "Varian", "acqu2": "Bruker"
    }
    zip_file_extention = [".7z",".ace", ".adf",".alz",".ape",".a",".arc", ".arj", ".bz2",".cab", ".Z",
                          ".cpio",".deb",".dms",".flac",".gz",".iso",".lrz", ".lha", ".lzh", ".lz", ".lzma", 
                          ".lzo", ".rpm", ".rar", ".rz", ".shn", ".tar", ".xz", ".zip", ".jar", ".zoo", ".zpaq"]
    

    def __init__(self, input_zip: str):
        self.error_message = []
        self.meta_info = self.find_meta(input_zip)

    def find_meta(self, input_zip: str) -> dict:
        input_zip = os.path.normpath(input_zip)
        submitted_zip = ZipFile(input_zip)
        path_in_zip = submitted_zip.namelist()
        
        # MAC zip saves some files that interrupt nmr glue, this script exclude the interrupting files
        all_path_list = []
        for path in path_in_zip:
            if re.search("__MACOSX",path) is None:
                all_path_list.append(path)
        
        meta_file_name_list = list(MetaFinder.meta_name_by_vendor.keys())

        # Search for a meta data file names
        vendor_list = []
        param_path_list= []
        core_path_dict = {}
        for name in meta_file_name_list:
            lst = self.param_file_finder(all_path_list, name, core_path_dict)
            if lst:
                vendor_list += lst
        param_path_list = list(core_path_dict.values())         
        meta_info = {"vendor_name": vendor_list, "meta_file": param_path_list}

        # If meta data file is not found, raise an assertion
        if not meta_info["meta_file"]:
            # Raise an error if known error are found
            self._vendor_not_found_error(path_in_zip)
            if not self.error_message:
                # No known error are found
                self.error_message.append("Only Varian, Jcampdx, Bruker files are accepted")

        # Based on found meta data, go through file validation
        # for vendor in meta_info["vendor_name"]:
        for i in range(len(meta_info["vendor_name"])):
            parent_dir = re.search("^(.+)/([^/]+)$", meta_info["meta_file"][i][0])
            target_exp = parent_dir[1] if parent_dir is not None else ""

          
            if meta_info["vendor_name"][i] == "Varian":
                self._varian_validation(all_path_list, target_exp)
            elif meta_info["vendor_name"][i] == "Bruker":
                self._bruker_validation(all_path_list, target_exp)
            elif meta_info["vendor_name"][i] == "Jcampdx":
                self._jcampdx_validation(all_path_list, target_exp)

        return meta_info

    @staticmethod
    def param_file_finder(path_list: str, keyword: str, core_path_dict: dict) -> list:
        vendor_list = []
        for path in path_list:
            if path.endswith(keyword):
                parent_dir = re.search("^(.+)/([^/]+)$", path)[1]
                try:
                    core_path_dict[parent_dir]
                    core_path_dict[parent_dir].append(path)
                except:
                    core_path_dict[parent_dir] = [path]
                    vendor_list.append(MetaFinder.meta_name_by_vendor[keyword])
                    
        return vendor_list

    @staticmethod
    def key_file_finder(path_list: str, keyword: str, start_with: str) -> list:
        key_path_list = []
        for path in path_list:
            if path.endswith(keyword) and path.startswith(start_with):
                key_path_list.append(path)
        return key_path_list

    def _varian_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        # assert fid_path, f"{individual_folder_path} : Fid file is missing"
        if not fid_path : self.error_message.append(f"{individual_folder_path} : Fid file is missing")
        

    def _bruker_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        ser_path = self.key_file_finder(all_path_list, "ser", individual_folder_path)
        # assert fid_path+ser_path, f"{individual_folder_path} : Fid/Ser file is missing"
        if not fid_path and not ser_path : self.error_message.append(f"{individual_folder_path} : Fid file is missing")


    def _jcampdx_validation(self, all_path_list: str, individual_folder_path: str):
        jdx_path = self.key_file_finder(all_path_list, "jdx", individual_folder_path)
        # assert jdx_path, f"{individual_folder_path} : .jdf is not supported. Please convert to .jdx file"
        if not jdx_path : self.error_message.append(f"{individual_folder_path} : .jdf is not supported. Please convert to .jdx files using the export function in JEOL Delta")

    def _vendor_not_found_error(self, all_path_list: str):
        self._invalid_file_detector(all_path_list, '.mnova', '.mnova is not currently supported. Please submit original raw NMR files.')
        self._invalid_file_detector(all_path_list, '.nmrML', '.nmrML is not currently supported. Please submit original raw NMR files or .jdx files instead.')
        for extention in self.zip_file_extention:
            self._invalid_file_detector(all_path_list, extention, f'Please make sure that the submission does not include nested {extention} file. You can put all original NMR files in the same zip folder')
        
    def _invalid_file_detector(self, all_path_list: str, keyword : str, error_message : str):
        self.error_message.extend([f"{path} : {error_message}" for path in all_path_list if path.endswith(keyword)])        

    

if __name__ == '__main__':
    res = MetaFinder(sys.argv[1]).meta_info
    print(res)