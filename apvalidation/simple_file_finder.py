import os
import sys
from zipfile import ZipFile
import re
from typing import Dict, List

class MetaFinder:

    """
    When a user submits a zip file of MNR data, script unzips it in memory and find path of parameter files.
    """

    meta_name_by_vendor = {
        ".jdf": "Jcampdx", ".jdx": "Jcampdx", "acqu": "Bruker", "procpar": "Varian", "acqu2": "Bruker"
    }
    

    def __init__(self, input_zip: str, zip_file_extention: list):
        self.zip_file_extention = zip_file_extention
        self.error_message = {}
        self.all_file_path = self.get_all_file_path(input_zip)
        self.meta_info = self.find_meta(self.all_file_path)
        self.validator(self.all_file_path)
    
    def append_error_message(self, error_type: str = None, file: str = None, message: str = None):
        if error_type not in self.error_message:
            self.error_message[error_type] = {}
        if file not in self.error_message[error_type]:
            self.error_message[error_type][file] = []
        self.error_message[error_type][file].append(message)
        
    # Get all file path list from the submitted zip file
    def get_all_file_path(self, input_zip: str) -> list:
        input_zip = os.path.normpath(input_zip)
        submitted_zip = ZipFile(input_zip)
        path_in_zip = submitted_zip.namelist()
        
        # MAC zip saves some files that interrupt nmr glue, this script exclude the interrupting files
        all_path_list = []
        for path in path_in_zip:
            if re.search("__MACOSX",path) is None:
                all_path_list.append(path)
                
        return all_path_list
    
    # Given the path list, determine vendor and param file path for each experiment(directory)
    def find_meta(self, all_path_list: str) -> Dict[str, list]:
        
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

        return meta_info
    
    def validator(self, all_path_list):
        # If meta data file is not found, raise an assertion
        if not self.meta_info["meta_file"]:
            # Raise an error if known error are found
            self._vendor_not_found_error(all_path_list)
            if not self.error_message:
                # No known error are found
                self.append_error_message('File Format Error', 'Uploaded File', 'Only Varian, Jcampdx, and Bruker files are accepted at this time')

        # Based on found meta data, go through file validation
        # for vendor in meta_info["vendor_name"]:
        for i in range(len(self.meta_info["vendor_name"])):
            parent_dir = re.search("^(.+)/([^/]+)$", self.meta_info["meta_file"][i][0])
            target_exp = parent_dir[1] if parent_dir is not None else ""
          
            if self.meta_info["vendor_name"][i] == "Varian":
                self._varian_validation(all_path_list, target_exp)
            elif self.meta_info["vendor_name"][i] == "Bruker":
                self._bruker_validation(all_path_list, target_exp)
            elif self.meta_info["vendor_name"][i] == "Jcampdx":
                self._jcampdx_validation(all_path_list, target_exp)


    @staticmethod
    def param_file_finder(path_list: str, keyword: str, core_path_dict: dict) -> List[str]:
        vendor_list = []
        for path in path_list:
            if path.endswith(keyword):
                parent_dir = os.path.dirname(path)
                try:
                    core_path_dict[parent_dir]
                    core_path_dict[parent_dir].append(path)
                except:
                    core_path_dict[parent_dir] = [path]
                    vendor_list.append(MetaFinder.meta_name_by_vendor[keyword])
                    
        return vendor_list

    @staticmethod
    def key_file_finder(path_list: str, keyword: str, start_with: str) -> List[str]:
        key_path_list = []
        for path in path_list:
            if path.endswith(keyword) and path.startswith(start_with):
                key_path_list.append(path)
        return key_path_list

    def _varian_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        # assert fid_path, f"{individual_folder_path} : Fid file is missing"
        if not fid_path:
            self.append_error_message('Failed to read NMR Data', individual_folder_path, 'Fid file is missing or in an invalid directory location')
        

    def _bruker_validation(self, all_path_list: str, individual_folder_path: str):
        fid_path = self.key_file_finder(all_path_list, "fid", individual_folder_path)
        ser_path = self.key_file_finder(all_path_list, "ser", individual_folder_path)
        # assert fid_path+ser_path, f"{individual_folder_path} : Fid/Ser file is missing"
        if not fid_path and not ser_path :
            self.append_error_message('Failed to read NMR Data', individual_folder_path, 'Fid file is missing or in an invalid directory location')

    def _jcampdx_validation(self, all_path_list: str, individual_folder_path: str):
        jdx_path = self.key_file_finder(all_path_list, "jdx", individual_folder_path)
        # assert jdx_path, f"{individual_folder_path} : .jdf is not supported. Please convert to .jdx file"
        if not jdx_path :
            self.append_error_message('File Format Error', individual_folder_path, '.jdf is not supported. Please convert to .jdx files using the export function in JEOL Delta or MestreNova')

    def _vendor_not_found_error(self, all_path_list: str):
        self._invalid_file_detector(all_path_list, '.mnova', '.mnova is not currently supported. Please submit original raw NMR files.')
        self._invalid_file_detector(all_path_list, '.nmrML', '.nmrML is not currently supported. Please submit original raw NMR files or .jdx files instead.')
        for extention in self.zip_file_extention:
            self._invalid_file_detector(all_path_list, extention, f'Please make sure that the submission does not include nested {extention} file. You can put all original NMR files in the same zip folder')
        
    def _invalid_file_detector(self, all_path_list: str, keyword : str, error_message : str):
        # self.error_message.extend([f"{path} : {error_message}" for path in all_path_list if path.endswith(keyword)])
        for path in all_path_list:
            if path.endswith(keyword):
                self.append_error_message('Vendor Not Supported Error', path, error_message)

    
if __name__ == '__main__':
    res = MetaFinder(sys.argv[1]).meta_info
    print(res)