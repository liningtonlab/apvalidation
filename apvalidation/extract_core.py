import re
from zipfile import ZipFile
import sys
import shutil
import os

varian = ["fid", "procpar", "log", "text", "log 2", "procpar 2", "text 2"]
bruker = ["fid", "ser", "acqu", "acqu2", "acqus", "acqu2s"]
jcamp = ["jdx"]
dir_path = os.path.dirname(os.path.realpath(__file__))


def extract_core_file(input_zip, target_exp, vendor, folder_name):
    with ZipFile(input_zip, 'r') as zipObject:
        all_paths = zipObject.namelist()

        # Below Python 3.10 doesn't support match case
        # match vendor:
        #     case "Varian":
        #         extract_files = search_keyword(all_paths, target_exp, varian)
        #     case "Bruker":
        #         extract_files = search_keyword(all_paths, target_exp, bruker)
        #     case "JEOL":
        #         extract_files = search_keyword(all_paths, target_exp, jcamp)

        # Get paths of necessary files

        if vendor == "Varian":
            extract_files = search_keyword(all_paths, target_exp, varian)
        elif vendor == "Bruker":
            extract_files = search_keyword(all_paths, target_exp, bruker)
        elif vendor == "JEOL":
            extract_files = search_keyword(all_paths, target_exp, jcamp)

        extract_to_folder(extract_files, zipObject, folder_name)


def extract_to_folder(extract_files, zipObject, folder_name):

    os.makedirs(os.path.join(dir_path, folder_name), exist_ok=True)
    target_path = os.path.join(dir_path, folder_name)

    for unzip in extract_files:
        member = zipObject.open(unzip)

        with open(os.path.join(target_path, os.path.basename(unzip)), 'wb') as out_file:
            # copy it directly to the output directory,
            # without creating the intermediate directory
            shutil.copyfileobj(member, out_file)


def search_keyword(all_paths, target_exp, vendor):
    unzip_path = []
    target_exp = re.search("^(.+)/([^/]+)$", target_exp)[1]
    for path in all_paths:
        for item in vendor:
            if path.endswith(item) and path.startswith(target_exp):
                unzip_path.append(path)
    return unzip_path


if __name__ == '__main__':
    extract_core_file(sys.argv[1], sys.argv[2])
