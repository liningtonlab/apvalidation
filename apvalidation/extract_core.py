import re
from zipfile import ZipFile
import sys
import shutil
import os

varian = ["fid", "procpar", "log", "text", "log 2", "procpar 2", "text 2"]
bruker = ["fid", "ser", "acqu", "acqu2", "acqus", "acqu2s"]
jcamp = ["jdx"]


def extract_core_file(input_zip : str, indiv_exp_path, vendor, folder_name, parent_dir):
    with ZipFile(input_zip, 'r') as zipObject:
        all_paths = zipObject.namelist()

        # Get paths of necessary files
        if vendor == "Varian":
            extract_files = search_keyword(all_paths, indiv_exp_path, varian)
        elif vendor == "Bruker":
            extract_files = search_keyword(all_paths, indiv_exp_path, bruker)
        elif vendor == "JEOL":
            extract_files = search_keyword(all_paths, indiv_exp_path, jcamp)

        extract_to_folder(extract_files, zipObject, folder_name, parent_dir)


def extract_to_folder(extract_files, zipObject, folder_name, parent_dir):

    target_path = create_dir(parent_dir, folder_name, -1)

    for unzip in extract_files:
        member = zipObject.open(unzip)

        with open(os.path.join(target_path, os.path.basename(unzip)), 'wb') as out_file:
            # copy it directly to the output directory,
            # without creating the intermediate directory
            shutil.copyfileobj(member, out_file)


def search_keyword(all_paths, indiv_exp_path, vendor):
    unzip_path = []
    # print(indiv_exp_path)
    # indiv_exp_path = re.search("^(.+)/([^/]+)$", indiv_exp_path)[1]
    for path in all_paths:
        for item in vendor:
            if path.endswith(item) and path.startswith(indiv_exp_path):
                unzip_path.append(path)
    return unzip_path

def create_dir(parent_dir, folder_name, index):
    index += 1
    target_path = os.path.join(parent_dir, folder_name)
    try:
       if(index != 0):
           folder_name = re.search(".*[^ \(\d\)]", folder_name)[0]
           folder_name = f"{folder_name} ({index})"
       os.makedirs(os.path.join(parent_dir, folder_name), exist_ok=False)
       target_path = os.path.join(parent_dir, folder_name)
    except:
        create_dir(parent_dir, folder_name, index)
    return target_path

if __name__ == '__main__':
    extract_core_file(sys.argv[1], sys.argv[2])
