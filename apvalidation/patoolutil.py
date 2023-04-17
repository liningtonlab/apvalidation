import patoolib
import sys
from io import StringIO
import os
from pathlib import Path
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def return_name_list(zip_file_path:str):
    tag_list = []
    save_stdout = sys.stdout
    result = StringIO()
    sys.stdout = result
    patoolib.list_archive(zip_file_path)
    sys.stdout = save_stdout
    tag_list.append(result.getvalue())
    file_list = tag_list[0].split('\n')
    file_list = list(filter(None, file_list))
    return file_list

def repack_to_zip(zip_file_path:str):
    parent_path = os.path.dirname(zip_file_path)
    path_without_extention = os.path.splitext(zip_file_path)[0]
    file_name = os.path.basename(path_without_extention)
    unzip_file_name = f"{parent_path}/{file_name}.zip"
    try:
        patoolib.repack_archive(zip_file_path, unzip_file_name)
    except:
        random_name = id_generator()
        unzip_file_name = f"{parent_path}/{file_name}-{random_name}.zip"
        patoolib.repack_archive(zip_file_path, unzip_file_name)
        
    return unzip_file_name

def is_zip(zip_file_path:str):
    return Path(zip_file_path).suffix.lower() == ".zip"

def is_compressed_but_not_zip(zip_file_path:str, compression_extensions: list):
    file_extension = Path(zip_file_path).suffix.lower()
    for extension in compression_extensions:
        if file_extension == extension:
            if file_extension != ".zip":
                print(f"file extension match as {file_extension}")
                return True
    print(f"file extension match not found")
    return False
    
if __name__ == "__main__":
    # res = repack_to_zip(sys.argv[1])
    # file_list = return_name_list(sys.argv[1])
    pass
