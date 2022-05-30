from posixpath import split
from apvalidation.extract import Varian, Bruker, JEOL, Jcampdx_Handler
from apvalidation.smiles_validation import validate_struct
from apvalidation.mnova_jdx_reader import separate_mnova_jdx, add_indents
from rdkit import Chem
import os
import pandas as pd
import nmrglue as ng
import numpy as np
import json
import re
import nmrglue as ng


def test_split():
    """
        take a combined jdx file from any of the three vendors and split it up and save it in a folder
    """
    save_location = separate_mnova_jdx("test_files/Sheets Testing/JEOl mnovatojdx combined/JEOLSheetTesting1.jdx", "test_files/testing_save_singles/test_folder_5")
    return save_location

def test_read(save_location):
    for filename in os.listdir(save_location):
        print(filename)
        read_output = Jcampdx_Handler.read([f"{save_location}/{filename}"])
        find_params_output = Jcampdx_Handler.find_params(read_output)
        print(find_params_output)

split_output = test_split()
test_read(split_output)

                

