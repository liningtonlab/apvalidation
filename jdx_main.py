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

def test_bruker_structure(filepath):

    with open(f"{filepath}", "r") as f:
        line_list = f.readlines()
        line_list_tabs = add_indents(line_list)
    with open(f"test_files/bruker_j.jdx", "w+") as f:
        for line in line_list_tabs:
            f.write(line)

def test_split():
    output = separate_mnova_jdx("test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx", "test_files/testing_save_singles/test_folder_2")
    return output

def test_read(split_output):
    read_output = Jcampdx_Handler.read(["test_files/testing_save_singles/test_folder_2/I1_85_17.jdx"])
    find_params_output = Jcampdx_Handler.find_params(read_output)
    return find_params_output

split_output = test_split()
print(test_read(split_output))
# test_bruker_structure(r"test_files/MNOVA_jdx/Combined JEOL jdx/combinedJEOLpart2.jdx")
                

