from apvalidation.extract import Varian, Bruker, JEOL, Jcampdx_Handler
from apvalidation.smiles_validation import validate_struct
from apvalidation.mnova_jdx_reader import separate_mnova_jdx
from rdkit import Chem
import os
import pandas as pd
import nmrglue as ng
import numpy as np
import json
import re
import nmrglue as ng


def test_jdx_single_varian():
    for compound in os.listdir("test_files/MNOVA_jdx/Varian"):
        if compound == ".DS_Store":
            continue
        
        for experiment in os.listdir(f"test_files/MNOVA_jdx/Varian/{compound}"):
            filepath = f"test_files/MNOVA_jdx/Varian/{compound}/{experiment}"
            read_output = Jcampdx_Handler.read([filepath])
            find_params_output = Jcampdx_Handler.find_params(read_output)
            print(find_params_output)

def test_jdx_combined_varian():
    read_output = Jcampdx_Handler.read(["test_files/MNOVA_jdx/Combined jdx Varian/combinedvarian.jdx"])
    print(read_output[0]['_datatype_LINK'][0].keys())
    find_params_output = Jcampdx_Handler.find_params(read_output)
    # print(find_params_output)
    pass

def print_combined_varian_struct():
    with open("test_files/MNOVA_jdx/Combined jdx Varian/combinedvarian.jdx", "r") as f:
        line_list = f.readlines()
        tabs=0
        with open("test_files/MNOVA_jdx/WriteFiles/tab_jdx.txt", "w") as fi:
            for line in line_list:
                if line.startswith("##END="):
                    tabs-=1
                if line.startswith("##TITLE="):
                    tabs+=1
                tab_line=tabs*"\t"+line
                fi.write(tab_line)



def test_recursion_sol():
    output = separate_mnova_jdx("test_files/MNOVA_jdx/Combined jdx Varian/combinedvarian.jdx")
    with open("test_files/MNOVA_jdx/WriteFiles/recurse_jdx.txt", "w") as fi:
        fi.write(json.dumps(output, sort_keys=True, indent=4))
    # print(len(output['nested_dicts'][0]['nested_dicts']))

# test_jdx_combined_varian()
# print_combined_varian_struct()
test_recursion_sol()
                

