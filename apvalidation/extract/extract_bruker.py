import os
import nmrglue as ng
import json
from collections import OrderedDict

from .remove_unicode import remove_unicode_from_paramter_file

current_dir = os.path.dirname(__file__)
one_level_up = os.path.dirname(current_dir)
submodule_dir = os.path.join(one_level_up, 'npmrd_data_exchange')

experiment_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'experiment_standardizer.json')
with open(experiment_standardizer_path, 'r') as file:
    exp_dict = json.load(file) 

solvent_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'solvent_standardizer.json')
with open(solvent_standardizer_path, 'r') as file:
    all_solvents = json.load(file) 
    
    
class Bruker:
    """
    ALTERATIONS TO BE MADE TO THE BRUKER FUNCTION IN ORDER TO BETTER HANDLE 2D EXPERIMENTS:

    1. Change the read function so that it may take a list of filenames instead of just a single filename.
        This is to accomodate the fact that 2D data will be analyzed by looking at both acqu and acqu2

    2. Change the find_dim function to simply look at the number of files that have been read in from the read function.
        If there is one file then the dimension is 1D and if there are two files then the dimension is 2D.

    3. Change the find_nuclei function to look at both the acqu and the acqu2 files for 2D data. For now we will only consider NUC1
        to be the nuclei relevant to the experiment. This can be changed under further inspection.

    4. Find experiment type should be the same as long as find_dim works correctly as the its just looking up the key. If find_dim
        works correctly and you find this is not the case then come back to here and write your findings.
    """

    """
    A class containing the methods to help with the extraction of
    experiment parameters from NMR data from Bruker
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath_list):
        """
        Given the raw file path to an acqu file, read the file
        and return a python dictionary with all the parameters.

        :param filepath: string formatted filepath to the acqu file
        :return: dictionary containing all parameters found in the acqu file
        """
        param_dict_list = []
        for filepath in filepath_list:
            assert os.path.isfile(filepath)
            try:
                param_dict = ng.bruker.read_jcamp(filename=filepath)
            # If exception due to non-unicode character error remove it and try again
            except Exception as e:
                try:
                    if str(e).startswith("'utf-8' codec can't decode byte"):
                        remove_unicode_from_paramter_file(filepath)
                        param_dict = ng.bruker.read_jcamp(filename=filepath)
                except:
                    pass
            param_dict_list.append(param_dict)

        return param_dict_list

    @staticmethod
    def find_params(param_dict_list):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2 frequency,
        solvent, and temperature.

        :param param_dict: dictionary containing all the parameters retrieved from the file
        :return: dictionary containing only those parameters that are preferred
        """
        param_dict = param_dict_list[0]

        exp_dim = Bruker.find_dim(param_dict_list)
        exp_type = Bruker.find_exp_type(param_dict, exp_dim)        
        exp_nuc1, exp_nuc2 = Bruker.find_nuc(param_dict_list, exp_dim)
        exp_freq = Bruker.find_freq(param_dict, exp_dim)
        exp_solv = Bruker.find_solvent(param_dict)
        exp_temp = Bruker.find_temp(param_dict)

        pref_params = {
            "experiment_type": exp_type,
            "nuc_1": exp_nuc1,
            "nuc_2": exp_nuc2,
            "frequency": exp_freq,
            "solvent": exp_solv,
            "temperature": exp_temp,
        }

        return pref_params

    @staticmethod
    def find_temp(param_dict):
        try:
            te_value = param_dict.get("TE")
            if te_value is None:
                return None

            temp_number = round(float(te_value))
            if temp_number >= 250:
                return temp_number
            else:
                return temp_number + 273
        except (ValueError, TypeError):
            return None

    @staticmethod
    def find_solvent(param_dict):
        """
        Helper function. Find the solvent involved in the experiment. This is done by checking possible
        english names for the solvent against their chemical formulas.

        :param param_dict: Dictionary containing all the parameters for the experiment.
        :return: The solvent in string format.
        """
        solv_str = param_dict.get("SOLVENT")
        
        if not solv_str:
            return "FAILED_TO_DETECT"

        solv_upper = solv_str.upper()
        if solv_upper in all_solvents:
            return all_solvents[solv_upper]
        elif solv_str in all_solvents:
            return solv_str
        else:
            return "FAILED_TO_DETECT"
        
    @staticmethod
    def find_dim(param_dict_list):
        """
        Changes to be made to this function:
        1. Change the argument to a list of param_dicts. If there are more than one in the list
            then the experiment is 2D if not then it's 1D
        """
        """
        Helper function.
        Find the dimension of the experiment given a dictionary of parameters.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :return: dimension of the experiment.
        """
        exp_dim = None
        if len(param_dict_list) == 2:
            exp_dim = "2D"
        elif len(param_dict_list) == 1:
            exp_dim = "1D"

        return exp_dim

    @staticmethod
    def find_exp_type(param_dict, exp_dim):
        """
        Helper function.
        Determine the type of experiment that was conducted. Ex. COSY, HSQC etc...

        :param param_dict: dictionary containing all the parameter data
        :param exp_dim: dimension of the experiment
        :return: type of experiment in string. (1D experiments are not given a type)
        """    
        possible_exp_str_1 = param_dict["EXP"]
        possible_exp_str_2 = param_dict["PULPROG"]

        exp_type = ""
        if exp_dim == "2D":
            for entry in exp_dict:
                if (
                    entry in possible_exp_str_1.upper()
                    or entry in possible_exp_str_2.upper()
                ):
                    exp_type = exp_dict[entry]
                    return exp_type
            exp_type = "FAILED_TO_DETECT"
            return exp_type
        else:
            for entry in exp_dict:
                if entry in possible_exp_str_1 or entry in possible_exp_str_2:
                    exp_type = exp_dict[entry]
                    break
                
            if exp_type:
                exp_type = f"1D_{exp_type}".strip()
            else:
                exp_type = "1D"
                
            return exp_type

    @staticmethod
    def find_freq(param_dict, exp_dim):
        """
        Helper function.
        Find the frequency parameter given a dictionary of parameters and the
        dimension of the experiment.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :param exp_dim: the dimension of the experiment
        :return: a single float or a tuple of floats depending on the dimension
        """
        try:
            if exp_dim == "1D":
                sfo1 = param_dict.get("SFO1")
                if sfo1 is not None:
                    return [round(float(sfo1), 9)]
                else:
                    return None
            else:
                sfo1 = param_dict.get("SFO1")
                sfo2 = param_dict.get("SFO2")
                if sfo1 is not None and sfo2 is not None:
                    return (round(float(sfo1), 9), round(float(sfo2), 9))
                else:
                    return None
        except (ValueError, TypeError):
            return None

    @staticmethod
    def find_nuc(param_dict_list, exp_dim):
        """
         Helper function.
         Determine the nuclei that are used in the experiment.

        :param param_dict: a large dictionary containing all the parameters.
                             Probably returned from the read method.
         :param exp_dim: the dimension of the experiment
         :return: nucleus 1 and nucleus 2 in string format
        """        
        param_dict_1D = param_dict_list[0]
        
        param_dict_2D = {}
        if exp_dim == "2D":
            param_dict_2D = param_dict_list[1]
        
        exp_nuc1 = "FAILED_TO_DETECT"
        exp_nuc2 = None

        if exp_dim == "1D":
            exp_nuc1 = param_dict_1D.get("NUC1", "FAILED_TO_DETECT")

        elif exp_dim == "2D":
            exp_nuc1 = param_dict_1D.get("NUC1", "FAILED_TO_DETECT")
            exp_nuc2 = param_dict_2D.get("NUC1", "FAILED_TO_DETECT")

        return exp_nuc1, exp_nuc2

