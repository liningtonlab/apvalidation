import os
import re
import nmrglue as ng
import pandas as pd
import fileinput
import sys
import shutil
import json

current_dir = os.path.dirname(__file__)
one_level_up = os.path.dirname(current_dir)
submodule_dir = os.path.join(one_level_up, 'npmrd_data_exchange')

experiment_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'experiment_standardizer.json')
with open(experiment_standardizer_path, 'r') as file:
    exp_dict = json.load(file) 

solvent_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'solvent_standardizer.json')
with open(solvent_standardizer_path, 'r') as file:
    all_solvents = json.load(file) 

class Varian:
    """
    A class that contains the methods to extract parameters from Varian NMR data
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath_list):
        """
        Read the file to retrieve a dictionary of experiment parameters.

        :param filepath: a string formatted filepath to the procpar file
        :return: dictionary containing all parameters found in the procpar file OR
                'Empty File' if the file is empty
        """
        param_dict_list = []
        for filepath in filepath_list:
            assert os.path.isfile(filepath)
            param_dict = ng.varian.read_procpar(filename=filepath)
            param_dict_list.append(param_dict)

        return param_dict_list

    @staticmethod
    def find_params(param_dict_list):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2 frequency,
        solvent, and temperature.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :return: dictionary containing only the preferred parameters
        """

        param_dict = param_dict_list[0]
        exp_dim = Varian.find_dim(param_dict)
        exp_type = Varian.find_exp_type(param_dict, exp_dim)
        exp_freq = Varian.find_freq(param_dict, exp_dim)
        exp_nuc1, exp_nuc2 = Varian.find_nuc(param_dict, exp_dim)
        exp_solv = Varian.find_solvent(param_dict)
        exp_temp = Varian.find_temp(param_dict)

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
    def remove_personal_info(filepath):
        """
        Runs repalce_propcar_property on a list of fields known to contain
        sensetive information in a propcar file. Replaces them with empty
        strings in the file directly.

        :param filepath: string formatted filepath to the propcar file
        :return: None
        """
        Varian.replace_procpar_property(filepath, "go_id", 1, '""')
        Varian.replace_procpar_property(filepath, "emailaddr", 1, '""')

    @staticmethod
    def replace_procpar_property(
        filepath, tag, replace_index_after_tag, replacement_line
    ):
        """
        Replaces specifed properties in a propcar filed with a provided string.

        :param filepath: string formatted filepath to the propcar file
        :param tag: string of the property to replace
        :param replace_index_after_tag: the index of the specific element in the property
        to replace
        :param replacement_line: the string to replace the specified property with element with

        :return: None
        """
        # Make copy of original file to preserve it if an error occurs
        file_split = os.path.splitext(filepath)
        copypath = f"{file_split[0]}_copy{file_split[1]}"
        shutil.copy(filepath, copypath)

        try:
            # Run replacement on copy
            replace = False
            for line in fileinput.input(copypath, inplace=1):
                leading_str = ""
                if line.startswith("$$"):
                    leading_str = f"$$ {str(replace_index_after_tag)}"
                elif re.match("^[0-9 ]+\W", line):
                    leading_str = str(replace_index_after_tag)

                if (
                    replace is True
                    and line.startswith(leading_str)
                    and leading_str != ""
                ):
                    replace = False
                    line = f"{leading_str} {replacement_line}\n"
                if tag in line:
                    replace = True

                sys.stdout.write(line)
        except:
            os.remove(copypath)
            return False

        os.remove(filepath)
        os.rename(copypath, filepath)

        return True

    @staticmethod
    def find_temp(param_dict):
        try:
            temp_val = param_dict["temp"]["values"][0]
            temp_number = float(temp_val)
            if temp_number >= 250:
                return temp_number
            else:
                return temp_number + 273
        except:
            return None

    @staticmethod
    def find_solvent(param_dict):
        """
        Helper function.
        Searches for the solvent key and then attempts to match it to a solvent in the dictionary
        of solvents. First attempts to match to all possible english spellings of the solvent, and then
        attempts to match to all possible chemical formulas.

        :param param_dict: a large dictionary containing all the parameters.
        :return: a string containing the chemical formula of the solvent used.
        """
        try:
            solv_str = param_dict["solvent"]["values"][0]
            if solv_str is not None:
                solv_str = solv_str.upper()
        except KeyError:
            solv_str = None

        if solv_str.upper() in all_solvents.keys():
            exp_solv = all_solvents[solv_str.upper()]
        elif solv_str in all_solvents.values():
            exp_solv = solv_str
        else:
            exp_solv = "FAILED_TO_DETECT"
        return exp_solv

    @staticmethod
    def find_dim(param_dict):
        """
        Helper function.
        Find the dimension of the experiment given a dictionary of parameters.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :return: dimension of the experiment
        """
        
        if "plotoption" in param_dict.keys():
            for option in param_dict['plotoption']['values']:
                if option.startswith("plot2D"):
                    exp_dim = "2D"
                    return exp_dim

        # if "plt2Darg" in param_dict.keys():
        #     exp_dim = "2D"
        #     return exp_dim

        try:
            exp_dim = param_dict["apptype"]["values"][0][-2:].upper()
        except KeyError:
            exp_dim = None

        if exp_dim in ["1D", "2D"]:
            return exp_dim

        elif exp_dim not in ["1D", "2D"]:
            try:
                exp_dim = str(param_dict["procdim"]["values"][0]) + "D"
            except KeyError:
                exp_dim = None

        if exp_dim not in ["1D", "2D", None]:
            exp_dim = None
        return exp_dim

    @staticmethod
    def find_exp_type(param_dict, exp_dim):
        """
        Helper function.
        Determine the type of experiment that was conducted. Ex. COSY, HSQC etc...
        The first bit of the function constructs a list of possible locations in the parameter
        file which could house the name of the experiment being conducted.
        These locations are then searched one of the keywords that denote an experiment.

        :param param_dict: dictionary containing all the parameter data
        :param exp_dim: dimension of the experiment
        :return: type of experiment in string. (1D experiments are not given a type)
        """
        """
        CHANGES WERE MADE TO THIS FUNCTION IN ORDER TO CATCH SOME MISSING EXPERIMENT TYPES:
        1. Lines 156 to 159 were added to add an alternative source from which to take the experiment type string.
            This was added due to the lack of information in the other key on some files. 
            An example of where this helps is in Lobosamide C for the HMBC and the HSQC. 

        2. IMPORTANT UPDATE: The above does not work. Try running the main and find that all the experiments are now labelled 
            as HSQCTCOSY even for those that are not. Confusion. Try again lol.
        """

        try:
            exp_type_loc1 = param_dict["explist"]["values"][0]
        except KeyError:
            exp_type_loc1 = ""

        try:
            exp_type_loc2 = param_dict["apptype"]["values"][0]
        except KeyError:
            exp_type_loc2 = ""
        try:
            long_string = param_dict["ap"]["values"][0]
            start_indicator = long_string.find("pwx:3;1:") + len("pwx:3;1:")
            end_indicator = long_string.find(":j1xh:")
            exp_type_loc3 = long_string[start_indicator:end_indicator]
        except KeyError:
            exp_type_loc3 = ""

        try:
            exp_type_loc4 = param_dict["pslabel"]["values"][0]
        except KeyError:
            exp_type_loc4 = ""

        exp_loc_list = [exp_type_loc1, exp_type_loc2, exp_type_loc3, exp_type_loc4]
        exp_loc_list = [x.upper() for x in exp_loc_list]

        exp_type = ""
        if exp_dim == "2D":
            for exp_loc in exp_loc_list:
                for type_str in exp_dict:
                    if type_str in exp_loc:
                        exp_type = exp_dict[type_str]
                        return exp_type
        elif exp_dim == "1D":
            for exp_loc in exp_loc_list:
                for type_str in exp_dict:
                    if type_str in exp_loc:
                        exp_type = exp_dict[type_str]
                        break
            
            if exp_type:
                exp_type = f"1D_{exp_type}".strip()
            else:
                exp_type = "1D"
                
            return exp_type
        else:
            exp_type = "FAILED_TO_DETECT"
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
            if exp_dim == "2D":
                freq1 = round(float(param_dict["reffrq"]["values"][0]), 9)
                freq2 = round(float(param_dict["reffrq1"]["values"][0]), 9)
                return (freq1, freq2)
            else:
                freq = round(float(param_dict["reffrq"]["values"][0]), 9)
                return [freq]
        except:
            return (None, None)

    @staticmethod
    def find_nuc(param_dict, exp_dim):
        """
        Find the nuclei included in the Varian experiment.

        :param param_dict: dict containing param data
        :param exp_dim: dimension of the experiment
        :return: dimension of the experiment
        """
        explist_keyword_dict = {"PROTON": "1H", "CARBON": "13C"}
        rewrite_nucs_dict = {"H1": "1H", "C13": "13C", "N15": "15N"}

        if exp_dim == "1D":
            try:
                atom_name = param_dict["explist"]["values"][0]
                exp_nuc1 = explist_keyword_dict[atom_name]
                exp_nuc2 = None
            except KeyError:
                atom_name = param_dict["tn"]["values"][0]
                exp_nuc1 = rewrite_nucs_dict[atom_name]
                exp_nuc2 = None

        elif exp_dim == "2D":
            exp_nuc1, exp_nuc2 = Varian.help_find_2d_nuc(param_dict)
        else:
            exp_nuc1 = None
            exp_nuc2 = None

        return exp_nuc1, exp_nuc2

    @staticmethod
    def help_find_2d_nuc(param_dict):
        """
        Helper function. A function to help the find_nuc function to calculate the correct nuclei in Varian's 2D case.
        In this case, finding the nuclei is especially complicated and therefore a second
        helper function is needed.

        :param param_dict: a dictionary containing the parameters of the experiment.
        :return: the two nuclei involved in the experiment.
        """
        nuc_dict = {
            "H1": [
                (0.99875, 1.00125),
                (3.971172962, 3.98111332),
                (9.864197531, 9.888888889),
                (1.061652936, 1.064310391),
                (2.467572576, 2.473749228),
            ],
            "C13": [
                (0.25025, 0.25275),
                (0.9950298211, 1.004970179),
                (2.471604938, 2.496296296),
                (0.2660111613, 0.2686686155),
                (0.6182828907, 0.6244595429),
            ],
            "N15": [
                (0.1, 0.1025),
                (0.3976143141, 0.407554672),
                (0.987654321, 1.012345679),
                (0.1062981664, 0.1089556205),
                (0.2470660902, 0.2532427424),
            ],
            "F": [
                (0.9395, 0.942),
                (3.735586481, 3.745526839),
                (9.279012346, 9.303703704),
                (0.9986712729, 1.001328727),
                (2.321185917, 2.327362569),
            ],
            "P": [
                (0.4035, 0.406),
                (1.604373757, 1.614314115),
                (3.985185185, 4.009876543),
                (0.4289131012, 0.4315705554),
                (0.9969116739, 1.003088326),
            ],
        }

        nuc_df = pd.DataFrame(nuc_dict)
        nuc_df.index = ["H1", "C13", "N15", "F", "P"]

        nuc_1 = param_dict["tn"]["values"][0]
        nuc_2 = ""
        freq_ratio = float(param_dict["reffrq"]["values"][0]) / float(
            param_dict["reffrq1"]["values"][0]
        )

        for frq_range in nuc_df[nuc_1]:
            if frq_range[0] < freq_ratio < frq_range[1]:
                nuc_2 = nuc_df.index[nuc_df[nuc_1] == frq_range][0]

        rewrite_nucs_dict = {"H1": "1H", "C13": "13C", "N15": "15N"}

        if nuc_1 in rewrite_nucs_dict.keys():
            nuc_1 = rewrite_nucs_dict[nuc_1]
        if nuc_2 in rewrite_nucs_dict.keys():
            nuc_2 = rewrite_nucs_dict[nuc_2]

        return nuc_1, nuc_2