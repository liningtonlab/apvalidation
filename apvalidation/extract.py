import os
import re
import nmrglue as ng
import pandas as pd
import numpy as np
import fileinput
import sys
import shutil
from collections import OrderedDict


all_solvents = {
    **dict.fromkeys(
        [
            "DMSO",
            "DIMETHYL SULFOXIDE",
            "DIMETHYL-SULFOXIDE",
            "DMSO-D6",
            "dmso-d6",
            "d6-dmso",
            "D6-DMSO",
        ],
        "C2D6OS",
    ),
    **dict.fromkeys(
        ["TETRAHYDROFURAN", "THF", "D6-THF", "THF-D6", "C2D6OS_SPE"], "C2D6OS"
    ),
    **dict.fromkeys(
        ["CHLOROFORM", "CHLOROFORM-D", "DEUTEROCHLOROFORM", "CDCL3_SPE", "CDCL3"],
        "CDCl3",
    ),
    **dict.fromkeys(
        [
            "DICHLOROMETHANE",
            "DEUTERATED-DICHLOROMETHANE",
            "METHYLENE CHLORIDE",
            "METHYLENE-CHLORIDE",
            "METHYLENE-CHLORI",
            "CD2CL2_SPE",
        ],
        "CD2Cl2",
    ),
    **dict.fromkeys(["ACETONE", "C3D6O_SPE", "ACETONE-D6"], "C3D6O"),
    **dict.fromkeys(["DEUTERATED-METHANOL", "MEOD", "CD3OD_SPE"], "CD3OD"),
    **dict.fromkeys(["TOLUENE", "C7D8_SPE"], "C7D8"),
    **dict.fromkeys(
        [
            "HEAVY-WATER",
            "OXIDANE",
            "DEUTERIUM-OXIDE",
            "H2O+D2O",
            "DEUTERIUM OXIDE",
            "D2O_SPE",
            "H2O+D2O_SALT_3MM",
        ],
        "D2O",
    ),
    **dict.fromkeys(
        ["TRIFLUOROACETIC-ACID", "TRIFLUOROACETIC ACID", "C2DF3O2_SPE"], "C2DF3O2"
    ),
    **dict.fromkeys(["PYRIDINE", "PYR", "PYRIDINE-D5", "C5D5N_SPE"], "C5D5N"),
    **dict.fromkeys(["ACETONITRILE", "C2D3N_SPE"], "C2D3N"),
    **dict.fromkeys(["ACETONITRILE-D3", "ACETONITRILE D3", "CD3CN_SPE"], "CD3CN"),
    **dict.fromkeys(["BENZENE", "C6D6_SPE", "BENZENE-D6"], "C6D6"),
    **dict.fromkeys(["METHANOL", "MEOH", "METHANOL-D3"], "CH3OH"),
    **dict.fromkeys(["DIMETHYLFORMAMIDE", "DMF"], "C3H7NO"),
    **dict.fromkeys(
        [
            "<no solvents available",
            "<<no solvents available>",
            "<no solvents available>",
            "<<no solvents available>>",
        ],
        None,
    ),
}


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
        temp_number = int(param_dict["temp"]["values"][0])
        if temp_number >= 250:
            return temp_number
        else:
            return temp_number + 273

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

        if solv_str in all_solvents.keys():
            exp_solv = all_solvents[solv_str].upper()
        elif solv_str in all_solvents.values():
            exp_solv = solv_str.upper()
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

        if "plt2Darg" in param_dict.keys():
            exp_dim = "2D"
            return exp_dim

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

        # exp_list = ['HSQCTOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'DOSY', 'ROESY', 'NOESY']
        ordered_exp_dict = OrderedDict(
            [
                ("HSQCTOCSY", "HSQCTOCSY"),
                ("COSY", "COSY"),
                ("HSQC", "HSQC"),
                ("HMQC", "HMQC"),
                ("HMBC", "HMBC"),
                ("HSQMBC", "HSQMBC"),
                ("TOCSY", "TOCSY"),
                ("HETLOC", "TOCSY"),
                ("MLEVPHSW", "TOCSY"),
                ("DOSY", "DOSY"),
                ("ROESY", "ROESY"),
                ("NOESY", "NOESY"),
                ("DEPT", "DEPT"),
                ("H2BC", "H2BC"),
                ("HOMO2DJ", "HOMO2DJ"),
            ]
        )

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
                for type_str in ordered_exp_dict:
                    if type_str in exp_loc:
                        exp_type = ordered_exp_dict[type_str]
                        return exp_type
        elif exp_dim == "1D":
            for exp_loc in exp_loc_list:
                for type_str in ordered_exp_dict:
                    if type_str in exp_loc:
                        exp_type = ordered_exp_dict[type_str]
                        break
            exp_type = f"1D {exp_type}".strip()
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

        if exp_dim == "2D":
            freq1 = round(float(param_dict["reffrq"]["values"][0]), 9)
            freq2 = round(float(param_dict["reffrq1"]["values"][0]), 9)
            freq_val = (freq1, freq2)
        else:
            freq_val = [round(float(param_dict["reffrq"]["values"][0]), 9)]
        return freq_val

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
            param_dict = ng.bruker.read_jcamp(filename=filepath)
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
        temp_number = round(float(param_dict["TE"]))
        if temp_number >= 250:
            return temp_number
        else:
            return temp_number + 273

    @staticmethod
    def find_solvent(param_dict):
        """
        Helper function. Find the solvent involved in the experiment. This is done by checking possible
        english names for the solvent against their chemical formulas.

        :param param_dict: Dictionary containing all the parameters for the experiment.
        :return: The solvent in string format.
        """

        solv_str = param_dict["SOLVENT"]

        if solv_str.upper() in all_solvents.keys():
            exp_solv = all_solvents[solv_str.upper()]
        elif solv_str.upper() in all_solvents.values():
            exp_solv = solv_str
        else:
            exp_solv = "FAILED_TO_DETECT"
        return exp_solv

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

        # exp_2d_list = ['HSQCTOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'ROESY', 'NOESY', 'DOSY']

        ordered_exp_dict = OrderedDict(
            [
                ("HSQCTOCSY", "HSQCTOCSY"),
                ("COSY", "COSY"),
                ("HSQC", "HSQC"),
                ("HMQC", "HMQC"),
                ("HMBC", "HMBC"),
                ("TOCSY", "TOCSY"),
                ("HETLOC", "TOCSY"),
                ("MLEVPHSW", "TOCSY"),
                ("ROESY", "ROESY"),
                ("NOESY", "NOESY"),
                ("DOSY", "DOSY"),
                ("DEPT", "DEPT"),
                ("H2BC", "H2BC")
            ]
        )

        exp_type = ""
        if exp_dim == "2D":
            for entry in ordered_exp_dict:
                if (
                    entry in possible_exp_str_1.upper()
                    or entry in possible_exp_str_2.upper()
                ):
                    exp_type = ordered_exp_dict[entry]
                    return exp_type
            exp_type = "FAILED_TO_DETECT"
            return exp_type
        else:
            for entry in ordered_exp_dict:
                if entry in possible_exp_str_1 or entry in possible_exp_str_2:
                    exp_type = ordered_exp_dict[entry]
                    break
            exp_type = f"1D {exp_type}".strip()
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
        if exp_dim == "1D":
            freq_val = round(float(param_dict["SFO1"]), 9)
            return [freq_val]
        else:
            freq1 = round(float(param_dict["SFO1"]), 9)
            freq2 = round(float(param_dict["SFO2"]), 9)
            freq_val = (freq1, freq2)
            return freq_val

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
        if exp_dim == "2D":
            param_dict_2D = param_dict_list[1]

        exp_nuc2 = None
        if exp_dim == "1D":
            exp_nuc1 = param_dict_1D["NUC1"]
        elif exp_dim == "2D":
            exp_nuc1 = param_dict_1D["NUC1"]
            exp_nuc2 = param_dict_2D["NUC1"]
        else:
            exp_nuc1 = None

        return exp_nuc1, exp_nuc2


class JEOL:
    """
    A class containing the methods to help with the extraction of
    experiment parameters from NMR data from JEOL format.
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath_list):
        """
        Given the raw file path to a parameter file, read the file
        and return a python dictionary with all the parameters.

        :param filepath: string formatted filepath to the acqu file
        :return: dictionary containing all parameters found in the acqu file
        """

        param_dict_list = []
        for filepath in filepath_list:
            assert os.path.isfile(filepath)
            param_dict = ng.jcampdx.read(filename=filepath)
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
        try:
            param_dict = param_dict_list[0][0]
        except:
            param_dict = param_dict_list[0]
        try:
            param_dict = param_dict["_datatype_NMRSPECTRUM"][0]
        except KeyError:
            try:
                param_dict = param_dict["_datatype_LINK"][0]
            except KeyError:
                param_dict = param_dict

        exp_dim = JEOL.find_dim(param_dict)
        exp_freq = JEOL.find_freq(param_dict, exp_dim)
        exp_nuc_1, exp_nuc_2 = JEOL.find_nuc(param_dict, exp_dim)
        exp_type = JEOL.find_exp_type(param_dict, exp_dim)
        exp_solv = JEOL.find_solvent(param_dict)
        exp_temp = JEOL.find_temp(param_dict)

        pref_params = {
            "experiment_type": exp_type,
            "nuc_1": exp_nuc_1,
            "nuc_2": exp_nuc_2,
            "frequency": exp_freq,
            "solvent": exp_solv.upper(),
            "temperature": exp_temp,
        }

        return pref_params

    @staticmethod
    def find_temp(param_dict):
        try:
            temp_number = int(param_dict["$TEMPSET"][0])
            if temp_number >= 250:
                return temp_number
            else:
                return temp_number + 273
        except TypeError:
            return None

    @staticmethod
    def find_solvent(param_dict):
        """
        Helper function. Find the solvent involved in the experiment. This is done by checking possible
        english names for the solvent against their chemical formulas.

        :param param_dict: Dictionary containing all the parameters for the experiment.
        :return: The solvent in string format.
        """
        solv_str = param_dict["$SOLVENT"][0]

        if solv_str in all_solvents.keys():
            exp_solv = all_solvents[solv_str]
        elif solv_str.upper() in all_solvents.values():
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
        :return: dimension of the experiment.
        """

        try:
            exp_dim = param_dict["NUMDIM"][0] + "D"
        except KeyError:
            exp_dim = None
        if exp_dim == None:
            try:
                exp_dim = param_dict["$DIMENSIONS"][0] + "D"
            except KeyError:
                exp_dim = None

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
        exp_str = param_dict[".PULSESEQUENCE"][0].upper()

        # exp_2d_list = ['HSQC-TOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'ROESY', 'NOESY', 'DOSY']
        ordered_exp_dict = OrderedDict(
            [
                ("HSQC-TOCSY", "HSQCTOCSY"),
                ("COSY", "COSY"),
                ("HSQC", "HSQC"),
                ("HMQC", "HMQC"),
                ("HMBC", "HMBC"),
                ("TOCSY", "TOCSY"),
                ("HETLOC", "TOCSY"),
                ("MLEVPHSW", "TOCSY"),
                ("ROESY", "ROESY"),
                ("NOESY", "NOESY"),
                ("DOSY", "DOSY"),
                ("DEPT", "DEPT"),
                ("H2BC", "H2BC"),
            ]
        )

        exp_type = "FAILED_TO_DETECT"

        if exp_dim == "2D":
            for entry in ordered_exp_dict:
                if entry in exp_str.upper():
                    exp_type = ordered_exp_dict[entry]
                    return exp_type
            return exp_type
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

        if exp_dim == "2D":
            freq1 = round(float(param_dict["$XFREQ"][0]), 9)
            freq2 = round(float(param_dict["$YFREQ"][0]), 9)
            freq = (freq1, freq2)
            return freq
        else:
            freq = round(float(param_dict["$XFREQ"][0]), 9)
            return [freq]

    @staticmethod
    def find_nuc(param_dict, exp_dim):
        """
         Helper function.
         Determine the nuclei that are used in the experiment.

        :param param_dict: a large dictionary containing all the parameters.
                             Probably returned from the read method.
         :param exp_dim: the dimension of the experiment
         :return: nucleus 1 and nucleus 2 in string format
        """
        if exp_dim == "2D":
            try:
                nuc_1 = param_dict[".NUCLEUS"][0].split(", ")[1]
                nuc_2 = param_dict[".NUCLEUS"][0].split(", ")[0]
            except:
                pass
            try:
                nuc_1 = param_dict[".NUCLEUS"][0].split(",")[1]
                nuc_2 = param_dict[".NUCLEUS"][0].split(",")[0]
            except:
                pass
        else:
            nuc_1 = param_dict[".OBSERVENUCLEUS"][0][1:]
            nuc_2 = None

        return nuc_1, nuc_2


class Jcampdx_Handler:
    """
    A class containing the methods to help with the extraction of
    experiment parameters from NMR data from jdx format. This class handles jcamps
    by following these steps:
    1. Classify the file as Bruker, Varian.
    2. Format the jdx file to match the nmrglue read outputs.
    3. Feed those formatted file into the respective Class (Varian, Bruker).
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath_list):
        """
        Given the raw file path to a parameter file, read the file
        and return a python dictionary with all the parameters.

        :param filepath: string formatted filepath to the acqu file
        :return: dictionary containing all parameters found in the acqu file
        """

        param_dict_list = []
        for filepath in filepath_list:
            assert os.path.isfile(filepath)
            param_dict = ng.jcampdx.read(filename=filepath)
            param_dict_list.append(param_dict)

        # Band-aid to make para_dict consistent since some .jdx files produce a
        # nested tuple/dict and some don't.
        try:
            param_dict = param_dict_list[0][0]["_datatype_NMRSPECTRUM"]
        except:
            param_dict = param_dict_list[0]

        return param_dict

    @staticmethod
    def find_manuf(param_dict):
        """
        find the manufacturer that produced the data file that is being inspected.

        :param param_dict_list: a list of dictionaries read in from the read function
        :return: the name of the manufacturer
        """

        try:
            manuf_name = param_dict[0]["_datatype_LINK"][0]["$ORIGINALFORMAT"][0]
        except (KeyError, TypeError):
            manuf_name = "Not found"
            try:
                manuf_name = param_dict[0]["$ORIGINALFORMAT"][0]
            except (KeyError, TypeError):
                manuf_name = "Not found"
                try:
                    manuf_name = param_dict[0]["ORIGIN"][0]
                except (KeyError, TypeError):
                    manuf_name = "Not found"

        if manuf_name == "Varian":
            return manuf_name
        elif "Bruker" in manuf_name:
            manuf_name = "Bruker"
            return manuf_name
        elif manuf_name in ["JCAMP-DX NMR", "JEOL Delta", "DELTA2_NMR", "DELTA_NMR"]:
            manuf_name = "JEOL"
            return manuf_name
        else:
            manuf_name = "Not found"
            return manuf_name

    @staticmethod
    def find_params(param_dict):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2, frequency,
        solvent, and temperature.

        :param param_dict: dictionary containing all the parameters retrieved from the file
        :return: dictionary containing only those parameters that are preferred
        """

        output_list = []
        manuf = Jcampdx_Handler.find_manuf(param_dict)

        if manuf == "Varian":
            varian_structured_dict_list = Jcampdx_Handler.format_varian(param_dict)
            output_list.append(Varian.find_params(varian_structured_dict_list))

        elif manuf == "Bruker":
            bruker_structured_dict_list = Jcampdx_Handler.format_bruker(param_dict)
            output_list.append(Bruker.find_params(bruker_structured_dict_list))

        elif manuf == "JEOL":
            jeol_structured_dict_list = Jcampdx_Handler.format_jeol_combined(param_dict)
            output_list.append(JEOL.find_params(jeol_structured_dict_list))

        return output_list

    @staticmethod
    def format_varian(jdx_read_output):
        """
        Take the output produced from Jcamp read and format it to be passed into the Varian Class methods.
        :param jdx_read_output: a nested list object, the output from the Jcamp read method.
        :return: list of dictionaries, these are formatted for the Varian Class methods.
        """
        errored = False
        try:
            line_list = jdx_read_output[0]["_datatype_LINK"][0]["_comments"]
        except KeyError or TypeError:
            errored = True
        if errored == True:
            try:
                line_list = jdx_read_output[0]["_comments"]
            except KeyError or TypeError:
                pass

        key_holder = None
        param_dict = {}
        value_stack = []

        for line in line_list:
            line = line.replace("\n", "")
            if line[0].isalpha():
                key = line.split(" ")[0]

                if key_holder is None:
                    key_holder = key
                else:
                    value_dict = {"values": value_stack}
                    param_dict[key_holder] = value_dict
                    value_stack = []
                    key_holder = key
            else:
                line = line[2:]
                line = line.replace('"', "")
                value_stack.append(line)

        return [param_dict]

    @staticmethod
    def format_bruker(read_jdx_output):
        """
        Take the output produced from Jcamp read and format it to be passed into the Bruker Class methods.
        This function needs to separate parts of the this output to find 2 acqu files if there are 2 of them.
        :param jdx_read_output: a nested list object, the output from the Jcamp read method.
        :return: list of dictionaries, these are formatted for the Bruker Class methods.
        """

        param_dict = read_jdx_output[0]

        line_list = []
        try:
            line_list = param_dict["_comments"]
        except KeyError:
            line_list = "Not found"
        if line_list == "Not found":
            try:
                line_list = param_dict["_datatype_LINK"][0]["_comments"]
            except KeyError:
                line_list = "Not found"

        file_seps = []

        for index, item in enumerate(line_list):
            if "##TITLE= " in item:
                file_seps.append([index, 0])
            if "##END=" in item:
                matching_list = file_seps[-1]
                matching_list[1] = index + 1

        file_list = []
        for curr_file in file_seps:
            file_list.append(line_list[curr_file[0] : curr_file[1]])

        if "Parameter file" in file_list[1][0]:
            dim = "2D"
        else:
            dim = "1D"

        if dim == "1D":
            line_list = file_list[0]
            param_dict = {}
            for line in line_list:
                if line.startswith("##"):
                    split_line = line.split("=")
                    key = split_line[0]
                    value = split_line[1]

                    key = key.replace("#", "")
                    key = key.replace("$", "")
                    value = re.sub("[< | >]*", "", value)

                    param_dict[key] = value

            return [param_dict]

        elif dim == "2D":
            dim_1_line_list = file_list[1]
            dim_2_line_list = file_list[0]
            param_dict_dim1 = {}
            param_dict_dim2 = {}

            for line in dim_1_line_list:
                if line.startswith("##"):
                    split_line = line.split("=")
                    key = split_line[0]
                    value = split_line[1]

                    key = key.replace("#", "")
                    key = key.replace("$", "")
                    value = re.sub("[< | >]*", "", value)

                    param_dict_dim1[key] = value
            for line in dim_2_line_list:
                if line.startswith("##"):
                    split_line = line.split("=")
                    key = split_line[0]
                    value = split_line[1]

                    key = key.replace("#", "")
                    key = key.replace("$", "")
                    value = re.sub("[< | >]*", "", value)

                    param_dict_dim2[key] = value
            return [param_dict_dim1, param_dict_dim2]

    @staticmethod
    def format_jeol_combined(jdx_read_output):
        """
        Take the output produced from Jcamp read and format it to be passed into the Varian Class methods.

        :param jdx_read_output: a nested list object, the output from the Jcamp read method.
        :return: list of dictionaries, these are formatted for the Varian Class methods.
        """
        try:
            param_dict = jdx_read_output[0]["_datatype_LINK"]
        except KeyError:
            try:
                param_dict = jdx_read_output
            except KeyError:
                param_dict = "None"

        # re-name and format frequency keys so they match the delta version of JEOL data.
        try:
            freq_list = param_dict[0][".OBSERVEFREQUENCY"][1]
            freq_list = freq_list.split(",")
            param_dict[0]["$XFREQ"] = [freq_list[0]]
            param_dict[0]["$YFREQ"] = [freq_list[1]]
        except:
            freq_list = param_dict[0][".OBSERVEFREQUENCY"][0]
            param_dict[0]["$XFREQ"] = [freq_list]

        # add SOLVENT keys
        param_dict[0]["$SOLVENT"] = param_dict[0][".SOLVENTNAME"]

        # set the temperature to None for now
        param_dict[0]["$TEMPSET"] = [None]

        return [jdx_read_output]
