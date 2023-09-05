import os
import nmrglue as ng
import json

with open('apvalidation/metadata_standardizers/experiment_standardizer.json', 'r') as file:
    exp_dict = json.load(file) 
    
with open('apvalidation/metadata_standardizers/solvent_standardizer.json', 'r') as file:
    all_solvents = json.load(file)
    

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

        exp_type = "FAILED_TO_DETECT"

        if exp_dim == "2D":
            for entry in exp_dict:
                if entry in exp_str.upper():
                    exp_type = exp_dict[entry]
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
        
        if isinstance(nuc_1, str):
            nuc_1 = nuc_1.strip()
        if isinstance(nuc_2, str):
            nuc_2 = nuc_2.strip()

        return nuc_1, nuc_2