import os
import nmrglue as ng
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
    def find_params(
        param_dict_list,
        json_nmr_data_dict:dict = {},
    ):
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

        exp_dim = JEOL.find_dim(param_dict, json_nmr_data_dict)
        exp_freq = JEOL.find_freq(param_dict, exp_dim, json_nmr_data_dict)
        exp_nuc_1, exp_nuc_2 = JEOL.find_nuc(param_dict, exp_dim, json_nmr_data_dict)
        exp_type = JEOL.find_exp_type(param_dict, exp_dim, json_nmr_data_dict)
        exp_solv = JEOL.find_solvent(param_dict, json_nmr_data_dict)
        exp_temp = JEOL.find_temp(param_dict, json_nmr_data_dict)

        pref_params = {
            "experiment_type": exp_type,
            "nuc_1": exp_nuc_1,
            "nuc_2": exp_nuc_2,
            "frequency": exp_freq,
            "solvent": exp_solv,
            "temperature": exp_temp,
        }

        return pref_params

    @staticmethod
    def find_temp(param_dict, json_nmr_data_dict={}):
        """
        Helper function to determine the experimental temperature.
        Checks both param_dict and json_nmr_data_dict as fallback.

        :param param_dict: Dictionary containing experiment parameters.
        :param json_nmr_data_dict: Dictionary containing JSON NMR data (optional).
        :return: Temperature in Kelvin as integer, or None if not found.
        """
        temp_number = None
        
        # First attempt: param_dict
        try:
            temp_number = float(param_dict["$TEMPSET"][0])
        except:
            try:
                temp_number = float(param_dict["$TEMPGET"][0])
            except:
                pass  # proceed to JSON fallback

        # Fallback: json_nmr_data_dict
        if temp_number is None and json_nmr_data_dict:
            try:
                json_temp = json_nmr_data_dict.get("SpecInfo", {}).get("Temperature", None)
                if json_temp is not None:
                    temp_number = int(float(json_temp))
            except (TypeError, ValueError):
                pass

        # If temp found, standardize to Kelvin
        if temp_number is not None:
            # if under 250 assume C. Otherwise assume K and adjust.
            return temp_number if temp_number >= 250 else temp_number + 273

        return None  # Failed to detect temperature

    @staticmethod
    def find_solvent(param_dict, json_nmr_data_dict={}):
        """
        Helper function. Find the solvent involved in the experiment. This is done by checking possible
        english names for the solvent against their chemical formulas.

        :param param_dict: Dictionary containing all the parameters for the experiment.
        :return: The solvent in string format.
        """
        exp_solv = "FAILED_TO_DETECT"

        # First attempt: param_dict
        solv_list = param_dict.get("$SOLVENT", [None])
        solv_str = solv_list[0] if solv_list else None
        if solv_str:
            solv_upper = solv_str.upper()
            if solv_upper in all_solvents:
                exp_solv = all_solvents[solv_upper]
            elif solv_str in all_solvents.values():
                exp_solv = solv_str

        # Fallback: json_nmr_data_dict
        if exp_solv == "FAILED_TO_DETECT" and json_nmr_data_dict:
            json_solv = json_nmr_data_dict.get("SpecInfo", {}).get("Solvent", "")
            if json_solv:
                json_solv_upper = json_solv.upper()
                if json_solv_upper in all_solvents:
                    exp_solv = all_solvents[json_solv_upper]
                elif json_solv in all_solvents.values():
                    exp_solv = json_solv

        return exp_solv

    @staticmethod
    def find_dim(param_dict, json_nmr_data_dict={}):
        """
        Helper function.
        Find the dimension of the experiment given a dictionary of parameters.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :return: dimension of the experiment.
        """
        exp_dim = None

        # Try to get from param_dict
        try:
            exp_dim = param_dict["NUMDIM"][0] + "D"
        except:
            pass
        if exp_dim == None:
            try:
                exp_dim = param_dict["$DIMENSIONS"][0] + "D"
            except:
                pass
        
        # FALLBACK: check json_nmr_data_dict
        if not exp_dim and json_nmr_data_dict:
            try:
                n_dim = json_nmr_data_dict.get("nDim", "")
                if n_dim and n_dim <= 2:
                    exp_dim = str(n_dim) + "D"
            except:
                pass

        return exp_dim

    @staticmethod
    def find_exp_type(param_dict, exp_dim, json_nmr_data_dict={}):
        """
        Helper function.
        Determine the type of experiment that was conducted. Ex. COSY, HSQC etc...

        :param param_dict: dictionary containing all the parameter data
        :param exp_dim: dimension of the experiment
        :param json_nmr_data_dict: dict containing JSON NMR data (optional)
        :return: experiment type as string ("1D", "COSY", etc.), or 'FAILED_TO_DETECT' if unknown
        """
        exp_type = "FAILED_TO_DETECT"

        def standardize_exp_type(exp_string):
            """Standardize the experiment string using exp_dict mappings."""
            if not exp_string:
                return None
            exp_string = exp_string.upper()
            for key, val in exp_dict.items():
                if key.upper() in exp_string:
                    return val
            return None

        # For 2D experiments: first try param_dict
        if exp_dim == "2D":
            pulse_list = param_dict.get(".PULSESEQUENCE", [None])
            pulse_seq = pulse_list[0] if pulse_list else None
            if pulse_seq:
                exp_type_from_pulse = standardize_exp_type(pulse_seq)
                if exp_type_from_pulse:
                    return exp_type_from_pulse

            # Fallback: try JSON NMR data
            json_exp_type = json_nmr_data_dict.get("SpecInfo", {}).get("ExperimentType.str", "")
            exp_type_from_json = standardize_exp_type(json_exp_type)
            if exp_type_from_json:
                return exp_type_from_json

            return exp_type  # Failed to detect

        # For 1D experiments
        return "1D"

    @staticmethod
    def find_freq(param_dict, exp_dim, json_nmr_data_dict={}):
        """
        Helper function.
        Find the frequency parameter given a dictionary of parameters and the
        dimension of the experiment.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :param exp_dim: the dimension of the experiment
        :return: a single float or a tuple of floats depending on the dimension
        """
        
        # print("param_dict is")
        # print(param_dict.keys())
        # for key, value in param_dict.items():
        #     if (key != "$PARAMETERFILE") and (key != "DATATABLE"):
        #         # print(f"{key} - {value}")
        #         pass
        def get_safe_freq(key):
            try:
                return round(float(param_dict[key][0]), 9)
            except (KeyError, IndexError, ValueError, TypeError):
                return None

        def get_json_freqs(n=2):
            freq_list = json_nmr_data_dict.get("SpecInfo", {}).get("SpectrometerFrequencies", [])
            return freq_list[:n] if freq_list else []

        if exp_dim == "2D":
            freq1 = get_safe_freq("$XFREQ")
            freq2 = get_safe_freq("$YFREQ")

            if freq1 is not None and freq2 is not None:
                return (freq1, freq2)

            # Try json_nmr_data_dict as fallback
            exp_freq = get_json_freqs(2)
            if exp_freq:
                return exp_freq

            return [freq1, freq2]  # Will contain None if missing

        else:  # "1D" or other
            freq = get_safe_freq("$XFREQ")

            if freq is None:
                exp_freq = get_json_freqs(1)
                return exp_freq[0] if exp_freq else None

            return [freq]

    @staticmethod
    def find_nuc(param_dict, exp_dim, json_nmr_data_dict={}):
        """
         Helper function.
         Determine the nuclei that are used in the experiment.

        :param param_dict: a large dictionary containing all the parameters.
                             Probably returned from the read method.
        :param exp_dim: the dimension of the experiment
        :return: tuple of (nuc_1, nuc_2)
        """
        nuc_1, nuc_2 = None, None

        # Primary method: param_dict
        try:
            if exp_dim == "2D":
                nuc_values = param_dict.get(".NUCLEUS", [None])[0]
                if nuc_values:
                    nuc_values = nuc_values.split(",")
                    nuc_1 = nuc_values[1].strip() if len(nuc_values) > 1 else None
                    nuc_2 = nuc_values[0].strip() if len(nuc_values) > 0 else None
            else:  # 1D
                nuc_1 = param_dict.get(".OBSERVENUCLEUS", [None])[0]
                if nuc_1:
                    nuc_1 = nuc_1[1:].strip()  # Skip first char
        except Exception:
            pass

        # Normalize param_dict values
        def normalize_nucleus(nuc):
            if isinstance(nuc, str):
                nuc = nuc.strip()
                if "proton" in nuc.lower():
                    return "1H"
                if "carbon13" in nuc.lower():
                    return "13C"
            return nuc

        nuc_1 = normalize_nucleus(nuc_1)
        nuc_2 = normalize_nucleus(nuc_2)

        # Fallback: use json_nmr_data_dict if needed
        if not nuc_1 or (exp_dim == "2D" and not nuc_2):
            nucleides = json_nmr_data_dict.get("SpecInfo", {}).get("Nucleides", [])

            # Filter out empty isotopes
            valid_nucleides = [nuc for nuc in nucleides if nuc.get("Isotope", 0) != 0]

            if valid_nucleides:
                if not nuc_1 and len(valid_nucleides) > 0:
                    nuc_1 = f"{valid_nucleides[0]['Isotope']}{valid_nucleides[0]['Name']}"
                if exp_dim == "2D" and not nuc_2 and len(valid_nucleides) > 1:
                    nuc_2 = f"{valid_nucleides[1]['Isotope']}{valid_nucleides[1]['Name']}"

        return nuc_1, nuc_2