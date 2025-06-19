import os
import re
import nmrglue as ng
import json
import traceback

from apvalidation.extract.extract_varian import Varian
from apvalidation.extract.extract_bruker import Bruker
from apvalidation.extract.extract_jeol import JEOL

current_dir = os.path.dirname(__file__)
one_level_up = os.path.dirname(current_dir)
submodule_dir = os.path.join(one_level_up, 'npmrd_data_exchange')

experiment_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'experiment_standardizer.json')
with open(experiment_standardizer_path, 'r') as file:
    exp_dict = json.load(file) 

solvent_standardizer_path = os.path.join(submodule_dir, 'standardization_files', 'solvent_standardizer.json')
with open(solvent_standardizer_path, 'r') as file:
    all_solvents = json.load(file) 

class Jcampdx:
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
    def read(filepath):
        """
        Given the raw file path to a parameter file, read the file
        and return a python dictionary with all the parameters and
        optionally embedded JSON NMR data.
        """
        assert os.path.isfile(filepath)

        read_dict, read_np_array = ng.jcampdx.read(filepath)

        # Defensive helper to flatten nested lists/tuples if necessary
        def flatten_first_entry(obj):
            while isinstance(obj, (list, tuple)) and len(obj) > 0:
                obj = obj[0]
            return obj

        flat_read_dict = flatten_first_entry(read_dict)

        # Try extracting param_dict following your old key hierarchy
        param_dict = {}
        try:
            param_dict = flat_read_dict["_datatype_NMRSPECTRUM"]
        except (KeyError, TypeError):
            try:
                param_dict = flat_read_dict["_datatype_NDNMRSPECTRUM"]
            except (KeyError, TypeError):
                param_dict = flat_read_dict

        # If param_dict is still nested list/tuple, flatten again
        param_dict = flatten_first_entry(param_dict)

        # Attempt to extract JSON NMR data from param_dict
        json_nmr_data_dict = {}
        try:
            # Try to get JSON NMR string
            json_nmr_str = read_dict['_datatype_LINK'][0]['$JASONNMRDATA'][0]

            # Extract JSON from string start
            json_nmr_start = json_nmr_str.find('{')
            json_nmr_str_clean = json_nmr_str[json_nmr_start:]
            json_nmr_data_dict = json.loads(json_nmr_str_clean)
            json_nmr_data_dict = flatten_first_entry(json_nmr_data_dict)
            print("Successfully extracted embedded JSON NMR data")

        except (KeyError, IndexError, TypeError, json.JSONDecodeError) as e:
            json_nmr_data_dict = {}

        print("json_nmr_data_dict.keys() is:", json_nmr_data_dict.keys())
        print("param_dict keys are:", getattr(param_dict, 'keys', lambda: [])())

        return param_dict, json_nmr_data_dict

    @staticmethod
    def find_manuf(param_dict, json_nmr_data_dict={}):
        """
        find the manufacturer that produced the data file that is being inspected.

        :param param_dict_list: a list of dictionaries read in from the read function
        :return: the name of the manufacturer
        """
        
        # Path lcoations to check for manufacturer name
        try:
            manuf_name = param_dict["_datatype_LINK"][0]["$ORIGINALFORMAT"][0]
        except (KeyError, TypeError):
            manuf_name = "Not found"
            try:
                manuf_name = param_dict["$ORIGINALFORMAT"][0]
            except (KeyError, TypeError):
                manuf_name = "Not found"
                try:
                    manuf_name = param_dict["ORIGIN"][0]
                except (KeyError, TypeError):
                    try:
                        if json_nmr_data_dict:
                            manuf_name = json_nmr_data_dict['SpecInfo']['OrigFileFormat.str']
                    except:
                        manuf_name = "Not found"

        # Try to match manufacturer name
        manuf_name_lower = manuf_name.lower()
        manufacturer_keywords = {
            "varian": "Varian",
            "bruker": "Bruker",
            "jeol": "JEOL",
            "delta": "JEOL"
        }
        for keyword, label in manufacturer_keywords.items():
            if keyword in manuf_name_lower:
                return label
        
        return "Not found"

    @staticmethod
    def find_params(
        param_dict,
        json_nmr_data_dict: dict = {},
        manuf: str = None
    ):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2, frequency,
        solvent, and temperature.

        :param param_dict: dictionary containing all the parameters retrieved from the file
        :return: dictionary containing only those parameters that are preferred
        """

        output_list = []
        if not manuf:
            manuf = Jcampdx.find_manuf(param_dict)

        if manuf == "Varian":
            varian_structured_dict_list = Jcampdx.format_varian(param_dict)
            output_list.append(Varian.find_params(varian_structured_dict_list))

        elif manuf == "Bruker":
            bruker_structured_dict_list = Jcampdx.format_bruker(param_dict) # THIS IS WHAT'S BREAKING THE .DX (probably)!!!!!!!
            output_list.append(Bruker.find_params(bruker_structured_dict_list))

        # Assume jeol as fallback if we can't find manufacturer
        elif manuf == "JEOL" or manuf == "Not found":
            jeol_structured_dict_list = Jcampdx.get_jeol_structured_dict_list(param_dict)
            output_list.append(JEOL.find_params(
                jeol_structured_dict_list,
                json_nmr_data_dict=json_nmr_data_dict
            ))

            # As final fallback try to extract from other other manufacturers
            if not output_list and manuf == "Not found":
                try:
                    varian_structured_dict_list = Jcampdx.format_varian(param_dict)
                    output_list.append(Varian.find_params(varian_structured_dict_list))
                except:
                    pass
                if not output_list:
                    try:
                        bruker_structured_dict_list = Jcampdx.format_bruker(param_dict)
                        output_list.append(Bruker.find_params(bruker_structured_dict_list))
                    except:
                        pass

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
            line_list = jdx_read_output["_datatype_LINK"][0]["_comments"]
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

        param_dict = read_jdx_output

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
    def get_jeol_structured_dict_list(jdx_read_output):
        """
        Take the output produced from Jcamp read and format it to be passed into the Varian Class methods.

        :param jdx_read_output: a nested list object, the output from the Jcamp read method.
        :return: list of dictionaries, these are formatted for the Varian Class methods.
        """
        try:
            param_dict = jdx_read_output["_datatype_LINK"][0]
        except:
            try:
                param_dict = jdx_read_output["_datatype_LINK"]
            except:
                try:
                    param_dict = jdx_read_output
                except:
                    param_dict = {}

        # re-name and format frequency keys so they match the delta version of JEOL data.
        try:
            freq_list = param_dict[".OBSERVEFREQUENCY"][1]
            freq_list = freq_list.split(",")
            param_dict["$XFREQ"] = [freq_list[0]]
            param_dict["$YFREQ"] = [freq_list[1]]
        except:
            try:
                freq_list = param_dict[".OBSERVEFREQUENCY"][0]
                param_dict["$XFREQ"] = [freq_list]
            except:
                param_dict.setdefault("$XFREQ", [None])
        
        # add SOLVENT keys
        try:
            param_dict["$SOLVENT"] = param_dict[".SOLVENTNAME"]
        except:
            param_dict.setdefault("$SOLVENT", [None])

        # Add Temperature values
        try:
            param_dict["$TEMPSET"] = param_dict[".TEMPSET"]
        except:
            param_dict.setdefault("$TEMPSET", [None])
        try:
            param_dict["$TEMPGET"] = param_dict[".TEMPGET"]
        except:
            param_dict.setdefault("$TEMPGET", [None])

        return [jdx_read_output]