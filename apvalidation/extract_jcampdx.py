import os
import re
import nmrglue as ng
import json

from apvalidation.extract_varian import Varian
from apvalidation.extract_bruker import Bruker
from apvalidation.extract_jeol import JEOL

current_dir = os.path.dirname(__file__)
submodule_dir = os.path.join(current_dir, 'npmrd_data_exchange')

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
        manuf = Jcampdx.find_manuf(param_dict)

        if manuf == "Varian":
            varian_structured_dict_list = Jcampdx.format_varian(param_dict)
            output_list.append(Varian.find_params(varian_structured_dict_list))

        elif manuf == "Bruker":
            bruker_structured_dict_list = Jcampdx.format_bruker(param_dict)
            output_list.append(Bruker.find_params(bruker_structured_dict_list))

        elif manuf == "JEOL":
            jeol_structured_dict_list = Jcampdx.format_jeol_combined(param_dict)
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