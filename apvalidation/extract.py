import os

import nmrglue as ng
import pandas as pd


class Varian:
    """
    A class that contains the methods to extract parameters from Varian NMR data
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath):
        """
        Read the file to retrieve a dictionary of experiment parameters.

        :param filepath: a string formatted filepath to the procpar file
        :return: dictionary containing all parameters found in the procpar file OR
                'Empty File' if the file is empty
        """
        assert os.path.isfile(filepath)
        param_dict = ng.varian.read_procpar(filename=filepath)

        if len(param_dict) > 0:
            return param_dict
        else:
            return "Empty File"

    @staticmethod
    def find_params(param_dict):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2 frequency,
        solvent, and temperature.

        :param param_dict: a large dictionary containing all the parameters.
                            Probably returned from the read method.
        :return: dictionary containing only the preferred parameters
        """
        exp_dim = Varian.find_dim(param_dict)
        exp_type = Varian.find_exp_type(param_dict, exp_dim)
        exp_freq = Varian.find_freq(param_dict, exp_dim)
        exp_nuc1, exp_nuc2 = Varian.find_nuc(param_dict, exp_dim)
        exp_solv = Varian.find_solvent(param_dict)
        exp_temp = Varian.find_temp(param_dict)

        pref_params = {'experiment_type': exp_type,
                       'nuc_1': exp_nuc1,
                       'nuc_2': exp_nuc2,
                       'frequency': exp_freq,
                       'solvent': exp_solv,
                       'temperature': exp_temp
                       }
        return pref_params

    @staticmethod
    def find_temp(param_dict):
        temp_number = int(param_dict['temp']['values'][0])
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
        all_solvents = {
            **dict.fromkeys(['DMSO'], 'C2D6OS'),
            **dict.fromkeys(['CHLOROFORM', 'CHLOROFORM-D'], 'CDCL3'),
            **dict.fromkeys(['TETRAHYDROFURAN'], 'C2D6OS'),
            **dict.fromkeys(['DICHLOROMETHANE', 'DEUTERATED-DICHLOROMETHANE'], 'CD2CL2'),
            **dict.fromkeys(['ACETONE'], 'C3D6O'),
            **dict.fromkeys(['METHANOL', 'DEUTERATED-METHANOL', 'MEOD'], 'CD3OD'),
            **dict.fromkeys(['TOLUENE'], 'C7D8'),
            **dict.fromkeys(['HEAVY-WATER', 'OXIDANE', 'DEUTERIUM-OXIDE', 'H2O+D2O'], 'D2O'),
            **dict.fromkeys(['TRIFLUOROACETIC-ACID'], 'C2DF3O2'),
            **dict.fromkeys(['PYRIDINE'], 'C5D5N'),
            **dict.fromkeys(['METHYLENE-CHLORIDE', 'METHYLENE-CHLORI'], 'CD2CL2'),
            **dict.fromkeys(['ACETONITRILE'], 'C2D3N'),
            **dict.fromkeys(['BENZENE'], 'C6D6'),
            **dict.fromkeys(['URINE_KEY'], 'URINE'),
        }
        try:
            solv_str = param_dict['solvent']['values'][0]
            if solv_str is not None:
                solv_str = solv_str.upper()
        except KeyError:
            solv_str = None

        if solv_str in all_solvents.keys():
            exp_solv = all_solvents[solv_str].upper()
        elif solv_str in all_solvents.values():
            exp_solv = solv_str.upper()
        else:
            exp_solv = None
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
        try:
            exp_dim = param_dict['apptype']['values'][0][-2:]
        except KeyError:
            exp_dim = None

        if exp_dim in ['1D', '2D']:
            return exp_dim
        else:
            try:
                exp_dim = str(param_dict['procdim']['values'][0]) + "D"
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

        :param param_dict: dictionary containing all the parameter data
        :param exp_dim: dimension of the experiment
        :return: type of experiment in string. (1D experiments are not given a type)
        """
        exp_list = ['HSQCTOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'DOSY', 'ROESY', 'NOESY']
        # assume the experiment type is 1D and change from there 
        exp_type = ''
        if exp_dim == '2D':
            try:
                exp_type = param_dict['explist']['values'][0]
                print(f"THe explist value is: {exp_type}")
            except KeyError:
                exp_type = ""
            if exp_type == "":
                exp_type = param_dict['apptype']['values'][0]
            for type_str in exp_list:
                if type_str in exp_type.upper():
                    exp_type = type_str
                    return exp_type
        elif exp_dim == '1D':
            for type_str in exp_list:
                if type_str in exp_type.upper():
                    exp_type = type_str
                    break
            exp_type = f'1D-{exp_type}'
            return exp_type
        else:
            exp_type = None
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
        if exp_dim == '2D':
            freq1 = round(float(param_dict['reffrq']['values'][0]), 2)
            freq2 = round(float(param_dict['reffrq1']['values'][0]), 2)
            freq_val = (freq1, freq2)
        else:
            freq_val = round(float(param_dict['reffrq']['values'][0]), 2)
        return [freq_val]

    @staticmethod
    def find_nuc(param_dict, exp_dim):
        """
        Find the nuclei included in the Varian experiment.

        :param param_dict: dict containing param data
        :param exp_dim: dimension of the experiment
        :return: dimension of the experiment
        """
        explist_keyword_dict = {'PROTON': '1H', 'CARBON': '13C'}
        rewrite_nucs_dict = {"H1": "1H", "C13": "13C", "N15": "15N"}

        if exp_dim == '1D':
            try:
                atom_name = param_dict['explist']['values'][0]
                exp_nuc1 = explist_keyword_dict[atom_name]
                exp_nuc2 = None
            except KeyError:
                atom_name = param_dict['tn']['values'][0]
                exp_nuc1 = rewrite_nucs_dict[atom_name]
                exp_nuc2 = None

        elif exp_dim == '2D':
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
            "H1": [(0.99875, 1.00125), (3.971172962, 3.98111332), (9.864197531, 9.888888889),
                   (1.061652936, 1.064310391), (2.467572576, 2.473749228)],
            "C13": [(0.25025, 0.25275), (0.9950298211, 1.004970179),
                    (2.471604938, 2.496296296), (0.2660111613, 0.2686686155), (0.6182828907, 0.6244595429)],
            "N15": [(0.1, 0.1025), (0.3976143141, 0.407554672),
                    (0.987654321, 1.012345679), (0.1062981664, 0.1089556205), (0.2470660902, 0.2532427424)],
            "F": [(0.9395, 0.942), (3.735586481, 3.745526839),
                  (9.279012346, 9.303703704), (0.9986712729, 1.001328727), (2.321185917, 2.327362569)],
            "P": [(0.4035, 0.406), (1.604373757, 1.614314115),
                  (3.985185185, 4.009876543), (0.4289131012, 0.4315705554), (0.9969116739, 1.003088326)]
        }

        nuc_df = pd.DataFrame(nuc_dict)
        nuc_df.index = ["H1", "C13", "N15", "F", "P"]

        nuc_1 = param_dict['tn']['values'][0]
        nuc_2 = ''
        freq_ratio = float(param_dict['reffrq']['values'][0]) / float(param_dict['reffrq1']['values'][0])

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

        pref_params = {'experiment_type': exp_type, 'nuc_1': exp_nuc1, 'nuc_2': exp_nuc2, 'frequency': exp_freq,
                       'solvent': exp_solv, 'temperature': exp_temp}

        return pref_params

    @staticmethod
    def find_temp(param_dict):
        temp_number = param_dict['TE']
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

        all_solvents = {
            **dict.fromkeys(['DMSO'], 'C2D6OS'),
            **dict.fromkeys(['CHLOROFORM', 'CHLOROFORM-D'], 'CDCL3'),
            **dict.fromkeys(['TETRAHYDROFURAN'], 'C2D6OS'),
            **dict.fromkeys(['DICHLOROMETHANE', 'DEUTERATED-DICHLOROMETHANE'], 'CD2CL2'),
            **dict.fromkeys(['ACETONE'], 'C3D6O'),
            **dict.fromkeys(['METHANOL', 'DEUTERATED-METHANOL', 'MEOD'], 'CD3OD'),
            **dict.fromkeys(['TOLUENE'], 'C7D8'),
            **dict.fromkeys(['HEAVY-WATER', 'OXIDANE', 'DEUTERIUM-OXIDE', 'H2O+D2O'], 'D2O'),
            **dict.fromkeys(['TRIFLUOROACETIC-ACID'], 'C2DF3O2'),
            **dict.fromkeys(['PYRIDINE'], 'C5D5N'),
            **dict.fromkeys(['METHYLENE-CHLORIDE', 'METHYLENE-CHLORI'], 'CD2CL2'),
            **dict.fromkeys(['ACETONITRILE'], 'C2D3N'),
            **dict.fromkeys(['BENZENE'], 'C6D6'),
            **dict.fromkeys(['URINE_KEY'], 'URINE'),
        }

        solv_str = param_dict['SOLVENT']
        if solv_str.upper() in all_solvents.keys():
            exp_solv = all_solvents[solv_str.upper()]
        elif solv_str.upper() in all_solvents.values():
            exp_solv = solv_str
        else:
            exp_solv = None
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
            exp_dim = '2D'
        elif len(param_dict_list) == 1:
            exp_dim = '1D'

        # exp_2d_list = ['HSQCTOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'ROESY', 'NOESY', 'DOSY']
        # exp_str = param_dict['EXP']
        # exp_dim = exp_str[0:2].upper()

        # if exp_dim not in ['1D', '2D']:
        #     for entry in exp_2d_list:
        #         if entry in exp_str.upper():
        #             exp_dim = '2D'
        #             return exp_dim
        #     exp_dim = '1D'
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
        possible_exp_str_1 = param_dict['EXP']
        possible_exp_str_2 = param_dict['PULPROG']

        exp_2d_list = ['HSQCTOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'ROESY', 'NOESY', 'DOSY']
        # assume exp_type is 1D and change from there 
        exp_type = ''
        if exp_dim == '2D':
            for entry in exp_2d_list:
                if entry in possible_exp_str_1.upper() or entry in possible_exp_str_2.upper():
                    exp_type = entry
                    return exp_type
            exp_type = None
            return exp_type
        else:
            print("This experiment is 1D!")
            for entry in exp_2d_list:
                if entry in possible_exp_str_1 or entry in possible_exp_str_2:
                    exp_type = entry
                    break
            exp_type = f'1D-{exp_type}'
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
        if exp_dim == '1D':
            freq_val = round(float(param_dict['SFO1']), 2)
            return [freq_val]
        else:
            freq1 = round(float(param_dict['SFO1']), 2)
            freq2 = round(float(param_dict['SFO2']), 2)
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
        if exp_dim == '2D':
            param_dict_2D = param_dict_list[1]

        exp_nuc2 = None
        if exp_dim == '1D':
            exp_nuc1 = param_dict_1D['NUC1']
        elif exp_dim == '2D':
            exp_nuc1 = param_dict_1D['NUC1']
            exp_nuc2 = param_dict_2D['NUC1']
        else:
            exp_nuc1 = None

        return exp_nuc1, exp_nuc2
        
        # exp_nuc2 = None
        # if exp_dim == '1D':
        #     exp_nuc1 = param_dict['NUC1']
        # elif param_dict['NUC2'] != 'off':
        #     exp_nuc1 = param_dict['NUC1']
        #     exp_nuc2 = param_dict['NUC2']
        # else:
        #     exp_nuc1 = param_dict['NUC1']
        #     exp_nuc2 = param_dict['NUC1']
        # return exp_nuc1, exp_nuc2


class Jcampdx:
    """
    A class containing the methods to help with the extraction of
    experiment parameters from NMR data from jdx format.
    """

    def __init__(self):
        pass

    @staticmethod
    def read(filepath):
        """
        Given the raw file path to a parameter file, read the file
        and return a python dictionary with all the parameters.

        :param filepath: string formatted filepath to the acqu file
        :return: dictionary containing all parameters found in the acqu file
        """
        assert os.path.isfile(filepath)
        param_dict = ng.jcampdx.read(filename=filepath)
        try:
            return param_dict[0]['_datatype_NMRSPECTRUM'][0]
        except KeyError:
            return param_dict[0]

    @staticmethod
    def find_params(param_dict):
        """
        Searches a dictionary of all the parameters from a given experiment to find specific parameters needed
        to fill the database. The parameters needed are experiment_type, nucleus 1 and 2 frequency,
        solvent, and temperature.

        :param param_dict: dictionary containing all the parameters retrieved from the file
        :return: dictionary containing only those parameters that are preferred
        """
        exp_dim = Jcampdx.find_dim(param_dict)
        exp_freq = Jcampdx.find_freq(param_dict, exp_dim)
        exp_nuc_1, exp_nuc_2 = Jcampdx.find_nuc(param_dict, exp_dim)
        exp_type = Jcampdx.find_exp_type(param_dict, exp_dim)
        exp_solv = Jcampdx.find_solvent(param_dict)
        exp_temp = Jcampdx.find_temp(param_dict)

        pref_params = {'experiment_type': exp_type, 'nuc_1': exp_nuc_1, 'nuc_2': exp_nuc_2, 'frequency': exp_freq,
                       'solvent': exp_solv.upper(), 'temperature': exp_temp}

        return pref_params

    @staticmethod
    def find_temp(param_dict):
        temp_number = int(param_dict['$TEMPSET'][0])
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
        all_solvents = {
            **dict.fromkeys(['DMSO'], 'C2D6OS'),
            **dict.fromkeys(['CHLOROFORM', 'CHLOROFORM-D'], 'CDCL3'),
            **dict.fromkeys(['TETRAHYDROFURAN'], 'C2D6OS'),
            **dict.fromkeys(['DICHLOROMETHANE', 'DEUTERATED-DICHLOROMETHANE'], 'CD2CL2'),
            **dict.fromkeys(['ACETONE'], 'C3D6O'),
            **dict.fromkeys(['METHANOL', 'DEUTERATED-METHANOL', 'MEOD'], 'CD3OD'),
            **dict.fromkeys(['TOLUENE'], 'C7D8'),
            **dict.fromkeys(['HEAVY-WATER', 'OXIDANE', 'DEUTERIUM-OXIDE', 'H2O+D2O'], 'D2O'),
            **dict.fromkeys(['TRIFLUOROACETIC-ACID'], 'C2DF3O2'),
            **dict.fromkeys(['PYRIDINE'], 'C5D5N'),
            **dict.fromkeys(['METHYLENE-CHLORIDE', 'METHYLENE-CHLORI'], 'CD2CL2'),
            **dict.fromkeys(['ACETONITRILE'], 'C2D3N'),
            **dict.fromkeys(['BENZENE'], 'C6D6'),
            **dict.fromkeys(['URINE_KEY'], 'URINE'),
        }
        solv_str = param_dict['$SOLVENT'][0]
        if solv_str in all_solvents.keys():
            exp_solv = all_solvents[solv_str]
        elif solv_str.upper() in all_solvents.values():
            exp_solv = solv_str
        else:
            exp_solv = None
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
        exp_dim = param_dict['NUMDIM'][0] + 'D'
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
        exp_str = param_dict['.PULSESEQUENCE'][0].upper()

        exp_2d_list = ['HSQC-TOCSY', 'COSY', 'HSQC', 'HMQC', 'HMBC', 'TOCSY', 'ROESY', 'NOESY', 'DOSY']
        exp_type = None

        if exp_dim == '2D':
            for entry in exp_2d_list:
                if entry in exp_str.upper():
                    exp_type = entry
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
        if exp_dim == '2D':
            freq1 = round(float(param_dict['$XFREQ'][0]), 2)
            freq2 = round(float(param_dict['$YFREQ'][0]), 2)
            freq = (freq1, freq2)
            return freq
        else:
            freq = round(float(param_dict['$XFREQ'][0]), 2)

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
        if exp_dim == '2D':
            nuc_1 = param_dict['.NUCLEUS'][0].split(', ')[1]
            nuc_2 = param_dict['.NUCLEUS'][0].split(', ')[0]
        else:
            nuc_1 = param_dict['.OBSERVENUCLEUS'][0][1:]
            nuc_2 = None

        return nuc_1, nuc_2
