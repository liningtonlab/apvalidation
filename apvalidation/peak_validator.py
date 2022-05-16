import re
import statistics as stats
import warnings

from rdkit import Chem
from rdkit.Chem import rdqueries, rdmolops



"""
Create three classes, one to validate Proton peaks and one to validate Carbon peaks
"""

class InvalidCharacters(Exception):
    pass

class NoSplit(Exception):
    pass

class InvalidValueType(Exception):
    pass

class InvalidAtomNumber(Exception):
    pass

class WarnBadRange(Exception):
    pass

class ErrorBadRange(Exception):
    pass

class Validate:

    @staticmethod
    def check_valid_characters(text_block):

        clean_text = text_block.replace(" ", "")
        accepted_pattern = re.compile(r"^[\d . - ( ) , ; \u002D \u05BE \u1806 \u2010 \u2011 \u2012 \u2013 \\\
                                         \u2014 \u2015 \u207B \u208B \u2212 \uFE58 \uFE63 \uFF0D]+$", re.UNICODE)
        if re.search(accepted_pattern, clean_text):
            return "Valid String"
        else:
            raise InvalidCharacters
    
    @staticmethod
    def parse_text_to_list(valid_text):
        try:
            split_text = valid_text.split(",")
            if len(split_text) == 1:
                raise NoSplit
        except NoSplit:
            try:
                split_text = valid_text.split(";")
                if len(split_text) == 1:
                    raise NoSplit
            except NoSplit:
                try:
                    split_text = valid_text.split("\n")
                    if len(split_text) == 1:
                        raise NoSplit
                except NoSplit:
                    try:
                        split_text = valid_text.split("\t")
                        if len(split_text) == 1:
                            raise NoSplit
                    except NoSplit:
                        raise NoSplit
        return split_text

    @staticmethod
    def is_valid_range(value):
        value = value.replace(" ", "")
        range_pattern = re.compile(r"^\({0,1}\d{1,3}.{0,1}\d{0,4} {0,5}[\u002D\u05BE\u1806\u2010\u2011\u2012\u2013 \\\
                                    \u2014\u2015\u207B\u208B\u2212\uFE58\uFE63\uFF0D]{1} {0,5}\d{1,3}.{0,1}\d{0,4}\){0,1}$", re.UNICODE)
        re_result = re.search(range_pattern, value)         
        if re_result:
            return True
        else:
            return False

    @staticmethod
    def check_data_type(value_list):

        value_check_list = [False]*len(value_list)
        value_list = [value.replace(" ","") for value in value_list]

        for index, value in enumerate(value_list):
            if Validate.is_valid_range(value):
                value_check_list[index] = True
                print(value_check_list)
                return value_check_list
            else:
                try:
                    if isinstance(float(value.replace(" ", "")), float):
                        value_check_list[index] = True
                        print(value_check_list)
                        return value_check_list
                except ValueError:
                    return value_check_list

    @staticmethod
    def check_number_atoms(value_list, atom_type, smiles):
        atoms = {"H": 1, "C": 6}
        atom_num = atoms[atom_type]
        mol = Chem.MolFromSmiles(smiles)
        mol = rdmolops.AddHs(mol)

        atom_query = rdqueries.AtomNumEqualsQueryAtom(atom_num)
        number_of_atoms = len(mol.GetAtomsMatchingQuery(atom_query))

        if 0 < len(value_list) <= number_of_atoms:
            return True
        else:
            raise InvalidAtomNumber

    @staticmethod
    def check_value_ranges(value_list, atom_type):

        unchecked_values = []
        for value in value_list:
            is_range = Validate.is_valid_range(value)
            if is_range:
                split_value = re.split("[\u002D\u05BE\u1806\u2010\u2011\u2012\u2013\u2014\u2015\u207B\u208B\u2212\uFE58\uFE63\uFF0D]",value)
                split_value = [unchecked_values.append(float(re.sub("[ ()]", "", sub_str))) for sub_str in split_value]
            else:
                unchecked_values.append(float(value))

        for value in unchecked_values:
            try:
                if atom_type == "H":
                    if value < -2 or value > 20:
                        raise ErrorBadRange
                    elif value < -0.5 or value > 14:
                        warnings.warn(f"Warning: {value} is out of typical bounds")
            except ErrorBadRange:
                return f"Error: {value} out of bounds"
            try:
                if atom_type == "C":
                    if value < -20 or value >250:
                        raise ErrorBadRange
                    elif value < -10 or value >230:
                        warnings.warn(f"Warning: {value} is out of typical bounds")
            except ErrorBadRange:
                return f"Error: {value} out of bounds"
        
        return True
        

    @staticmethod
    def validate(H_text_block, C_text_block, smiles):
        # Check that characters in the text block are valid 
        try:
            Validate.check_valid_characters(H_text_block)
        except InvalidCharacters:
            return "Error Invalid Characters in H List: Please make sure that only contains the following allowed characters 0-9 , . - ; ()"
        try:
            Validate.check_valid_characters(C_text_block)
        except InvalidCharacters:
            return "Error Invalid Characters in C List: Please make sure that only contains the following allowed characters 0-9 , . - ; ()"

        # Parse the text blocks into lists based on the seporators
        try:
            H_list = Validate.parse_text_to_list(H_text_block)
        except NoSplit:
            return "Failed to split H list, please check your seporators."
        try:
            C_list = Validate.parse_text_to_list(C_text_block)
        except NoSplit:
            return "Failed to split C list, please check your seporators."
        
        # Check if each element in the parsed lists are either floats or ranges
        try:
            error_list = Validate.check_data_type(H_list)
            if sum(error_list) != len(H_list):
                error_value_list = []
                for index, b_val in enumerate(error_list):
                    if b_val is False:
                        error_value_list.append(H_list[index])
                    else:
                        continue
                raise InvalidValueType

        except InvalidValueType:
            return f"Error: The following values are of invalid data type {error_value_list}"
        try:
            if sum(Validate.check_data_type(C_list)) != len(H_list):
                error_value_list = []
                for index, b_val in enumerate(H_list):
                    if b_val is False:
                        error_value_list.append(H_list[index])
                    else:
                        continue
                raise InvalidValueType
        except InvalidValueType:
            return f"Error: The following values are of invalid data type {error_value_list}"

        try:
            Validate.check_number_atoms(H_list, "H", smiles)
        except InvalidAtomNumber:
            return "Error: Invalid number of H atoms in the peak list"
        try:
            Validate.check_number_atoms(C_list, "C", smiles)
        except InvalidAtomNumber:
            return "Error: Invalid number of C atoms in the peak list"

        try:
            output = Validate.check_value_ranges(H_list, "H")
        except ErrorBadRange:
            return output
        try:
            output = Validate.check_value_ranges(C_list, "C")
        except ErrorBadRange:
            return output

        return "Both lists are valid"

        































# class Validate:
#     """
#     This will be the class that contains the validate_peak_lists function.
#     This function will be the one called by the user in order to validate the peaks.
#     """
#     @staticmethod
#     def check_lengths(peak_dict):
#         h_peaks = peak_dict['H']
#         c_peaks = peak_dict['C']
#         if len(h_peaks) != len(c_peaks):
#             return False
#         return True

#     @staticmethod
#     def validate_peak_lists(peak_dict, smiles):

#         return_bool = True
#         return_messages = []

#         if Validate.check_lengths(peak_dict) is False:
#             return_bool = False
#             return_messages.append("Length of carbon and proton list do not match.")

#         # currently returns the number of unique protons doesn't work properly I need to know what to check
#         num_protons = Proton.check_atom_counts(smiles)
#         print(num_protons)
#         if num_protons is False:
#             return_bool = False
#             return_messages.append("Invalid number of unique protons.")

#         # currently returns the number of carbons doesn't work properly I need to know what to check
#         num_carbons = Carbon.check_atom_counts(smiles)
#         print(num_carbons)
#         if num_carbons is False:
#             return_bool = False
#             return_messages.append("Invalid number of unique carbons.")  

#         h_list = peak_dict['H']
#         if Proton.check_value_ranges(h_list) is False:
#             return_messages.append("WARNING: Proton values out of range.")

#         c_list = peak_dict['C']
#         if Carbon.check_values_ranges(c_list) is False:
#             return_messages.append("WARNING: Carbon values out of range.")

#         return return_bool, return_messages

# class Proton:

#     @staticmethod
#     def check_value_ranges(h_list):
#         mean = stats.mean(h_list)
#         median = stats.median(h_list)
#         if mean > 10:
#             return False
#         if median > 10:
#             return False
        
#         for item in h_list:
#             if item > 20:
#                 return False
#         return True

#     def check_atom_counts(smiles):
#         mol = Chem.MolFromSmiles(smiles)
#         mol_h = Chem.AddHs(mol)
#         ordered_atoms = Chem.CanonicalRankAtoms(mol_h, breakTies=False)
#         unique_hs = set()
#         all_hs = []
#         for atom, sym_classes in zip(mol_h.GetAtoms(), ordered_atoms):
#             if atom.GetAtomicNum() == 1:
#                 all_hs.append(atom)
#                 unique_hs.add(sym_classes)

#         return len(unique_hs)


# class Carbon:

#     @staticmethod
#     def check_values_ranges(c_list):
#         mean = stats.mean(c_list)
#         median = stats.median(c_list)
#         if mean < 10:
#             return False
#         if median < 10:
#             return False
        
#         for item in c_list:
#             if item > 250:
#                 return False
#         return True

#     @staticmethod
#     def check_atom_counts(smiles):
#         mol = Chem.MolFromSmiles(smiles)
#         all_cs = []
#         for atom in mol.GetAtoms():
#             if atom.GetAtomicNum() == 6:
#                 all_cs.append(atom)

#         return len(all_cs)



    