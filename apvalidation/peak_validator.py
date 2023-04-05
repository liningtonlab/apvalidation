import re
import stat
import statistics as stats
import warnings

from rdkit import Chem
from rdkit.Chem import rdqueries, rdmolops



"""
These are custom exception classes to help annotate what error we are handling.
"""

# class EmptyList(Exception):
#     def __init__(self):
#         self.error_type = "error"
#         self.error = "EmptyList"

class InvalidCharacters(Exception):
    def __init__(self):
        self.error_type = "error"
        self.error = "InvalidCharacters"

class NoSplit(Exception):
    pass

class InvalidValueType(Exception):
    pass

class InvalidAtomNumber(Exception):
    pass

class WarnBadRange(Exception):
    def __init__(self, bad_value):
        self.error_type = "warning"
        self.bad_value = bad_value

class ErrorBadRange(Exception):
    def __init__(self, bad_value):
        self.error_type = "error"
        self.bad_value = bad_value



class Validate:
    @staticmethod
    def check_valid_characters(text_block):
        """
            Check that only accepted characters are present in the text taken from the input text box. 

            :param text_block: The text taken from the user input
            :return: If valid, return a string confirming valid chars. If invalid raise InvalidCharacters exception.
        """
        if not text_block:
            return "Valid String"

        clean_text = text_block.replace(" ", "")
        accepted_pattern = re.compile(r"^[\d . - ( ) , ; \u002D \u05BE \u1806 \u2010 \u2011 \u2012 \u2013 \\\
                                         \u2014 \u2015 \u207B \u208B \u2212 \uFE58 \uFE63 \uFF0D \t \n]+$", re.UNICODE)
        if re.search(accepted_pattern, clean_text):
            return "Valid String"
        else:
            raise InvalidCharacters
    
    @staticmethod
    def parse_text_to_list(valid_text):
        """
            Split the text from the peak list text box into a list of values and ranges. There are multiple
            different list formats to try.

            :param valid_text: The text to split
            :return: If works, return the split list. If split does not work, raise a NoSplit exception.
        """
        char_list = [",", ";", "\n", "\t", "\\t", "    "]
        split = False
        for split_char in char_list:
            split_text = valid_text.split(split_char)
            if split_text[-1] == "":
                split_text.pop()
            if len(split_text) == 1:
                continue
            else:
                split = True
                break
 
        if split == False:
            raise NoSplit
        return split_text

    @staticmethod
    def is_valid_range(value):
        """
            Check if the value is of a valid range format. The format is specified by the regex expression where the unicode characters
            are different representations of dashes.

            :param value: a value from the peak list that may be a valid range
            :return: True or False depending on if the value is a valid range.
        """
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
        """
            Check the datatype for each value in the peak value list. Ensure that the type is either 
            a float value or a properly formatted range of values, ex. (13.2 - 14.5).

            :param value_list: a list of the values and ranges in the peak list.
            :return: A list of booleans which represents which values in value_list are valid datatypes and which are not.
                    Ex. If all are valid return [True, True,...], If first is invalid and rest are valid return [False, True, True, ...]
        """
        value_check_list = [False]*len(value_list)
        # value_list = [value.replace(" ","") for value in value_list]

        for index, value in enumerate(value_list):

            value = value.replace(" ", "")

            if Validate.is_valid_range(value):
                value_check_list[index] = True
            else:
                try:
                    if isinstance(float(value.replace(" ", "")), float):
                        value_check_list[index] = True

                except ValueError:
                    continue
        return value_check_list

    @staticmethod
    def check_number_atoms(value_list, atom_type, smiles):
        """
            Check that the number of atoms in the list is between 0 and the amount of atoms in the
            molecule structure provided by the smiles string.

            :param value_list: a list of the values and ranges in the peak list.
            :param atom_type: the type of atom being checked either H or C.
            :param smiles: the smiles string for the corresponding compound. 
            :return: If input is valid, return True. If input is invalid, raise InvalidAtomNumer exception.
        """
        if not value_list:
            return True
            
        atoms = {"H": 1, "C": 6}
        atom_num = atoms[atom_type]
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        atom_query = rdqueries.AtomNumEqualsQueryAtom(atom_num)
        number_of_atoms_1 = len(mol.GetAtomsMatchingQuery(atom_query))

        # number_of_atoms_2 = 0
        # for atom in mol.GetAtoms():
        #     if atom.GetSymbol() == atom_type:
        #         number_of_atoms_2 += 1

        # patt = Chem.MolFromSmarts(f"[{atom_type}]")
        # number_of_atoms_3 = len(mol.GetSubstructMatches(patt))

        # sanity check that all 3 methods of finding the number of atoms in the molecule are returning the same amount
        # assert number_of_atoms_1 == number_of_atoms_2 == number_of_atoms_3, "UnEqual number of atoms"\

        if 0 < len(value_list) <= number_of_atoms_1:
            return True
        else:
            raise InvalidAtomNumber

    @staticmethod
    def check_value_ranges_C_H(value_list, atom_type):
        """
            Ensure the numerical range for each value makes sense for the type of atom. Ex. H > -2 etc..

            :param value_list: a list of the values and ranges in the peak list.
            :param atom_type: the type of atom being checked either H or C
            :return: If valid input, return True. If invalid input, return error message.
        """
        unchecked_values = []
        for value in value_list:
            is_range = Validate.is_valid_range(value)
            if is_range:
                split_value = re.split("[\u002D\u05BE\u1806\u2010\u2011\u2012\u2013\u2014\u2015\u207B\u208B\u2212\uFE58\uFE63\uFF0D]",value)
                split_value = [unchecked_values.append(float(re.sub("[ ()]", "", sub_str))) for sub_str in split_value]
            else:
                unchecked_values.append(float(value))

        error_values = []
        warning_values = []
        has_error = False
        has_warning = False
        for value in unchecked_values:
            if atom_type == "H":
                try:
                    if value < -2 or value > 20:
                        has_error = True
                        error_values.append(value)
                    elif value < -1 or value > 16:
                        has_warning = True
                        warning_values.append(value)
                except ErrorBadRange as exc:
                    raise exc
            if atom_type == "C":
                try:
                    if value < -10 or value >250:
                        has_error = True
                        error_values.append(value)
                    elif value < 0 or value >230:
                        has_warning = True
                        warning_values.append(value)
                except ErrorBadRange as exc:
                    raise exc
        
        if has_error:
            raise ErrorBadRange(bad_value=error_values)
        if has_warning:
            raise WarnBadRange(bad_value=warning_values)

        return "Valid"

    @staticmethod
    def check_value_ranges_other(value, feature):
        """
            Ensure the numerical value for the specified feature (temperature, frequency)

            :param value: feature value.
            :param feature: the feature to validate ("temperature" or "frequency")
            :return: If valid input, return True. If invalid input, return error message.
        """
        if feature == "temperature":
            if value:
                try:
                    if value < 200 or value > 400:
                        raise ErrorBadRange(bad_value=value)
                    elif value < 250 or value > 350:
                        raise WarnBadRange(bad_value=value)
                except ErrorBadRange as exc:
                    raise exc

        elif feature == "frequency":
            try:
                if value < 50 or value > 1400:
                    raise ErrorBadRange(bad_value=value)
                elif value < 75 or value > 900:
                    raise WarnBadRange(bad_value=value)
            except ErrorBadRange as exc:
                raise exc
        
        return "Valid"
        

    @staticmethod
    def validate(
        smiles,
        solvent,
        reference,
        H_text_block=None,
        C_text_block=None,
        h_frequency=None,
        h_temperature=None,
        c_frequency=None,
        c_temperature=None,
    ):
        """
            Check that the peak lists given check some basic validity checks before accepting them into the DB.

            :param H_text_block: The text entered into the H peak list text box.
            :param C_text_block: The text entered into the C peak list text box.
            :param smiles: The smiles string for the corresponding compound.
            :return: If the input is valid, return a string confirming its validity. If the input is not valid, 
                    raise an Exception specifying the problem with the input.
        """
        # Tuple to append warning messages to (all in first string since that is easier to display)
        warning_message = ["", "Warning"]

        # Check that characters in the text block are valid
        if not solvent:
            return ("No solvent provided", "Error")
        if not reference:
            return ("No reference provided", "Error")
        if not H_text_block and not C_text_block:
            return ("At least one of H Values or C Values is required", "Error")
        if H_text_block and (not h_frequency):
            return ("No hydrogen frequency provided", "Error")
        if C_text_block and (not c_frequency):
            return ("No carbon frequency provided", "Error")
        if (h_frequency or h_temperature) and (not H_text_block):
            return ("No hydrogen values provided. Either add a list of values or remove the hydrogen frequency and temperature values.", "Error")
        if (c_frequency or c_temperature) and (not C_text_block):
            return ("No carbon values provided. Either add a list of values or remove the carbon frequency and temperature values.", "Error")

        if H_text_block:
            try:
                Validate.check_valid_characters(H_text_block)
            except Exception as exc:
                if exc.error == "InvalidCharacters":
                    return ("Invalid Characters in H List: Please make sure that only contains the following allowed characters 0-9 , . - ; ()", "Error")
            # Parse the text blocks into lists based on the separators
            try:
                H_list = Validate.parse_text_to_list(H_text_block)
            except NoSplit:
                return ("Failed to split H list, please check your separators.", "Error")
        else:
            H_list = None
                
        if C_text_block:
            try:
                Validate.check_valid_characters(C_text_block)
            except Exception as exc:
                if exc.error == "InvalidCharacters":
                    return ("Invalid Characters in C List: Please make sure that only contains the following allowed characters 0-9 , . - ; ()", "Error")
            try:
                C_list = Validate.parse_text_to_list(C_text_block)
            except NoSplit:
                return ("Failed to split C list, please check your separators.", "Error")
        else:
            C_list = None
        
        
        # Check if each element in the parsed lists are either floats or ranges
        if H_list:
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
                return (f"Error: The following values in the H list are of invalid data type {error_value_list}", "Error")
            
            # Check the number of atoms in each list do not exceed amount of atoms in struct
            try:
                Validate.check_number_atoms(H_list, "H", smiles)
            except InvalidAtomNumber:
                return ("Error: Invalid number of H atoms in the peak list", "Error")
            
            if not h_temperature:
                warning_message[0] += "No hydrogen temperature provided.\n"
                
            try:    
                Validate.check_value_ranges_C_H(H_list, "H")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Hydrogen peak value(s) {exc.bad_value} out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Hydrogen peak value(s) {exc.bad_value} outside of the typical H value range.\n"
            
            try:
                Validate.check_value_ranges_other(h_temperature, "temperature")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Hydrogen Temperature value {exc.bad_value} K is out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Hydrogen Temperature {exc.bad_value} K is outside of the typical temperature value range.\n"
                    
            try:
                Validate.check_value_ranges_other(h_frequency, "frequency")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Hydrogen Frequency value {exc.bad_value} MHz is out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Hydrogen Frequency {exc.bad_value} MHz is outside of a the typical frequency value range.\n"
        
        if C_list:
            # Check the datatype for each of the entries in each list
            try:
                error_list = Validate.check_data_type(C_list)
                if sum(error_list) != len(C_list):
                    error_value_list = []
                    for index, b_val in enumerate(H_list):
                        if b_val is False:
                            error_value_list.append(H_list[index])
                        else:
                            continue
                    raise InvalidValueType
            except InvalidValueType:
                return (f"Error: The following values in the C list are of invalid data type {error_value_list}", "Error")
            
            try:
                Validate.check_number_atoms(C_list, "C", smiles)
            except InvalidAtomNumber:
                return ("Error: Invalid number of C atoms in the peak list", "Error")

            if not c_temperature:
                warning_message[0] += "No carbon temperature provided.\n"

            try:
                Validate.check_value_ranges_C_H(C_list, "C")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Carbon peak value(s) {exc.bad_value} out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Carbon peak value(s) {exc.bad_value} outside of the typical C value range.\n"
        
            try:
                Validate.check_value_ranges_other(c_temperature, "temperature")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Carbon Temperature value {exc.bad_value} K is out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Carbon Temperature {exc.bad_value} K is outside of the typical temperature value range.\n"
        
            try:
                Validate.check_value_ranges_other(c_frequency, "frequency")
            except (ErrorBadRange, WarnBadRange) as exc:
                if exc.error_type == "error":
                    return (f"Carbon Frequency value {exc.bad_value} MHz is out of the accepted range", "Error")
                elif exc.error_type == "warning":
                    warning_message[0] += f"Carbon Frequency {exc.bad_value} MHz is outside of a the typical frequency value range.\n"
        
        if H_list and not C_list:
            warning_message[0] += f"No carbon list provided. Only hydrogen values will be submitted.\n"
        if C_list and not H_list:
            warning_message[0] += f"No hydrongen list provided. Only carbon values will be submitted.\n"

        # empty_message = ["", "Empty"]
        if not H_list and not C_list:
            return (f'Empty: Both Lists contain no peaks. At least one list is required (you may also submit no peak lists for this compound by clicking "remove").', "Error")
        # if not c_frequency or c_frequency=="\n":
        #     empty_message[0] += "Empty C list: If you do not wish to submit a peak list for this compound please click the remove button before clicking the submit button"
        # if not h_frequency or h_frequency=="\n":
        #     empty_message[0] += "Empty H list: If you do not wish to submit a peak list for this compound please click the remove button before clicking the submit button"

        # if empty_message[0]:
        #     empty_message[0] = empty_message[0].rsplit('\n', 1)[0]
        #     return tuple(empty_message)

        if warning_message[0]:
            warning_message[0] = warning_message[0].rsplit('\n', 1)[0]
            return tuple(warning_message)

        return ("Both lists are valid", "No Errors")


class Convert:

    @staticmethod
    def convert_to_float_list(peak_string):
        peak_list = Validate.parse_text_to_list(peak_string)
        output_list = []
        for value in peak_list:
            if Validate.is_valid_range(value) is True:
                split_values = re.split("[\u002D\u05BE\u1806\u2010\u2011\u2012\u2013\u2014\u2015\u207B\u208B\u2212\uFE58\uFE63\uFF0D]",value)
                split_values = [float(re.sub("[ ()]", "", sub_str)) for sub_str in split_values]
                range_tuple = (split_values[0], split_values[1])
                output_list.append(range_tuple)
            else:
                value = float(value)
                output_list.append(value)
        return output_list

    @staticmethod
    def sort_list_desc(peak_list):
        map_dict = {}
        list_to_sort = []
        for value in peak_list:
            if type(value) is tuple:
                key_value = value[0]
            else:
                key_value = value
            if key_value not in map_dict:
                map_dict[key_value] = []
            map_dict[key_value].append(value)
            if key_value not in list_to_sort:
                list_to_sort.append(key_value)

        list_to_sort.sort(reverse=True)

        output_list = []
        for number in list_to_sort:
            output_list.extend(map_dict[number])

        return output_list

    @staticmethod
    def convert(peak_string):
        peak_list = Convert.convert_to_float_list(peak_string=peak_string)
        sorted_peak_list = Convert.sort_list_desc(peak_list=peak_list)
        return sorted_peak_list
        

