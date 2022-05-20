import re
import statistics as stats
import warnings

from rdkit import Chem
from rdkit.Chem import rdqueries, rdmolops



"""
These are custom exception classes to help annotate what error we are handling.
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
        """
            Check that only accepted characters are present in the text taken from the input text box. 

            :param text_block: The text taken from the user input
            :return: If valid, return a string confirming valid chars. If invalid raise InvalidCharacters exception.
        """
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
        value_list = [value.replace(" ","") for value in value_list]

        for index, value in enumerate(value_list):
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
        
        return "Valid"
        

    @staticmethod
    def validate(H_text_block, C_text_block, smiles):
        """
            Check that the peak lists given check some basic validity checks before accepting them into the DB.

            :param H_text_block: The text entered into the H peak list text box.
            :param C_text_block: The text entered into the C peak list text box.
            :param smiles: The smiles string for the corresponding compound.
            :return: If the input is valid, return a string confirming its validity. If the input is not valid, 
                    raise an Exception specifying the problem with the input.
        """

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
            return f"Error: The following values are of invalid data type {error_value_list}"

        # Check the number of atoms in each list do not exceed amount of atoms in struct
        try:
            Validate.check_number_atoms(H_list, "H", smiles)
        except InvalidAtomNumber:
            return "Error: Invalid number of H atoms in the peak list"
        try:
            Validate.check_number_atoms(C_list, "C", smiles)
        except InvalidAtomNumber:
            return "Error: Invalid number of C atoms in the peak list"

        # Check the values to ensure they are real H or C values
        try:
            output = Validate.check_value_ranges(H_list, "H")
        except ErrorBadRange:
            return output
        try:
            output = Validate.check_value_ranges(C_list, "C")
        except ErrorBadRange:
            return output

        return "Both lists are valid"

