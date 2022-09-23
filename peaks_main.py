from ctypes import c_long
import apvalidation.peak_validator as peaks

unicode_dashes = " - ־ ᠆ ‐ ‑ ‒ - – — ― ⁻ ₋ − ﹘ ﹣ －"

def test_validate():
    smiles_string = "O=C1CCC[C@@]2([H])C(CCC[C@@]21[H])=O"

    H_valid = "13.0, (13.5-13.9), 15.5, 15.3, 12.4 ₋ 14.6"
    H_bug = "10, (10-11), 10"
    H_invalid_characters = "1H NMR (CD3OD, 500 MHz) δ 1.03 (d, J = 7.0 Hz, 3H, H-8)"
    H_invalid_list_format = "13 (13.5-14.2) 34.5 67.3 12.4 ₋ 14.6"
    H_warning_range = "40, (13.5-13.9), 15.5, 15.3, 12.4 ₋ 14.6"
    H_invalid_range = "13, (13.5 14.2), 34.5, 67.3, (12.4, 14.6)"
    H_invalid_number = "10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10"
    temp_warning = 360
    freq_warning = 1500
    temp_error = -20
    freq_error = -20

    C_valid = "20, (20.5-25.2), 64.7, 76.3, 21.4 ₋ 22.6"
    C_invalid_characters = "1H NMR (CD3OD, 500 MHz) δ 1.03 (d, J = 7.0 Hz, 3H, H-8)"
    C_invalid_list_format = "11 (11.5-15.2) 64.7 76.3 1.4 ₋ 19.6"
    C_warning_range = "15, (20.5-25.2), 64.7, 76.3, 21.4 ₋ 22.6"
    C_invalid_range = "11, (11.5,15.2), 64.7, 76.3, 1.4 ₋ 19.6"
    
    temp_valid = 300
    freq_valid = 800

    output_valid = peaks.Validate.validate(
        H_text_block = H_valid,
        C_text_block = C_valid, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_valid,
        h_temperature = temp_valid,
        c_frequency = freq_valid,
        c_temperature = temp_valid,
        reference = "test_reference"
    )
    print(f"Valid: {output_valid}")
    invalid_number = peaks.Validate.validate(
        H_text_block = H_invalid_number,
        C_text_block = C_valid, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_valid,
        h_temperature = temp_valid,
        c_frequency = freq_valid,
        c_temperature = temp_valid,
        reference = "test_reference"
    )
    print(f"Invalid Number of H atoms: {invalid_number}")
    output_warning_ranges = peaks.Validate.validate(
        H_text_block = H_warning_range,
        C_text_block = C_warning_range, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_warning,
        h_temperature = temp_warning,
        c_frequency = freq_warning,
        c_temperature = temp_warning,
        reference = "test_reference",
    )
    print(f"Warning Ranges: {output_warning_ranges}")
    invalid_temp_freq = peaks.Validate.validate(
        H_text_block = H_valid,
        C_text_block = C_valid, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_valid,
        h_temperature = temp_error,
        c_frequency = freq_valid,
        c_temperature = temp_error,
        reference = "test_reference"
    )
    print(f"Invalid Temp Freq: {invalid_temp_freq}")
    output_invalid_characters = peaks.Validate.validate(
        H_text_block = H_invalid_range,
        C_text_block = C_invalid_range, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_warning,
        h_temperature = temp_warning,
        c_frequency = freq_warning,
        c_temperature = temp_warning,
        reference = "test_reference",
    )
    print(f"Invalid Characters: {output_invalid_characters}")
    output_invalid_list_format = peaks.Validate.validate(
        H_text_block = H_invalid_list_format,
        C_text_block = C_invalid_list_format, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = temp_valid,
        h_temperature = freq_valid,
        c_frequency = temp_valid,
        c_temperature = freq_valid,
        reference = "test_reference",
    )
    print(f"Invalid List: {output_invalid_list_format}")
    output_invalid_range = peaks.Validate.validate(
        H_text_block = H_invalid_range,
        C_text_block = C_invalid_range, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = temp_valid,
        h_temperature = freq_valid,
        c_frequency = temp_valid,
        c_temperature = freq_valid,
        reference = "test_reference",
    )
    print(f"Invalid Range: {output_invalid_range}")
    output_no_temp = peaks.Validate.validate(
        H_text_block = H_valid,
        C_text_block = C_valid, 
        smiles = smiles_string,
        solvent = "test_solvent",
        h_frequency = freq_valid,
        h_temperature = None,
        c_frequency = freq_valid,
        c_temperature = None,
        reference = "test_reference"
    )
    print(f"No Temperature: {output_no_temp}")

    test_parse = peaks.Convert.convert_to_float_list(H_valid)
    test_parse_sort = peaks.Convert.sort_list_desc(test_parse)
    print(f"TEST CONVERT VALID: {test_parse_sort}")

    test_parse_bug = peaks.Convert.convert_to_float_list(H_bug)
    print(f"TEST CONVERT BUG: {test_parse_bug}")
    test_parse_bug_sort = peaks.Convert.sort_list_desc(test_parse_bug)
    print(f"TEST CONVERT BUG ORDER: {test_parse_bug_sort}")


    legacy_ver_test = peaks.Validate.legacy_validate(
        H_text_block = H_valid,
        C_text_block = C_valid, 
        smiles = smiles_string,
    )
    print(f"Legacy validate: {legacy_ver_test}")

# def test_sanghoon():

#     H_test1 = "11.87, 8.85, 8.83, 8.83, 8.52, 8.50, 8.14, 8.10, 8.09, 8.08, 8.08, 8.07, 8.07, 7.62, 7.62, 7.61, 7.61, 7.61, 7.61, 7.60, 7.60, 4.18, 2.46"
#     C_test1 = "166.83, 155.98, 154.49, 152.69, 149.21, 146.77, 137.35, 124.64, 121.07, 120.78, 102.91, 56.34, 17.82"
#     smiles_test1 = "COC1=CC(C2=CC=CC=N2)=NC(/C=N/O)=C1SC"
#     output_test1 = peaks.Validate.validate(H_test1, C_test1, smiles_test1, )
#     print(f"Test1: {output_test1}")
#     if output_test1 == "Both lists are valid":
#         converted_H_list1 = peaks.Convert.convert(H_test1)
#         converted_C_list1 = peaks.Convert.convert(C_test1)
#         print(f"Converted H test1: {converted_H_list1}")
#         print(f"Converted C test1: {converted_C_list1}")


# def test_convert():
#     H_invalid = "13.0, (13.5-13.9), 34.5, 67.3, 12.4 ₋ 14.6"
#     C_valid = "140, 25, 75, 100, 54"
#     smiles_test = "COC1=CC(C2=CC=CC=N2)=NC(/C=N/O)=C1SC"
#     validation_test = peaks.Validate.validate(H_invalid, C_valid, smiles_test)
#     print(validation_test)
#     if validation_test == "Both lists are valid":

#         converted_list = peaks.Convert.convert(H_invalid)
#         original_list = H_invalid


#         print(f'Converted: {converted_list}')
#         print(f'Original {original_list}')


# Test Validate
test_validate()
# test_convert()
# test_sanghoon()





