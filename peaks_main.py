import apvalidation.peak_validator as peaks

unicode_dashes = " - ־ ᠆ ‐ ‑ ‒ - – — ― ⁻ ₋ − ﹘ ﹣ －"

def test_validate():
    smiles_string = "O=C1CCC[C@@]2([H])C(CCC[C@@]21[H])=O"

    H_valid = "13.0, (13.5-13.9), 34.5, 67.3, 12.4 ₋ 14.6"
    H_invalid_characters = "1H NMR (CD3OD, 500 MHz) δ 1.03 (d, J = 7.0 Hz, 3H, H-8)"
    H_invalid_list_format = "13 (13.5-14.2) 34.5 67.3 12.4 ₋ 14.6"
    H_invalid_range = "13, (13.5 14.2), 34.5, 67.3, (12.4, 14.6)"

    C_valid = "11, (11.5-15.2), 64.7, 76.3, 1.4 ₋ 19.6"
    C_invalid_characters = "1H NMR (CD3OD, 500 MHz) δ 1.03 (d, J = 7.0 Hz, 3H, H-8)"
    C_invalid_list_format = "11 (11.5-15.2) 64.7 76.3 1.4 ₋ 19.6"
    C_invalid_range = "11, (11.5,15.2), 64.7, 76.3, 1.4 ₋ 19.6"

    output_valid = peaks.Validate.validate(H_valid, C_valid, smiles_string)
    print(f"Valid: {output_valid}")
    output_invalid_characters = peaks.Validate.validate(H_invalid_characters, C_invalid_characters, smiles_string)
    print(f"Invalid Characters: {output_invalid_characters}")
    output_invalid_list_format = peaks.Validate.validate(H_invalid_list_format, C_invalid_list_format, smiles_string)
    print(f"Invalid List: {output_invalid_list_format}")
    output_invalid_range = peaks.Validate.validate(H_invalid_range, C_invalid_range, smiles_string)
    print(f"Invalid Range: {output_invalid_range}")


# Test Validate
test_validate()