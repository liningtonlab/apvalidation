import apvalidation.peak_validator as peaks

peak_dict = {"H": [9.4, 7.3, 16.3, 15], "C": [75.3, 153.8, 45.8, 200.4]}

smiles_string = "CC(=O)NCCC1=CNc2c1cc(OC)cc2"
num_unique_h = peaks.Validate.validate_peak_lists(peak_dict=peak_dict, smiles=smiles_string)
print(f"The number of unique hydrogens is: {num_unique_h}")
