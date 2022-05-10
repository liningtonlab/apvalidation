import statistics as stats
from rdkit import Chem

"""
Create three classes, one to validate Proton peaks and one to validate Carbon peaks
"""

class Validate:
    """
    This will be the class that contains the validate_peak_lists function.
    This function will be the one called by the user in order to validate the peaks.
    """
    @staticmethod
    def check_lengths(peak_dict):
        h_peaks = peak_dict['H']
        c_peaks = peak_dict['C']
        if len(h_peaks) != len(c_peaks):
            return False
        return True

    @staticmethod
    def validate_peak_lists(peak_dict, smiles):

        return_bool = True
        return_messages = []

        if Validate.check_lengths(peak_dict) is False:
            return_bool = False
            return_messages.append("Length of carbon and proton list do not match.")

        # currently returns the number of unique protons doesn't work properly I need to know what to check
        num_protons = Proton.check_atom_counts(smiles)
        print(num_protons)
        if num_protons is False:
            return_bool = False
            return_messages.append("Invalid number of unique protons.")

        # currently returns the number of carbons doesn't work properly I need to know what to check
        num_carbons = Carbon.check_atom_counts(smiles)
        print(num_carbons)
        if num_carbons is False:
            return_bool = False
            return_messages.append("Invalid number of unique carbons.")  

        h_list = peak_dict['H']
        if Proton.check_value_ranges(h_list) is False:
            return_messages.append("WARNING: Proton values out of range.")

        c_list = peak_dict['C']
        if Carbon.check_values_ranges(c_list) is False:
            return_messages.append("WARNING: Carbon values out of range.")

        return return_bool, return_messages

class Proton:

    @staticmethod
    def check_value_ranges(h_list):
        mean = stats.mean(h_list)
        median = stats.median(h_list)
        if mean > 10:
            return False
        if median > 10:
            return False
        
        for item in h_list:
            if item > 20:
                return False
        return True

    def check_atom_counts(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol_h = Chem.AddHs(mol)
        ordered_atoms = Chem.CanonicalRankAtoms(mol_h, breakTies=False)
        unique_hs = set()
        all_hs = []
        for atom, sym_classes in zip(mol_h.GetAtoms(), ordered_atoms):
            if atom.GetAtomicNum() == 1:
                all_hs.append(atom)
                unique_hs.add(sym_classes)

        return len(unique_hs)


class Carbon:

    @staticmethod
    def check_values_ranges(c_list):
        mean = stats.mean(c_list)
        median = stats.median(c_list)
        if mean < 10:
            return False
        if median < 10:
            return False
        
        for item in c_list:
            if item > 250:
                return False
        return True

    @staticmethod
    def check_atom_counts(smiles):
        mol = Chem.MolFromSmiles(smiles)
        all_cs = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                all_cs.append(atom)

        return len(all_cs)



    