from itertools import count
from collections import defaultdict
import re
import numpy as np
import pandas as pd
from tree import molecule_breaker
from configparser import ConfigParser


glyco_mass_native = {
    'N': 203.079372533,
    'H': 162.0528234315,
    'F': 146.0579088094,
    'A': 291.09541652769997,
    'G': 307.09033114979997,
    'X': 132.0422587452,
    'pH': 242.0191544315,
    'aH': 179.0793635,
    'J': 203.079372533 + 18.01056,
    'Q': 162.0528234315 + 18.01056
}
glyco_mass_reduced = {
    'N': 203.079372533,
    'H': 162.0528234315,
    'F': 146.0579088094,
    'A': 291.09541652769997,
    'G': 307.09033114979997,
    'X': 132.0422587452,
    'pH': 242.0191544315,
    'aH': 179.0793635,
    'J': 203.079372533 + 20.02621,
    'Q': 162.0528234315 + 20.02621
}
ATOM_MAPPING = {
    "H": {"C": 6, "H": 10, "O": 5, "N": 0},
    "Q": {"C": 6, "H": 10, "O": 5, "N": 0},
    "N": {"C": 8, "H": 13, "O": 5, "N": 1},
    "J": {"C": 8, "H": 13, "O": 5, "N": 1},
    "F": {"C": 6, "H": 10, "O": 4, "N": 0},
    "A": {"C": 11, "H": 17, "O": 8, "N": 1},
    "G": {"C": 11, "H": 17, "O": 9, "N": 1},
    "X": {"C": 5, "H": 8, "O": 4, "N": 0},
}
molecule_mapping = {
    "H": {"C": 6, "H": 10, "O": 5, "N": 0},
    "Q": {"C": 6, "H": 10, "O": 5, "N": 0},
    "N": {"C": 8, "H": 13, "O": 5, "N": 1},
    "J": {"C": 8, "H": 13, "O": 5, "N": 1},
    "F": {"C": 6, "H": 10, "O": 4, "N": 0},
    "A": {"C": 11, "H": 17, "O": 8, "N": 1},
    "G": {"C": 11, "H": 17, "O": 9, "N": 1},
    "X": {"C": 5, "H": 8, "O": 4, "N": 0},
}

def calc_glycan_mass(glycan, glyco_condition=glyco_mass_native):
    glycos = glycan.strip(')').split(')')
    mass = 0
    for glyco in glycos:
        glyco, n = glyco.split('(')
        mass += glyco_condition[glyco] * int(n)
    return mass

def parse_nested_structure(structure):
    """
    Parse the  molecular formula structure and return the total count of each atom.
    """
    stack = []
    current_atoms = defaultdict(int)
    i = 0

    while i < len(structure):
        char = structure[i]

        if char in ATOM_MAPPING:
            for atom, count in ATOM_MAPPING[char].items():
                current_atoms[atom] += count

        elif char == '(':
            stack.append(current_atoms)
            current_atoms = defaultdict(int)

        elif char == ')':
            nested_atoms = current_atoms
            current_atoms = stack.pop()
            for atom, count in nested_atoms.items():
                current_atoms[atom] += count

        i += 1

    return current_atoms

def format_molecule_formula(atom_count):
    """
    Format the atom counts into a molecular formula string.
    """
    formula = ""
    for atom in ["C", "H", "O", "N"]:
        if atom_count[atom] > 0:
            formula += f"{atom}{atom_count[atom]}"
    return formula

def calculate_molecule_formula(molecule_name):
    """
    Calculate the molecular formula based on the molecular structure.
    """
    atom_count = parse_nested_structure(molecule_name)
    return format_molecule_formula(atom_count)

def calculate_accurate_mass(formula):
    """
    accurate mass
    """
    isotope_masses = {
        "C": 12.000000,
        "H": 1.007825,
        "O": 15.994915,
        "N": 14.003074
    }

    pattern = r"([A-Z][a-z]?)(\d*)"
    matches = re.findall(pattern, formula)
    accurate_mass = 0.0
    for element, count in matches:
        count = int(count) if count else 1
        if element in isotope_masses:
            accurate_mass += isotope_masses[element] * count
        else:
            raise ValueError(f"Unknown element: {element}")
    return accurate_mass

def parse_glycan_composition(formula):
    """
    Parse the molecular structure and count the number of specific glycosyl components.
    """
    element_mapping = {
        "H": "H",
        "Q": "H",  # Q ---> H
        "N": "N",
        "J": "N",  # J ---> N
        "F": "F",
        "A": "A",
        "G": "G",
        "X": "X"
    }
    target_elements = ["H", "N", "F", "A", "G", "X"]

    def parse_structure(structure):
        stack = []
        current_count = defaultdict(int)
        i = 0
        while i < len(structure):
            char = structure[i]
            if char in element_mapping:
                mapped_char = element_mapping[char]
                current_count[mapped_char] += 1
            elif char == '(':
                stack.append(current_count)
                current_count = defaultdict(int)
            elif char == ')':
                nested_count = current_count
                current_count = stack.pop()
                for element, count in nested_count.items():
                    current_count[element] += count
            i += 1
        return current_count

    total_count = parse_structure(formula)
    result = ""
    for element in target_elements:
        if total_count[element] > 0:
            result += f"{element}({total_count[element]})"
    return result

def calculate_product_formula(composition):
    """
     J(1)N(1)H(1)--> formula。
    """
    total_atoms = defaultdict(int)
    pattern = r"([A-Z])\((\d+)\)"
    matches = re.findall(pattern, composition)

    for symbol, count in matches:
        count = int(count)
        if symbol not in molecule_mapping:
            raise ValueError(f"unknown symbol: {symbol}")
        formula = molecule_mapping[symbol]
        for atom, atom_count in formula.items():
            total_atoms[atom] += atom_count * count

    final_formula = "".join(f"{atom}{total_atoms[atom]}" for atom in sorted(total_atoms.keys()))
    return final_formula

def molecule_product_results(df):
    """
    process DataFrame，calculate molecular formula, mass
    """
    df['Molecule Formula'] = df['Molecule Structure'].apply(calculate_molecule_formula)

    df['Molecule Mass (native)'] = df['Molecule Formula'].apply(calculate_accurate_mass) + 18.01056
    df['Molecule Mass (reduced)'] = df['Molecule Formula'].apply(calculate_accurate_mass) + 20.02621

    df['Product Charge'] = 1

    df['Molecule Name'] = df['Molecule Structure'].apply(parse_glycan_composition)

    df['Product Formula'] = df['Product Name'].apply(calculate_product_formula)

    symbols = ['H', 'N', 'F', 'A', 'G', 'X']

    def parse_molecule_name(name):
        counts = {symbol: 0 for symbol in symbols}
        matches = re.findall(r"([A-Z])\((\d+)\)", name)
        for symbol, count in matches:
            if symbol in counts:
                counts[symbol] = int(count)
        if counts['G'] == 0 and counts['X'] == 0:
            return "_".join(str(counts[symbol]) for symbol in ['H', 'N', 'F', 'A'])
        return "_".join(str(counts[symbol]) for symbol in symbols)

    #  Formula Base
    df['Formula Base'] = df['Molecule Name'].apply(parse_molecule_name)

    # Assign suffixes based on the molecular structure.
    def assign_suffix(group):
        if group['Molecule Structure'].nunique() == 1:
            group['Processed Formula'] = group['Formula Base'] + "_a"
        else:
            suffix_map = {structure: chr(97 + i) for i, structure in enumerate(group['Molecule Structure'].unique())}
            group['Processed Formula'] = group['Formula Base'] + "_" + group['Molecule Structure'].map(suffix_map)
        return group

    df = df.groupby('Formula Base', group_keys=False).apply(assign_suffix)
    df.drop(columns=['Formula Base'], inplace=True)

    df["Molecule Name"] = df['Processed Formula']

    # reorder
    output_columns = [
        "Molecule Name",
        "Molecule Structure",
        "Molecule Formula",
        "Molecule Mass (native)",
        "Molecule Mass (reduced)",
        "Product Name",
        "Product Formula",
        "Product m/z (native)",
        "Product m/z (reduced)",
        "Product Charge"
    ]
    df = df[output_columns]

    return df

def generate_charge_ms(df, windows_start=600, windows_end=1200):
    """
    generate charge and m/z under specified windows
    """
    rows = []
    for _, row in df.iterrows():
        charge_num = 3
        while charge_num:
            charge = charge_num
            mass_native = row['Molecule Mass (native)']
            mass_reduced = row['Molecule Mass (reduced)']
            mz_native = (mass_native + charge * 1.00783) / charge
            mz_reduced = (mass_reduced + charge * 1.00783) / charge
            if windows_start <= mz_native <= windows_end:
                new_row = row.copy()
                new_row['Precursor charge'] = charge
                new_row['Precursor m/z (native)'] = mz_native
                new_row['Precursor m/z (reduced)'] = mz_reduced
                rows.append(new_row)
            charge_num -= 1

    if rows:
        df_new = pd.DataFrame(rows)
        new_columns = [
            "Molecule Name",
            "Molecule Structure",
            "Molecule Formula",
            "Precursor charge",
            "Precursor m/z (native)",
            "Precursor m/z (reduced)",
            "Product Name",
            "Product Formula",
            "Product m/z (native)",
            "Product m/z (reduced)",
            "Product Charge"
        ]
        df_new = df_new[new_columns]
        return df_new
    else:
        return pd.DataFrame(columns=[
            "Molecule Name",
            "Molecule Structure",
            "Molecule Formula",
            "Precursor charge",
            "Precursor m/z (native)",
            "Precursor m/z (reduced)",
            "Product Name",
            "Product Formula",
            "Product m/z (native)",
            "Product m/z (reduced)",
            "Product Charge"
        ])

def calc_glycan_ions_mass(df, count=2):
    """
    calculate fragments mass
    """
    records = []
    for _, row in df.iterrows():
        molecule_structure = row['Molecule Structure']
        ions = molecule_breaker(molecule_structure, count)
        ion_masses_native = [calc_glycan_mass(glycan, glyco_mass_native) + 1.00783 for glycan in ions]
        ion_masses_reduced = [calc_glycan_mass(glycan, glyco_mass_reduced) + 1.00783 for glycan in ions]
        for glycan, mass_n, mass_r in zip(ions, ion_masses_native, ion_masses_reduced):
            records.append({
                "Molecule Structure": molecule_structure,
                "Product Name": glycan,
                "Product m/z (native)": mass_n,
                "Product m/z (reduced)": mass_r
            })
    return pd.DataFrame(records)

def process_glyco_data(input_file, output_file, count=3, windows_start=600, windows_end=1200):

    # step 1：read input file
    initial_df = pd.read_csv(input_file, header=None, names=["Molecule Structure"], sep="\t")

    # step 2：calculate ions mass
    ions_df = calc_glycan_ions_mass(initial_df, count=count)

    # step 3：process precursor
    processed_df = molecule_product_results(ions_df)

    # step 4：Generate charges and m/z values within the specified window.
    final_df = generate_charge_ms(processed_df, windows_start=windows_start, windows_end=windows_end)

    # step 5：save results
    final_df.to_csv(output_file, index=False, sep="\t")
    print("done, exit")


if __name__ == "__main__":
    #process_glyco_data("glyco_library/N-Human.txt", "glyco_library/N-Human-glycan-ions.txt",count=2,windows_start=600, windows_end=1200)
    #process_glyco_data("glyco_library/N-Mouse.txt", "glyco_library/N-Mouse-glycan-ions.txt",count=2,windows_start=600, windows_end=1200)
    #process_glyco_data("glyco_library/N-Plant.txt", "glyco_library/N-Plant-glycan-ions.txt",count=2,windows_start=600, windows_end=1200)
    #process_glyco_data("glyco_library/O-Glycan.txt", "glyco_library/O-Glycan-ions.txt",count=2,windows_start=600, windows_end=1200)
   # process_glyco_data("glyco_library/test.txt", "glyco_library/output_final.txt",count=3,windows_start=600, windows_end=1200)

    cfg = ConfigParser()
    cfg.read("./glycan_ions_config.ini")
    cfg_dict = dict(cfg.items("config"))
    print("Configuration loaded:", cfg_dict)

    input_path = cfg_dict["input_path"]
    output_path = cfg_dict["output_path"]
    count = int(cfg_dict["count"])
    windows_start = int(cfg_dict["windows_start"])
    windows_end = int(cfg_dict["windows_end"])

    process_glyco_data(
        input_file=input_path,
        output_file=output_path,
        count=count,
        windows_start=windows_start,
        windows_end=windows_end
    )