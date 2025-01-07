import glob
import os
from fileinput import filename
from numpy import dot
from numpy.linalg import norm
import pandas as pd
import re
from tqdm import tqdm
import matplotlib.pyplot as plt
from collections import defaultdict
from pythoms.molecule import IPMolecule
from configparser import ConfigParser

def extract_mgf_data(input_folder):
    """
    Recursively search a specified folder and its subfolders for all .mgf files,
    extract the TITLE, PEPMASS, CHARGE, and MS1_IDX information, and compile the results into a Pandas DataFrame.

    Parameters:
        input_folder (str): The path to the main folder to be searched.

    Return:
        pandas.DataFrame: dataframe ['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS']。
    """

    mgf_files = glob.glob(os.path.join(input_folder, '**', '*.mgf'), recursive=True)

    if not mgf_files:
        print(f"No .mgf files found in '{input_folder}' and its subdirectories.")
        return pd.DataFrame(columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    data_list = []

    title_pattern = re.compile(r'^TITLE=(.*)', re.IGNORECASE)
    pep_mass_pattern = re.compile(r'^PEPMASS=(\S+)', re.IGNORECASE)
    charge_pattern = re.compile(r'^CHARGE=(\S+)', re.IGNORECASE)
    ms1_idx_pattern = re.compile(r'^MS1_IDX=(\d+)', re.IGNORECASE)

    # Iterate  each *.mgf file.


    for file_path in mgf_files:
        try:
            with open(file_path, 'r') as file:
                content = file.read()
        except Exception as e:
            print(f"Error reading file '{file_path}': {e}")
            continue

        if not content.strip():
            print(f"Warning: '{file_path}' is empty. Skipping.")
            continue

        # split content by  BEGIN IONS and END IONS
        blocks = content.split("BEGIN IONS")
        for block in blocks[1:]:
            block = block.strip()
            if not block:
                continue
            # 分割到 END IONS
            block_content = block.split("END IONS")[0]
            record = {}

            lines = block_content.splitlines()
            for line in lines:
                line = line.strip()
                # TITLE
                title_match = title_pattern.match(line)
                if title_match:
                    record['TITLE'] = title_match.group(1)
                    continue

                # PEPMASS
                pep_mass_match = pep_mass_pattern.match(line)
                if pep_mass_match:
                    record['PEPMASS'] = pep_mass_match.group(1)
                    continue

                #CHARGE
                charge_match = charge_pattern.match(line)
                if charge_match:
                    record['CHARGE'] = charge_match.group(1)
                    continue

                # MS1_IDX
                ms1_idx_match = ms1_idx_pattern.match(line)
                if ms1_idx_match:
                    record['MS1_IDX'] = ms1_idx_match.group(1)
                    continue

            if 'TITLE' in record:
                record.setdefault('PEPMASS', None)
                record.setdefault('CHARGE', None)
                record.setdefault('MS1_IDX', None)
                data_list.append(record)
            else:
                print(f"Warning: 'TITLE' not found in a block of '{file_path}'. Skipping this block.")
                continue

    if not data_list:
        print("No valid data extracted from the .mgf files.")
        return pd.DataFrame(columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    df = pd.DataFrame(data_list, columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    return df
def extract_from_debug(mgf_df, debug_file_path):
    """
    integrate information from mgf_df and debug_file

    Parameters:
        mgf_df (pd.DataFrame): contain column of  'TITLE' and  'MS1_IDX''
        debug_file_path (str): debug file path。

    Returns:
        pd.DataFrame: columns：
            ['TITLE', 'MS1_IDX', ‘charge’，'spectrum_scan_num',
             'isotopic1', 'isotopic2', 'isotopic3',
             'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    """
    # Step 1:  debug_file
    data_entries = []
    current_cpd_addon = None
    current_spectrum = {}

    cpd_addon_pattern = re.compile(r'^cpd_addon\s+(.+)$', re.IGNORECASE)
    spectrum_ms1_idx_pattern = re.compile(r'^specturm MS1 idx:\s*(\d+)', re.IGNORECASE)
    spectrum_scan_num_pattern = re.compile(r'^spectrum scan num:\s*(\d+)', re.IGNORECASE)
    df_pattern = re.compile(r'^df(\d+)\s+\[(\d+\.\d+)\]\s+(\d+\.\d+)', re.IGNORECASE)

    try:
        with open(debug_file_path, 'r') as file:
            for line in tqdm(file, desc="Parsing debug_file"):
                line = line.strip()

                # Check if it is a cpd_addon row.
                cpd_match = cpd_addon_pattern.match(line)
                if cpd_match:
                    current_cpd_addon = cpd_match.group(1).strip()
                    continue

                # Check if it is a ----processing Spectrum---- row.
                if line.startswith('----processing Spectrum----'):
                    current_spectrum = {}
                    continue

                # spectrum MS1 idx
                ms1_match = spectrum_ms1_idx_pattern.match(line)
                if ms1_match and current_cpd_addon is not None:
                    current_spectrum['MS1_IDX'] = ms1_match.group(1)
                    continue

                #  spectrum scan num
                scan_num_match = spectrum_scan_num_pattern.match(line)
                if scan_num_match and 'MS1_IDX' in current_spectrum:
                    current_spectrum['spectrum_scan_num'] = scan_num_match.group(1)
                    continue

                # df1, df2, df3
                df_match = df_pattern.match(line)
                if df_match and 'MS1_IDX' in current_spectrum:
                    df_num = df_match.group(1)  # df1, df2, df3  m/z
                    isotopic = df_match.group(2)  # isotopic value
                    peak = df_match.group(3)  # peak intensity

                    # assign isotopic and peak
                    if df_num == '1':
                        current_spectrum['isotopic1'] = isotopic
                        current_spectrum['isotopic1_peak'] = peak
                    elif df_num == '2':
                        current_spectrum['isotopic2'] = isotopic
                        current_spectrum['isotopic2_peak'] = peak
                    elif df_num == '3':
                        current_spectrum['isotopic3'] = isotopic
                        current_spectrum['isotopic3_peak'] = peak

                        if current_cpd_addon is not None:
                            record = {
                                'cpd_addon': current_cpd_addon,
                                'MS1_IDX': current_spectrum.get('MS1_IDX'),
                                'spectrum_scan_num': current_spectrum.get('spectrum_scan_num'),
                                'isotopic1': current_spectrum.get('isotopic1'),
                                'isotopic2': current_spectrum.get('isotopic2'),
                                'isotopic3': current_spectrum.get('isotopic3'),
                                'isotopic1_peak': current_spectrum.get('isotopic1_peak'),
                                'isotopic2_peak': current_spectrum.get('isotopic2_peak'),
                                'isotopic3_peak': current_spectrum.get('isotopic3_peak'),
                            }
                            data_entries.append(record)
                            current_spectrum = {}
                    continue

        debug_df = pd.DataFrame(data_entries)
        if debug_df.empty:
            print("Warning: No valid data was extracted from the debug_file.")
    except FileNotFoundError:
        print(f"Warning:  '{debug_file_path}' was not found。")
        return pd.DataFrame()
    except Exception as e:
        print(f"An error occurred while parsing the debug file: {e}")
        return pd.DataFrame()

    if debug_df.empty:
        return pd.DataFrame(columns=[
            'TITLE', 'MS1_IDX', 'spectrum_scan_num',
            'isotopic1', 'isotopic2', 'isotopic3',
            'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak'
        ])

    # Step 2: Create a lookup key.

    mgf_df = mgf_df.copy()
    #  'cpd_addon'
    mgf_df['cpd_addon'] = mgf_df['TITLE'].str.extract(r'^(.+?)(?:\(\d+\))?$', expand=False)
    mgf_df['cpd_addon'] = mgf_df['cpd_addon'].fillna('')
    mgf_df['cpd_addon'] = mgf_df['cpd_addon'].str.strip()

    mgf_df['MS1_IDX'] = mgf_df['MS1_IDX'].astype(str)
    debug_df['MS1_IDX'] = debug_df['MS1_IDX'].astype(str)

    mgf_df['merge_key'] = mgf_df['cpd_addon'] + '_' + mgf_df['MS1_IDX']
    debug_df['merge_key'] = debug_df['cpd_addon'] + '_' + debug_df['MS1_IDX']

    # Step 3: merge mgf_df and debug_df
    merged_df = pd.merge(mgf_df, debug_df, on='merge_key', how='left', suffixes=('', '_debug'))

    # Step 4: select and  rename colums and deduplicates
    summary_df = merged_df.groupby(['TITLE', 'MS1_IDX']).agg({
        'spectrum_scan_num': 'first',
        'CHARGE': 'first',
        'isotopic1': 'first',
        'isotopic2': 'first',
        'isotopic3': 'first',
        'isotopic1_peak': 'first',
        'isotopic2_peak': 'first',
        'isotopic3_peak': 'first'
    }).reset_index()
    return summary_df
def normalize_peaks(df, peak_cols, reference_col, new_prefix='norm'):
    """
        Normalize the specified list of peaks based on the reference column.
    Parameters:
        df (pd.DataFrame): input dataframe
        peak_cols (list of str): List of column names to be normalized.
        reference_col (str): The column name to be used as the normalization reference.
        new_prefix (str): The prefix for the new normalized columns.

    Returns:
        pd.DataFrame: The normalized DataFrame with additional columns.
    """
    df[reference_col] = pd.to_numeric(df[reference_col], errors='coerce')

    for col in peak_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        new_col = f"{new_prefix}_{col}"
        df[new_col] = (df[col] / df[reference_col])*100

    return df
isotope_data = {
    'C': [(12.0000, 0.9893), (13.0034, 0.0107)],
    'H': [(1.0078, 0.999885), (2.0141, 0.000115)],
    'O': [(15.9949, 0.99757), (16.9991, 0.00038), (17.9992, 0.00205)],
    'N': [(14.0031, 0.99632), (15.0001, 0.00368)],
}
def calculate_isotope_peaks_theoretical(formula, charge=1, threshold=0, mass_precision=0.01):
    """
    Calculate the m/z values and intensities of the isotope peaks for a chemical formula, considering the charge state.

    Parameters:
        formula (str):  "C6H12O6"。
        charge (int): charge number。
        threshold (float): Intensity threshold, peaks with intensities below this value will be ignored (default is 0.001).
        mass_precision (float): Mass precision is used to round mass values to a specific number of decimal places (default is 0.01).

    Return:
        peaks (list of tuples): peak (m/z, intensity)，sort by intensity。
    """
    if charge == 0:
        raise ValueError("charge num can't be 0.")

    #  "C6H12O6" -> {'C': 6, 'H': 12, 'O': 6}
    parsed_formula = {}
    pattern = r'([A-Z][a-z]*)(\d*)'
    for (element, count) in re.findall(pattern, formula):
        count = int(count) if count else 1
        parsed_formula[element] = parsed_formula.get(element, 0) + count

    overall = {0.0: 1.0}
    for element, count in parsed_formula.items():
        if element not in isotope_data:
            raise ValueError(f"不支持的元素: {element}")

        isotopes = isotope_data[element]

        element_dist = defaultdict(float)
        element_dist[0.0] = 1.0

        for _ in range(count):
            temp_dist = defaultdict(float)
            for mass1, intensity1 in element_dist.items():
                for mass2, intensity2 in isotopes:
                    new_mass = round(mass1 + mass2, 2)
                    new_intensity = intensity1 * intensity2
                    temp_dist[new_mass] += new_intensity
            element_dist = temp_dist

        temp_overall = defaultdict(float)
        for mass1, intensity1 in overall.items():
            for mass2, intensity2 in element_dist.items():
                new_mass = round(mass1 + mass2, 2)
                new_intensity = intensity1 * intensity2
                temp_overall[new_mass] += new_intensity
        overall = temp_overall

    peaks = [(mass+18.01056, intensity) for mass, intensity in overall.items() if intensity >= threshold]
    peaks = sorted(peaks, key=lambda x: x[1], reverse=True)

    adjusted_peaks = []
    for mass, intensity in peaks:
        mz = (mass + charge * 1.007276466812) / abs(charge)
        mz = round(mz, 4)
        adjusted_peaks.append((mz, intensity))
    return adjusted_peaks
def generate_isotope_peaks_theoretical_pythoms(formula: str, charge: int, ipmethod: str = 'isospec', resolution: int = 20000):
    """
    Generate a list of isotope pattern peak intensities for a given molecular formula and charge state.

    Parameters:
        formula (str): such as  "C62H1045"。
        charge (int): charge number。
        ipmethod (str): Isotope Pattern Calculation Method， 'isospec'（ 'combinatorics' or 'multiplicative'）。
        resolution (int): MS resolution, 20000。

    返回:
        peaks (list of tuples): peak (m/z, intensity)。
    """
    #  IPMolecule class
    mol = IPMolecule(
        formula,
        charge=charge,
        ipmethod=ipmethod,
        resolution=resolution,
        verbose=True
    )

    # obtain isotopic pattern
    mz_values, intensities = mol.raw_isotope_pattern

    # Merge adjacent m/z values with a difference less than 0.1.
    merged_mz = []
    merged_intensities = []

    temp_mz = mz_values[0]
    temp_intensity = intensities[0]

    for i in range(1, len(mz_values)):
        if abs(mz_values[i] - temp_mz) < 0.1:
            # Merge  m/z values and intensity
            temp_intensity += intensities[i]
        else:
            merged_mz.append(temp_mz)
            merged_intensities.append(temp_intensity)
            temp_mz = mz_values[i]
            temp_intensity = intensities[i]

    merged_mz.append(temp_mz)
    merged_intensities.append(temp_intensity)

    peaks = list(zip(merged_mz, merged_intensities))
    return peaks
def parse_molecular_formula(strings):
    # composition of the glycans
    monomers = {
        'H': {'C': 6, 'H': 10, 'O': 5, 'N': 0},
        'N': {'C': 8, 'H': 13, 'O': 5, 'N': 1},
        'F': {'C': 6, 'H': 10, 'O': 4, 'N': 0},
        'A': {'C': 11, 'H': 17, 'O': 8, 'N': 1}
    }
    results = []
    for s in strings:
        #  number of H, N, F, A
        counts = list(map(int, s.split('-')[0].split('_')[:4]))
        if len(counts) != 4:
            raise ValueError(f"error: {s}")

        total_elements = {'C': 0, 'H': 0, 'O': 0, 'N': 0}

        # total number of elements
        for monomer, count in zip(monomers.keys(), counts):
            for element, number in monomers[monomer].items():
                total_elements[element] += number * count
        # total number + 2H + O
        total_elements['H'] += 2
        total_elements['O'] += 1

        formula = ''.join(f"{el}{total_elements[el]}" for el in sorted(total_elements) if total_elements[el] > 0)
        results.append(formula)

    return results
def plot_isotopic_peaks(data):
    """
    Plot the barplot of experimental and theoretical isotope peaks.

    Parameters：
        data (pd.Series): contain isotopic peak information from plot_dat.iloc[0]。
    """
    plt.rcParams.update({'font.size': 14})
    #exp_mz = [data['isotopic1'], data['isotopic2'], data['isotopic3']]
    exp_intensity = [data['norm_isotopic_isotopic1_peak'], data['norm_isotopic_isotopic2_peak'],
                     data['norm_isotopic_isotopic3_peak']]


    theo_intensity = [data['T_isotopic1_peaks'], data['T_isotopic2_peaks'],
                      data['T_isotopic3_peaks']]
    exp_mz = [float(data['isotopic1']), float(data['isotopic2']), float(data['isotopic3'])]
    theo_mz = [float(data['T_isotopic1']), float(data['T_isotopic2']), float(data['T_isotopic3'])]

    fig, ax = plt.subplots(figsize=(6, 6))
    # experimental plot
    ax.bar(exp_mz, exp_intensity, width=0.01, label='from experiment', alpha=0.8, color='black')
    # theoretical plot
    ax.bar(theo_mz, theo_intensity, width=0.01, label='from theoretical', alpha=0.8, color='orange')

    label = f"{data['TITLE']}({data['cosine_similarity']})"
    ax.set_title(label)
    ax.set_xticks(exp_mz)
    ax.set_xticklabels([f'{mz:.4f}' for mz in exp_mz], rotation=45, ha='right')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Relative Intensity')
    ax.set_ylim(1, 100)

    ax.legend()
    ax.grid(False)

    # show the figure
    plt.tight_layout()
    filename = f"{data['TITLE']}_isotopic_peak_distribution.svg"
    folder = "distribution_plot"
    file_path = os.path.join(folder, filename)
    plt.savefig(file_path, format='svg', dpi=300)
    plt.show()
def calculate_cosine_similarity(row):
    # vector1
    vector1 = [
        row['norm_isotopic_isotopic1_peak'],
        row['norm_isotopic_isotopic2_peak'],
        row['norm_isotopic_isotopic3_peak']
    ]
    # vector2
    vector2 = [
        row['T_isotopic1_peaks'],
        row['T_isotopic1_peaks'],
        row['T_isotopic1_peaks']
    ]
    cosine_similarity = dot(vector1, vector2) / (norm(vector1) * norm(vector2))
    cosine_similarity = round(cosine_similarity, 3)
    return cosine_similarity
def isotopic_peak_distribution(input_folder,debug_file):
    # 1. mgf_df
    mgf_df = extract_mgf_data(input_folder)
    # 2. debug_df
    debug_df = extract_from_debug(mgf_df, debug_file)

    formulas = parse_molecular_formula(debug_df['TITLE'])
    debug_df['formula'] = formulas
    debug_df['T_isotopic1'] = None
    debug_df['T_isotopic1_peaks'] = None
    debug_df['T_isotopic2']= None
    debug_df['T_isotopic2_peaks']= None
    debug_df['T_isotopic3']= None
    debug_df['T_isotopic3_peaks'] =None
    charges = debug_df['CHARGE']
    charge = []
    for s in charges:
        match = re.search(r'(\d+)\+', s)
        if match:
            number = int(match.group(1))
            charge.append(number)
        else:
            charge.append(None)
    peaks_list = [generate_isotope_peaks_theoretical_pythoms(formula, charge=charge[index])[:3] for index, formula in
                  enumerate(formulas)]
    print(peaks_list)
    for i, peaks in enumerate(peaks_list):
        debug_df ['T_isotopic1'][i] = peaks[0][0] + 1.007276466812
        debug_df['T_isotopic1_peaks'][i] = peaks[0][1]
        debug_df['T_isotopic2'][i] = peaks[1][0] + 1.007276466812
        debug_df['T_isotopic2_peaks'][i] = peaks[1][1]
        debug_df['T_isotopic3'][i] = peaks[2][0] + 1.007276466812
        debug_df['T_isotopic3_peaks'][i] = peaks[2][1]

    # 3. Extract peaks for plotting and normalize the peaks.
    columns_to_extract = [
        'TITLE', 'formula', 'isotopic1', 'isotopic2', 'isotopic3',
        'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak',
        'T_isotopic1', 'T_isotopic2', 'T_isotopic3',
        'T_isotopic1_peaks', 'T_isotopic2_peaks', 'T_isotopic3_peaks'
    ]

    plot_dat = debug_df.loc[:, columns_to_extract]
    isotope_cols = [ 'isotopic1', 'isotopic2', 'isotopic3']
    T_isotopic_cols = ['T_isotopic1', 'T_isotopic2', 'T_isotopic3']
    isotopic_peak_cols = ['isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    T_isotopic_peak_cols = ['T_isotopic1_peaks', 'T_isotopic2_peaks', 'T_isotopic3_peaks']

    # covert to the numeric format
    for col in isotopic_peak_cols + T_isotopic_peak_cols + isotope_cols + T_isotopic_cols:
        plot_dat[col] = pd.to_numeric(plot_dat[col], errors='coerce')

    # standardize isotopic_peak
    isotopic_peak_cols = ['isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    plot_dat = normalize_peaks(plot_dat, isotopic_peak_cols, 'isotopic1_peak', new_prefix='norm_isotopic')

    # standardize T_isotopic_peak
    plot_dat['cosine_similarity'] = plot_dat.apply(calculate_cosine_similarity, axis=1)
    plot_dat.to_csv('distribution_plot/isotopic_peak_from_exp_solico.csv', index=False)
    print(plot_dat.iloc[0,])
    plot_dat.apply(lambda row: plot_isotopic_peaks(row), axis=1)

if __name__ == "__main__":
    cfg = ConfigParser()
    cfg.read("./isotopic_peak_distribution_config.ini")
    cfg_dict = dict(cfg.items("config"))
    print("Configuration loaded:", cfg_dict)
    input_folder = cfg_dict["input_folder"]
    debug_file = cfg_dict["debug_file"]
    isotopic_peak_distribution(input_folder = input_folder,debug_file =debug_file)
