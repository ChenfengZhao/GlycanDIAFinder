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

def extract_mgf_data(input_folder):
    """
    递归搜索指定文件夹及其子文件夹中的所有 .mgf 文件，提取 TITLE、PEPMASS、CHARGE 和 MS1_IDX 信息，
    并将结果汇总到一个 Pandas DataFrame 中。

    参数:
        input_folder (str): 要搜索的主文件夹的路径。

    返回:
        pandas.DataFrame: 包含所有提取信息的 DataFrame，列名为 ['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS']。
    """

    # 使用 glob 递归查找所有 .mgf 文件
    mgf_files = glob.glob(os.path.join(input_folder, '**', '*.mgf'), recursive=True)

    # 如果没有找到任何 .mgf 文件，提示用户并返回空的 DataFrame
    if not mgf_files:
        print(f"No .mgf files found in '{input_folder}' and its subdirectories.")
        return pd.DataFrame(columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    # 初始化一个空的列表，用于存储所有提取的数据
    data_list = []

    # 定义正则表达式模式
    title_pattern = re.compile(r'^TITLE=(.*)', re.IGNORECASE)
    pep_mass_pattern = re.compile(r'^PEPMASS=(\S+)', re.IGNORECASE)
    charge_pattern = re.compile(r'^CHARGE=(\S+)', re.IGNORECASE)
    ms1_idx_pattern = re.compile(r'^MS1_IDX=(\d+)', re.IGNORECASE)

    # 遍历每个 .mgf 文件
    for file_path in mgf_files:
        try:
            with open(file_path, 'r') as file:
                content = file.read()
        except Exception as e:
            print(f"Error reading file '{file_path}': {e}")
            continue  # 跳过无法读取的文件

        # 判断文件是否为空（去除空白字符后）
        if not content.strip():
            print(f"Warning: '{file_path}' is empty. Skipping.")
            continue  # 跳过空文件

        # 按照 BEGIN IONS 和 END IONS 分割文件内容
        blocks = content.split("BEGIN IONS")
        for block in blocks[1:]:  # 第一个分割结果可能是空的
            block = block.strip()
            if not block:
                continue
            # 分割到 END IONS
            block_content = block.split("END IONS")[0]

            # 初始化一个字典来存储当前块的信息
            record = {}

            # 按行遍历当前块的内容
            lines = block_content.splitlines()
            for line in lines:
                line = line.strip()

                # 搜索 TITLE
                title_match = title_pattern.match(line)
                if title_match:
                    record['TITLE'] = title_match.group(1)
                    continue

                # 搜索 PEPMASS
                pep_mass_match = pep_mass_pattern.match(line)
                if pep_mass_match:
                    record['PEPMASS'] = pep_mass_match.group(1)
                    continue

                # 搜索 CHARGE
                charge_match = charge_pattern.match(line)
                if charge_match:
                    record['CHARGE'] = charge_match.group(1)
                    continue

                # 搜索 MS1_IDX
                ms1_idx_match = ms1_idx_pattern.match(line)
                if ms1_idx_match:
                    record['MS1_IDX'] = ms1_idx_match.group(1)
                    continue

            # 如果找到了 TITLE，则添加记录
            if 'TITLE' in record:
                # 如果某些字段未找到，可以设置为 None
                record.setdefault('PEPMASS', None)
                record.setdefault('CHARGE', None)
                record.setdefault('MS1_IDX', None)

                data_list.append(record)
            else:
                print(f"Warning: 'TITLE' not found in a block of '{file_path}'. Skipping this block.")
                continue

    # 检查是否有任何数据被提取
    if not data_list:
        print("No valid data extracted from the .mgf files.")
        return pd.DataFrame(columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    # 将收集到的数据转换为 Pandas DataFrame
    df = pd.DataFrame(data_list, columns=['TITLE', 'MS1_IDX', 'CHARGE', 'PEPMASS'])

    return df
def extract_from_debug(mgf_df, debug_file_path):
    """
    汇总 mgf_df 和 debug_file 中的信息，返回一个包含所需字段的 DataFrame。

    参数:
        mgf_df (pd.DataFrame): 包含至少 'TITLE' 和 'MS1_IDX' 两列的 DataFrame。
        debug_file_path (str): 调试文件的路径。

    返回:
        pd.DataFrame: 汇总后的 DataFrame，包含以下列：
            ['TITLE', 'MS1_IDX', ‘charge’，'spectrum_scan_num',
             'isotopic1', 'isotopic2', 'isotopic3',
             'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    """
    # Step 1: 解析 debug_file
    data_entries = []  # 用于存储所有提取的光谱信息
    current_cpd_addon = None  # 当前的 cpd_addon 标识符
    current_spectrum = {}  # 当前光谱的信息

    # 定义正则表达式模式
    cpd_addon_pattern = re.compile(r'^cpd_addon\s+(.+)$', re.IGNORECASE)

    spectrum_ms1_idx_pattern = re.compile(r'^specturm MS1 idx:\s*(\d+)', re.IGNORECASE)
    spectrum_scan_num_pattern = re.compile(r'^spectrum scan num:\s*(\d+)', re.IGNORECASE)
    df_pattern = re.compile(r'^df(\d+)\s+\[(\d+\.\d+)\]\s+(\d+\.\d+)', re.IGNORECASE)

    try:
        with open(debug_file_path, 'r') as file:
            for line in tqdm(file, desc="Parsing debug_file"):
                line = line.strip()

                # 检查是否为 cpd_addon 行
                cpd_match = cpd_addon_pattern.match(line)
                if cpd_match:
                    current_cpd_addon = cpd_match.group(1).strip()
                    continue  # 继续下一行

                # 检查是否为 '----processing Spectrum----' 行
                if line.startswith('----processing Spectrum----'):
                    current_spectrum = {}  # 重置当前光谱信息
                    continue  # 继续下一行

                # 检查是否为 spectrum MS1 idx 行
                ms1_match = spectrum_ms1_idx_pattern.match(line)
                if ms1_match and current_cpd_addon is not None:
                    current_spectrum['MS1_IDX'] = ms1_match.group(1)
                    continue  # 继续下一行

                # 检查是否为 spectrum scan num 行
                scan_num_match = spectrum_scan_num_pattern.match(line)
                if scan_num_match and 'MS1_IDX' in current_spectrum:
                    current_spectrum['spectrum_scan_num'] = scan_num_match.group(1)
                    continue  # 继续下一行

                # 检查是否为 df1, df2, df3 行
                df_match = df_pattern.match(line)
                if df_match and 'MS1_IDX' in current_spectrum:
                    df_num = df_match.group(1)  # df1, df2, df3 的数字部分
                    isotopic = df_match.group(2)  # isotopic 值
                    peak = df_match.group(3)  # 峰值

                    # 根据 df_num 分配 isotopic 和 peak
                    if df_num == '1':
                        current_spectrum['isotopic1'] = isotopic
                        current_spectrum['isotopic1_peak'] = peak
                    elif df_num == '2':
                        current_spectrum['isotopic2'] = isotopic
                        current_spectrum['isotopic2_peak'] = peak
                    elif df_num == '3':
                        current_spectrum['isotopic3'] = isotopic
                        current_spectrum['isotopic3_peak'] = peak
                        # 假设 df3 是一个光谱条目的结束，保存当前条目
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
                            current_spectrum = {}  # 重置当前光谱信息
                    continue  # 继续下一行

        # 将数据条目转换为 DataFrame
        debug_df = pd.DataFrame(data_entries)
        if debug_df.empty:
            print("警告: 从 debug_file 中未提取到任何有效数据。")
    except FileNotFoundError:
        print(f"错误: 调试文件 '{debug_file_path}' 未找到。")
        return pd.DataFrame()
    except Exception as e:
        print(f"解析调试文件时发生错误: {e}")
        return pd.DataFrame()

    # 如果 debug_df 为空，返回空的 DataFrame
    if debug_df.empty:
        return pd.DataFrame(columns=[
            'TITLE', 'MS1_IDX', 'spectrum_scan_num',
            'isotopic1', 'isotopic2', 'isotopic3',
            'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak'
        ])

    # Step 2: 创建查找键
    # 从 'TITLE' 中提取标识符，忽略括号及其中的内容
    mgf_df = mgf_df.copy()
    # 提取 'cpd_addon'，忽略括号及其中的内容
    mgf_df['cpd_addon'] = mgf_df['TITLE'].str.extract(r'^(.+?)(?:\(\d+\))?$', expand=False)
    # 处理可能的空值
    mgf_df['cpd_addon'] = mgf_df['cpd_addon'].fillna('')
    mgf_df['cpd_addon'] = mgf_df['cpd_addon'].str.strip()

    # 将 'MS1_IDX' 转换为字符串类型
    mgf_df['MS1_IDX'] = mgf_df['MS1_IDX'].astype(str)
    debug_df['MS1_IDX'] = debug_df['MS1_IDX'].astype(str)

    # 创建合并键
    mgf_df['merge_key'] = mgf_df['cpd_addon'] + '_' + mgf_df['MS1_IDX']
    debug_df['merge_key'] = debug_df['cpd_addon'] + '_' + debug_df['MS1_IDX']

    # Step 3: 合并 mgf_df 与 debug_df
    merged_df = pd.merge(mgf_df, debug_df, on='merge_key', how='left', suffixes=('', '_debug'))

    # Step 4: 选择并重命名所需的列
    # 如果存在重复的 'merge_key'，可以根据需要进行聚合

    summary_df = merged_df.groupby(['TITLE', 'MS1_IDX']).agg({
        'spectrum_scan_num': 'first',  # 选择第一个扫描号
        'CHARGE': 'first',
        'isotopic1': 'first',
        'isotopic2': 'first',
        'isotopic3': 'first',
        'isotopic1_peak': 'first',
        'isotopic2_peak': 'first',
        'isotopic3_peak': 'first'
    }).reset_index()
   # summary_df = merged_df.drop(columns=['cpd_addon_x', 'cpd_addon_y', 'merge_key'])
    return summary_df
# 手动定义每种元素的同位素质量和丰度
def normalize_peaks(df, peak_cols, reference_col, new_prefix='norm'):
    """
    标准化指定的峰列，以参考列为基准。

    参数:
        df (pd.DataFrame): 输入的数据框。
        peak_cols (list of str): 需要标准化的列名列表。
        reference_col (str): 作为标准化参考的列名。
        new_prefix (str): 新标准化列的前缀。

    返回:
        pd.DataFrame: 标准化后的数据框，新增列。
    """
    # 确保 reference_col 为数值类型
    df[reference_col] = pd.to_numeric(df[reference_col], errors='coerce')

    for col in peak_cols:
        # 确保需要标准化的列为数值类型
        df[col] = pd.to_numeric(df[col], errors='coerce')

        # 创建新列
        new_col = f"{new_prefix}_{col}"
        # 避免除以零或 NaN
        df[new_col] = df[col] / df[reference_col]

    return df
isotope_data = {
    'C': [(12.0000, 0.9893), (13.0034, 0.0107)],
    'H': [(1.0078, 0.999885), (2.0141, 0.000115)],
    'O': [(15.9949, 0.99757), (16.9991, 0.00038), (17.9992, 0.00205)],
    'N': [(14.0031, 0.99632), (15.0001, 0.00368)],
    # 添加其他元素的数据
}

def calculate_isotope_peaks_theoretical(formula, charge=1, threshold=0, mass_precision=0.01):
    """
    计算化学式的同位素峰的 m/z 和强度，考虑电荷数。

    参数:
        formula (str): 化学式，例如 "C6H12O6"。
        charge (int): 电荷数，默认为 1。
        threshold (float): 强度阈值，低于该值的峰将被忽略 (默认值为 0.001)。
        mass_precision (float): 质量精度，用于将质量值固定到某一小数位 (默认值为 0.01)。

    返回:
        peaks (list of tuples): 每个峰的 (m/z, intensity)，按照强度降序排列。
    """
    if charge == 0:
        raise ValueError("电荷数不能为零。")

    # 解析化学式，例如 "C6H12O6" -> {'C': 6, 'H': 12, 'O': 6}
    parsed_formula = {}
    pattern = r'([A-Z][a-z]*)(\d*)'
    for (element, count) in re.findall(pattern, formula):
        count = int(count) if count else 1
        parsed_formula[element] = parsed_formula.get(element, 0) + count

    # 初始化整体同位素分布
    overall = {0.0: 1.0}

    # 对每种元素，计算其同位素分布并与整体分布进行卷积
    for element, count in parsed_formula.items():
        if element not in isotope_data:
            raise ValueError(f"不支持的元素: {element}")

        # 获取该元素的同位素数据
        isotopes = isotope_data[element]

        # 计算该元素的同位素分布
        element_dist = defaultdict(float)
        element_dist[0.0] = 1.0

        for _ in range(count):
            temp_dist = defaultdict(float)
            for mass1, intensity1 in element_dist.items():
                for mass2, intensity2 in isotopes:
                    new_mass = round(mass1 + mass2, 2)  # 固定精度
                    new_intensity = intensity1 * intensity2
                    temp_dist[new_mass] += new_intensity
            element_dist = temp_dist

        # 将当前元素的同位素分布与整体分布进行卷积
        temp_overall = defaultdict(float)
        for mass1, intensity1 in overall.items():
            for mass2, intensity2 in element_dist.items():
                new_mass = round(mass1 + mass2, 2)  # 固定精度
                new_intensity = intensity1 * intensity2
                temp_overall[new_mass] += new_intensity
        overall = temp_overall

    # 过滤并收集最终的峰
    peaks = [(mass+18.01056, intensity) for mass, intensity in overall.items() if intensity >= threshold]
    peaks = sorted(peaks, key=lambda x: x[1], reverse=True)

    # 考虑电荷数，计算 m/z
    adjusted_peaks = []
    for mass, intensity in peaks:
        mz = (mass + charge * 1.007276466812) / abs(charge)  # 1.007276466812 为质子的质量
        mz = round(mz, 4)  # 固定 m/z 的精度
        adjusted_peaks.append((mz, intensity))
    return adjusted_peaks


def parse_molecular_formula(strings):
    # 定义单体的元素组成
    monomers = {
        'H': {'C': 6, 'H': 10, 'O': 5, 'N': 0},
        'N': {'C': 8, 'H': 13, 'O': 5, 'N': 1},
        'F': {'C': 6, 'H': 10, 'O': 4, 'N': 0},
        'A': {'C': 11, 'H': 17, 'O': 8, 'N': 1}
    }

    # 存储结果的列表
    results = []

    # 处理每个字符串
    for s in strings:
        # 提取前四个数字，表示 H, N, F, A 的数量
        counts = list(map(int, s.split('-')[0].split('_')[:4]))
        if len(counts) != 4:
            raise ValueError(f"字符串格式错误: {s}")

        # 初始化元素计数
        total_elements = {'C': 0, 'H': 0, 'O': 0, 'N': 0}

        # 计算总的元素数量
        for monomer, count in zip(monomers.keys(), counts):
            for element, number in monomers[monomer].items():
                total_elements[element] += number * count

        # 生成化学分子式
        formula = ''.join(f"{el}{total_elements[el]}" for el in sorted(total_elements) if total_elements[el] > 0)
        results.append(formula)

    return results

def plot_isotopic_peaks(data):
    """
    绘制实验和理论同位素峰的条形图。

    参数：
        data (pd.Series): 包含同位素峰信息的 Series，来自 plot_dat.iloc[0]。
    """
    plt.rcParams.update({'font.size': 14})
    #exp_mz = [data['isotopic1'], data['isotopic2'], data['isotopic3']]
    exp_intensity = [data['norm_isotopic_isotopic1_peak'], data['norm_isotopic_isotopic2_peak'],
                     data['norm_isotopic_isotopic3_peak']]

    # 提取理论同位素数据
    #theo_mz = [data['T_isotopic1'], data['T_isotopic2'], data['T_isotopic3']]
    theo_intensity = [data['norm_T_isotopic_T_isotopic1_peaks'], data['norm_T_isotopic_T_isotopic2_peaks'],
                      data['norm_T_isotopic_T_isotopic3_peaks']]
    exp_mz = [float(data['isotopic1']), float(data['isotopic2']), float(data['isotopic3'])]
    theo_mz = [float(data['T_isotopic1']), float(data['T_isotopic2']), float(data['T_isotopic3'])]

    # 创建图形和轴
    fig, ax = plt.subplots(figsize=(6, 6))
    # 绘制实验同位素峰的条形图
    ax.bar(exp_mz, exp_intensity, width=0.01, label='from experiment', alpha=1, color='black')
    # 绘制理论同位素峰的条形图
    ax.bar(theo_mz, theo_intensity, width=0.01, label='from theoretical', alpha=1, color='orange')


    # 设置标题和标签
    label = f"{data['TITLE']}({data['cosine_similarity']})"
    ax.set_title(label)
    ax.set_xticks(exp_mz)
    ax.set_xticklabels([f'{mz:.4f}' for mz in exp_mz], rotation=45, ha='right')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_ylim(0, 1)

    # 设置 x 轴刻度
    #all_mz = sorted(set(exp_mz + theo_mz))
    #ax.set_xticks(all_mz)
    #ax.set_xticklabels([f'{mz:.4f}' for mz in all_mz], rotation=45, ha='right')

    # 显示图例
    ax.legend()

    # 显示网格
    ax.grid(False)

    # 显示图形
    plt.tight_layout()
    filename = f"{data['TITLE']}_isotopic_peak_distribution.svg"
    folder = "distribution_plot"
    file_path = os.path.join(folder, filename)
    plt.savefig(file_path, format='svg', dpi=300)
    plt.show()

# 计算余弦相似性
def calculate_cosine_similarity(row):
    # 提取向量1
    vector1 = [
        row['norm_isotopic_isotopic1_peak'],
        row['norm_isotopic_isotopic2_peak'],
        row['norm_isotopic_isotopic3_peak']
    ]
    # 提取向量2
    vector2 = [
        row['norm_T_isotopic_T_isotopic1_peaks'],
        row['norm_T_isotopic_T_isotopic2_peaks'],
        row['norm_T_isotopic_T_isotopic3_peaks']
    ]
    cosine_similarity = dot(vector1, vector2) / (norm(vector1) * norm(vector2))
    cosine_similarity = round(cosine_similarity, 3)
    return cosine_similarity
def isotopic_peak_distribution(input_folder,debug_file):
    # 1. 获取 mgf_df
    mgf_df = extract_mgf_data(input_folder)

    # 2. 获取 debug_df
    debug_df = extract_from_debug(mgf_df, debug_file)
    #print(mgf_df)
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
    # 遍历列表中的每个元素
    for s in charges:
        match = re.search(r'(\d+)\+', s)  # 匹配 "+" 之前的数字
        if match:
            number = int(match.group(1))  # 提取并转换为整数
            charge.append(number)  # 将结果添加到列表中
        else:
            charge.append(None)  # 如果没有匹配，存储 None
    peaks_list = [calculate_isotope_peaks_theoretical(formula, charge=charge[index])[:3] for index,formula in enumerate(formulas)]

    for i, peaks in enumerate(peaks_list):
        debug_df ['T_isotopic1'][i] = peaks[0][0]
        debug_df['T_isotopic1_peaks'][i] = peaks[0][1]
        debug_df['T_isotopic2'][i] = peaks[1][0]
        debug_df['T_isotopic2_peaks'][i] = peaks[1][1]
        debug_df['T_isotopic3'][i] = peaks[2][0]
        debug_df['T_isotopic3_peaks'][i] = peaks[2][1]
    #print(debug_df.iloc[0])
    # 3. 提取peaks绘图并对peaks进行标准化

    columns_to_extract = [
        'TITLE', 'formula', 'isotopic1', 'isotopic2', 'isotopic3',
        'isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak',
        'T_isotopic1', 'T_isotopic2', 'T_isotopic3',
        'T_isotopic1_peaks', 'T_isotopic2_peaks', 'T_isotopic3_peaks'
    ]

    # 使用 loc 提取指定的列
    plot_dat = debug_df.loc[:, columns_to_extract]
    # 定义需要转换的列
    isotope_cols = [ 'isotopic1', 'isotopic2', 'isotopic3']
    T_isotopic_cols = ['T_isotopic1', 'T_isotopic2', 'T_isotopic3']
    isotopic_peak_cols = ['isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    T_isotopic_peak_cols = ['T_isotopic1_peaks', 'T_isotopic2_peaks', 'T_isotopic3_peaks']

    # 转换为数值类型，使用 pd.to_numeric 并设置 errors='coerce' 将无法转换的值设为 NaN
    for col in isotopic_peak_cols + T_isotopic_peak_cols + isotope_cols + T_isotopic_cols:
        plot_dat[col] = pd.to_numeric(plot_dat[col], errors='coerce')
    #print(plot_dat.iloc[0,])
    # 标准化 isotopic_peak 列
    isotopic_peak_cols = ['isotopic1_peak', 'isotopic2_peak', 'isotopic3_peak']
    plot_dat = normalize_peaks(plot_dat, isotopic_peak_cols, 'isotopic1_peak', new_prefix='norm_isotopic')

    # 标准化 T_isotopic_peak 列
    T_isotopic_peak_cols = ['T_isotopic1_peaks', 'T_isotopic2_peaks', 'T_isotopic3_peaks']
    plot_dat = normalize_peaks(plot_dat, T_isotopic_peak_cols, 'T_isotopic1_peaks', new_prefix='norm_T_isotopic')
    plot_dat['cosine_similarity'] = plot_dat.apply(calculate_cosine_similarity, axis=1)
    plot_dat.to_csv('distribution_plot/isotopic_peak_from_exp_solico.csv', index=False)
    print(plot_dat.iloc[0,])
    plot_dat.apply(lambda row: plot_isotopic_peaks(row), axis=1)





# 示例：计算化学式 C6H12O6 的同位素峰
if __name__ == "__main__":
    isotopic_peak_distribution(input_folder = "ExampleDataset/Results/stagger_Nglycan_ExampleData/",
                               debug_file = "debug.txt")








# 绘制同位素峰
#plot_isotope_peaks(peaks)