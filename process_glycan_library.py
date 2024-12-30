def replace_first_non_bracket_character(line):
    """
    替换每一行的第一个非括号字符：
    - 如果是 'N'，替换为 'J'
    - 如果是 'H'，替换为 'Q'
    """
    result = []
    replaced = False  # 标记是否已经替换过第一个非括号字符
    for char in line:
        if not replaced and char not in "()":  # 找到第一个非括号字符
            if char == "N":
                result.append("J")  # 替换 N 为 J
            elif char == "H":
                result.append("Q")  # 替换 H 为 Q
            else:
                result.append(char)  # 保留其他非括号字符
            replaced = True
        else:
            result.append(char)  # 保留括号和其他字符
    return "".join(result)


def process_gdb_file(input_file, output_file):
    """
    读取 gdb 文件，解析括号平衡的部分，替换第一个非括号字符，并写入新文件。
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        next(infile)  # 跳过第一行
        for line in infile:
            line = line.strip()  # 去掉行首尾的空格

            buffer = ""  # 缓存当前未处理的内容
            left_count = 0  # 左括号计数
            right_count = 0  # 右括号计数

            for char in line:
                buffer += char
                if char == "(":
                    left_count += 1
                elif char == ")":
                    right_count += 1

                # 当左右括号数目相等时，提取完整部分
                if left_count == right_count and left_count > 0:
                    # 替换第一个非括号字符
                    processed_line = replace_first_non_bracket_character(buffer)
                    outfile.write(processed_line + "\n")  # 写入新文件
                    buffer = ""  # 清空缓存
                    left_count = 0  # 重置计数
                    right_count = 0

            # 如果还有未处理的内容，继续检查
            if buffer:
                while buffer:
                    new_buffer = ""
                    left_count = 0
                    right_count = 0
                    for char in buffer:
                        new_buffer += char
                        if char == "(":
                            left_count += 1
                        elif char == ")":
                            right_count += 1
                        if left_count == right_count and left_count > 0:
                            # 替换第一个非括号字符
                            processed_line = replace_first_non_bracket_character(new_buffer)
                            outfile.write(processed_line + "\n")
                            buffer = buffer[len(new_buffer):]
                            break
                    else:
                        # 如果无法匹配，跳出循环
                        break


if __name__ == "__main__":
    # 输入文件和输出文件路径
    O_glycan_input_file = "glyco_library/pGlyco-O-Glycan.gdb"
    O_glycan_output_file = "glyco_library/O-Glycan.txt"
    Human_N_glycan_input_file = "glyco_library/pGlyco-N-Human.gdb"
    Human_N_glycan_output_file = "glyco_library/N-Human.txt"
    Mouse_N_glycan_input_file = "glyco_library/pGlyco-N-Mouse.gdb"
    Mouse_N_glycan_output_file = "glyco_library/N-Mouse.txt"
    Plant_N_glycan_input_file = "glyco_library/pGlyco-N-Plant.gdb"
    Plant_N_glycan_output_file = "glyco_library/N-Plant.txt"

    process_gdb_file(O_glycan_input_file, O_glycan_output_file)
    process_gdb_file(Human_N_glycan_input_file, Human_N_glycan_output_file)
    process_gdb_file(Mouse_N_glycan_input_file, Mouse_N_glycan_output_file)
    process_gdb_file(Plant_N_glycan_input_file, Plant_N_glycan_output_file)