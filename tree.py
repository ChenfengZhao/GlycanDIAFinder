from anytree import Node, RenderTree
from copy import deepcopy


class TreeNode:
    def __init__(self, value):
        self.value = value  # 原子名称
        self.children = []  # 子节点列表
        self.parent = None  # 父节点

    def add_child(self, node):
        node.parent = self
        self.children.append(node)

    def __repr__(self, level=0):
        ret = " " * (level * 4) + repr(self.value) + "\n"
        for child in self.children:
            ret += child.__repr__(level + 1)
        return ret

    def to_anytree_node(self):
        anytree_node = Node(self.value)
        for child in self.children:
            child_anytree = child.to_anytree_node()
            child_anytree.parent = anytree_node
        return anytree_node

    def to_expression(self):
        """
        将树结构转换回分子表达式。
        """
        expression = f"({self.value}"
        for child in self.children:
            expression += child.to_expression()
        expression += ")"
        return expression


def parse_molecule(expression):
    """
    解析分子表达式并构建多叉树。

    Args:
        expression (str): 分子表达式，如 "(N(F)(N(H(H)(H))))"

    Returns:
        TreeNode: 构建好的多叉树的根节点
    """

    def parse_helper(s, index):
        if index[0] >= len(s):
            return None

        # 当前字符应为 '('
        if s[index[0]] != '(':
            raise ValueError(f"Expected '(', found '{s[index[0]]}' at position {index[0]}")

        index[0] += 1  # 跳过 '('

        # 提取原子名称，支持多字符（如 'Cl', 'Br' 等）
        value = ''
        while index[0] < len(s) and s[index[0]].isalpha():
            value += s[index[0]]
            index[0] += 1

        if not value:
            raise ValueError(f"Expected atom name at position {index[0]}")

        node = TreeNode(value)

        # 解析所有子节点
        while index[0] < len(s) and s[index[0]] == '(':
            child = parse_helper(s, index)
            node.add_child(child)

        # 当前字符应为 ')'
        if index[0] >= len(s) or s[index[0]] != ')':
            raise ValueError(
                f"Expected ')', found '{s[index[0]] if index[0] < len(s) else 'EOF'}' at position {index[0]}")
        index[0] += 1  # 跳过 ')'

        return node

    index = [0]
    root = parse_helper(expression, index)

    if index[0] != len(expression):
        raise ValueError(f"Unexpected characters at the end of the expression starting at position {index[0]}")

    return root


def list_all_bonds(root):
    """
    列出分子中所有的化学键，每个键由父节点到子节点的路径表示。

    Args:
        root (TreeNode): 树的根节点

    Returns:
        list of list of int: 每个键的路径表示，如 [[0], [1], [1, 0], ...]
    """
    bonds = []

    def traverse(node, path):
        for idx, child in enumerate(node.children):
            child_path = path + [idx]
            bonds.append(child_path)
            traverse(child, child_path)

    traverse(root, [])
    return bonds


def clone_tree(node):
    """
    递归克隆树。

    Args:
        node (TreeNode): 要克隆的树节点

    Returns:
        TreeNode: 克隆后的树节点
    """
    new_node = TreeNode(node.value)
    for child in node.children:
        cloned_child = clone_tree(child)
        new_node.add_child(cloned_child)
    return new_node


def find_node_by_path(root, path):
    """
    根据路径查找节点。

    Args:
        root (TreeNode): 树的根节点
        path (list of int): 子节点索引路径，如 [1, 0]

    Returns:
        TreeNode: 目标节点，如果未找到则返回 None
    """
    current = root
    for idx in path:
        if idx < 0 or idx >= len(current.children):
            return None
        current = current.children[idx]
    return current


def break_bond(root, path_to_child):
    """
    根据路径断裂指定的键，将分子分裂成两个独立的部分。

    Args:
        root (TreeNode): 分子的根节点
        path_to_child (list of int): 子节点的路径索引，如 [1]

    Returns:
        tuple: (分裂后的第一部分根节点, 分裂后的第二部分根节点)
    """
    child = find_node_by_path(root, path_to_child)

    if child is None:
        raise ValueError(f"Path {path_to_child} not found in the molecule.")

    if child.parent is None:
        raise ValueError("Cannot break the bond at the root node.")

    # 移除子节点与其父节点的连接
    parent = child.parent
    parent.children.remove(child)
    child.parent = None

    return root, child  # 返回主树（已断裂）和被断裂的子树


def generate_all_products(root):
    """
    遍历所有可能的键断裂点，生成每个断裂后的产物。

    Args:
        root (TreeNode): 分子的根节点

    Returns:
        list of tuple: 每个断裂点对应的产物，如 [(part1, part2), ...]
    """
    bonds = list_all_bonds(root)
    products = []
    seen = set()

    for bond_path in bonds:
        # 克隆原始树
        cloned_root = clone_tree(root)
        # 断裂指定的键
        try:
            part1, part2 = break_bond(cloned_root, bond_path)
            # 将产物转换为字符串表示，以便去重
            part1_str = part1.to_expression()
            part2_str = part2.to_expression()
            # 为了避免顺序不同导致的重复，排序
            product = tuple(sorted([part1_str, part2_str]))
            if product not in seen:
                seen.add(product)
                products.append((part1, part2))
        except ValueError as e:
            print(f"错误断裂路径 {bond_path}: {e}")
            continue

    return products


def visualize_tree(anytree_node):
    """
    使用 anytree 可视化树结构。

    Args:
        anytree_node (anytree.Node): anytree 的节点对象
    """
    #for pre, fill, node in RenderTree(anytree_node):
      #  print(f"{pre}{node.name}")


def tree_to_expression(node):
    """
    将树结构转换回分子表达式。

    Args:
        node (TreeNode): 树的根节点

    Returns:
        str: 分子表达式
    """
    expression = f"({node.value}"
    for child in node.children:
        expression += tree_to_expression(child)
    expression += ")"
    return expression


def trees_to_expressions(trees):
    """
    将多个树结构转换回分子表达式。

    Args:
        trees (list of TreeNode): 树的根节点列表

    Returns:
        list of str: 分子表达式列表
    """
    return [tree.to_expression() for tree in trees]



def break_molecule(molecule_str):
    # 分子表达式
    # molecule_str = "(N(F)(N(H(H)(H))))"

    # 解析分子并构建多叉树
    root = parse_molecule(molecule_str)
    """
    print("原始分子结构:")
    print(root)
    """

    # 列出所有可能的键断裂路径
    bonds = list_all_bonds(root)

    """
    print("所有可能的键断裂路径 (子节点路径索引):")
    for idx, bond in enumerate(bonds):
        print(f"键 {idx + 1}: {bond}")
    """

    # 生成所有可能的产物
    products = generate_all_products(root)
   # print("\n所有可能的键断裂产物:")

    result = []
    for idx, (part1, part2) in enumerate(products):
        # print(f"\n产物 {idx + 1}:")
        # 转换树为分子表达式
        expr1 = part1.to_expression()
        expr2 = part2.to_expression()
        # print(f"分子 A: {expr1}")
        # print(f"分子 B: {expr2}")
        result.append([expr1, expr2])

    # 如果需要可视化，可以选择打印树形结构
    # 例如，对于每个产物，打印树形结构
    """
    print("\n详细分子产物结构:")
    for idx, (part1, part2) in enumerate(products):
        print(f"\n产物 {idx + 1}:")

        # 可视化分子 A
        anytree_a = part1.to_anytree_node()
        print("分子 A:")
        visualize_tree(anytree_a)

        # 可视化分子 B
        anytree_b = part2.to_anytree_node()
        print("分子 B:")
        visualize_tree(anytree_b)
    """
    return result

def check_exit(count, current_product, max_count):
    if not count:
        #print("done, exit")
        return True, get_total_elments(current_product, max_count)
    else:
        return False, []

def get_total_elments(current_product, max_count):
    total_products = []
    for k,v in current_product.items():
        if k == str(max_count):
            continue
        for i in v:
            for j in i:
                total_products.append(j)
    total_products = sorted(list(set(total_products)), key = lambda x: [len(x), x])

    total_elements = []
    from collections import Counter
    for i in total_products:
        element = ""
        freq = dict(Counter(i))
        for k,v in freq.items():
            if k == "(" or k == ")":
                continue
            element += f"{k}({v})"
        total_elements.append(element)
    return  total_elements

def molecule_breaker(molecule_str, count):
    if not  count <= molecule_str.count("(") - 1:
        raise ValueError("break counts larger than keys")
    max_count = count
    current_molecule = molecule_str
    current_product = dict({str(count):[molecule_str]})
   # print(f"the {max_count-count} round")
   # print(current_product[str(count)][0])
   # print("="*30)

    count -=1
    current_product[str(count)] = []
    current_product[str(count)] = break_molecule(current_product[str(count+1)][0])
    current_product[str(count)] = list(set(tuple(sublist) for sublist in current_product[str(count)]))
    #print(f"the {max_count-count} round")
    for i in current_product[str(count)]:
        #print(i)
        #print("="*30)
        status, total_element = check_exit(count, current_product, max_count)
    if status:
        return total_element

    count -= 1
    current_product[str(count)] = []
    for raw_subproducts in current_product[str(count+1)]:
        subproducts = sorted(raw_subproducts, key = lambda x: len(x))
        res = []
        for i in range(len(subproducts)):
            if subproducts[i].count("(") == 1:
                continue
            for subres in break_molecule(subproducts[i]):
                if 0 <= i < len(subproducts):
                    others = subproducts[:i] + subproducts[i + 1:]
                res.append(sorted(subres + others))
        current_product[str(count)].extend(res)
    current_product[str(count)] = list(set(tuple(sublist) for sublist in current_product[str(count)]))
   # print(f"the {max_count - count} round")
    for i in current_product[str(count)]:
        #print(i)
        #print("="*30)
        status, total_element = check_exit(count, current_product, max_count)
    if status:
        return total_element

    while count :
        count -= 1
        current_product[str(count)] = []
        for raw_subproducts in current_product[str(count+1)]:
            subproducts = sorted(raw_subproducts, key = lambda x: len(x))
            res = []
            for i in range(len(subproducts)):
                if subproducts[i].count("(") == 1:
                    continue
                for subres in break_molecule(subproducts[i]):
                    if 0 <= i < len(subproducts):
                        others = subproducts[:i] + subproducts[i + 1:]
                    res.append(sorted(subres + others))
            current_product[str(count)].extend(res)
        current_product[str(count)] = list(set(tuple(sublist) for sublist in current_product[str(count)]))
        #print(f"the {max_count - count} round")
        for i in current_product[str(count)]:
            #print(i)
            #print("="*30)
            status, total_element = check_exit(count, current_product, max_count)
    if status:
        return total_element

import argparse

def argparsor():
    parser = argparse.ArgumentParser(description="一个简单的计算器，支持加法和减法。")
    parser.add_argument("--count", type=int, default=3, required=False, help="")  # 第一个数字
    parser.add_argument("--molecule", type=str, default="(N(F)(N(H(H)(H))))", required=False, help="第二个数字")  # 第二个数字
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = argparsor()
    molecule_str = args.molecule
    count = args.count

    res = molecule_breaker(molecule_str, count)
   # print(res)

    """
    python tree.py --count 5 -- molecule "(N(F)(N(H(H)(H))))"
    """
