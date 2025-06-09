import pandas as pd


def build_tree(paths):
    """
    Builds taxonomy tree.

    :param paths:
    :return:
    """
    root = {}
    for path in paths:
        try:
            parts = [p.strip() for p in path.split(";") if p.strip()]
        except AttributeError as e:
            print(path)
        current = root
        for part in parts:
            if part not in current:
                current[part] = {}
            current = current[part]
    return root


def print_tree(tree, prefix="", file=None):
    """
    Prints taxonomy tree to a standard output or a file.

    :param tree:
    :param prefix:
    :param file:
    :return:
    """
    keys = list(tree.keys())
    for i, key in enumerate(keys):
        connector = "└── " if i == len(keys) - 1 else "├── "
        line = prefix + connector + key
        if file:
            file.write(line + "\n")
        else:
            print(line)
        next_prefix = prefix + ("    " if i == len(keys) - 1 else "│   ")
        print_tree(tree[key], next_prefix, file)


def create_taxonomy(input_path, output_path=None):
    """
    Creates file with taxonomy from raw taxonomy file.

    :param input_path:
    :param output_path:
    :return:
    """
    if output_path is None:
        output_path = f"{input_path.split('.tsv')[0]}_print.txt"

    df = pd.read_csv(input_path, sep="\t", names=["taxid", "taxonomy"])
    raw_taxa = df.taxonomy.tolist()
    tree = build_tree(raw_taxa)

    with open(output_path, "w", encoding="utf-8") as outfile:
        print_tree(tree, file=outfile)
