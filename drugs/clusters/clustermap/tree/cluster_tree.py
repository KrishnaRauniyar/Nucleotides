import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import argparse

def load_and_process_data(file_path):
    df = pd.read_csv(file_path, header=None)
    df1 = df[0].str.split(';', expand=True)
    df = df.drop(0, axis=1)
    merged_df = pd.concat([df1, df], axis=1)
    ff = merged_df.set_index(merged_df[0]).T.reset_index(drop=True).rename_axis(None, axis=1)
    ff1 = ff.drop(index=0)
    ff2 = ff1.reset_index(drop=True)
    ff2.index = ff2.columns
    return ff2

def perform_clustering(data, sample_size=100):
    sampled_data = data.sample(n=sample_size, random_state=42)
    Z = linkage(sampled_data, method='average')
    return Z, sampled_data

def create_newick(node, parent_dist, leaf_names, newick=''):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = create_newick(node.get_left(), node.dist, leaf_names, newick)
        newick = create_newick(node.get_right(), node.dist, leaf_names, ",%s" % (newick))
        newick = "(%s" % (newick)
        return newick

def get_group(leaf_name, group_colors):
    suffixes = sorted(group_colors.keys(), key=len, reverse=True)
    for suffix in suffixes:
        if leaf_name.endswith(suffix):
            return suffix
    return 'U'

def create_tree(Z, data, group_colors):
    tree = to_tree(Z)
    newick_str = create_newick(tree, tree.dist, list(data.index)) + ";"
    ete_tree = Tree(newick_str, format=1)

    ts = TreeStyle()
    ts.mode = "c"
    ts.arc_start = -180
    ts.arc_span = 360

    for n in ete_tree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        n.set_style(nstyle)

    for leaf in ete_tree:
        leaf_group = get_group(leaf.name, group_colors)
        leaf_style = NodeStyle()
        leaf_style["fgcolor"] = group_colors.get(leaf_group)
        leaf_style["size"] = 10
        leaf.set_style(leaf_style)

        name_face = TextFace(leaf.name, fsize=10, fgcolor=group_colors.get(leaf_group))
        leaf.add_face(name_face, column=0, position="branch-right")

        leaf.name = ''

    return ete_tree, ts

def main(input):
    file_path = input
    group_colors = {
        'A': 'teal',
        'DA': 'orange',
        'C': 'violet',
        'DC': 'green',
        'G': 'blue',
        'DG': 'purple',
        'DT': 'red',
        'U': 'brown'
    }
    
    data = load_and_process_data(file_path)
    Z, sampled_data = perform_clustering(data)
    ete_tree, ts = create_tree(Z, sampled_data, group_colors)
    ete_tree.render("tree.png", w=4000, h=4000, tree_style=ts)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot tree for hierarchial clustering")
    parser.add_argument('-p', '--input_path', type=str, required=True, help="Generalized CSV input path")
    args = parser.parse_args()
    main(args.input_path)
