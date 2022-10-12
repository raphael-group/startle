import pandas as pd
import numpy as np

import json
import argparse

from collections import deque
from utilities import *

"""
Computes the equivalence classes of the rows
of the input character matrix. 

It outputs both the set of equivalence classes
and pruned character matrix.
"""
def compute_equivalence_classes(df_character_matrix):
    df_character_matrix = df_character_matrix.copy()
    df_character_matrix = df_character_matrix.astype(int)

    ncells = len(df_character_matrix)
    binary_col_dict = {}
    for column in df_character_matrix.columns:
        state_list = list(df_character_matrix[column].unique())
        for s in state_list:
            if s != -1 and s != 0:
                state_col = np.zeros((ncells))
                state_col[df_character_matrix[column] == s] = 1
                state_col[df_character_matrix[column] == -1] = -1

                binary_col_dict[f'{column}_{s}'] = state_col

    df_binary = pd.DataFrame(binary_col_dict, index = df_character_matrix.index, dtype=int)
     
    one_set_dict = {}
    zero_set_dict = {}
    for cell in df_binary.index:
        one_set_dict[cell] = set([x for x, y in df_binary.loc[cell].items() if y == 1])
        zero_set_dict[cell] = set([x for x, y in df_binary.loc[cell].items() if y == 0])

    eq_class_list = []
    for cell in df_binary.index:
        # check if it belongs to an existing equivalence class
        flag = 0
        for eq_class in eq_class_list:
            leader_cell = eq_class[0]
            if (one_set_dict[leader_cell].issuperset(one_set_dict[cell]) and
                zero_set_dict[leader_cell].issuperset(zero_set_dict[cell])):
                eq_class.append(cell)
                flag = 1
                break
            elif (one_set_dict[cell].issuperset(one_set_dict[leader_cell]) and
                  zero_set_dict[cell].issuperset(zero_set_dict[leader_cell])):
                eq_class.appendleft(cell)
                flag = 1
                break

        if flag == 0:
            eq_class_list.append(deque([cell]))

    # map from leader to equivalence class
    eq_class_dict = {eq_class[0]: list(eq_class) for eq_class in eq_class_list}

    # drop non-leader cells in each equivalence class
    for eq_class in eq_class_list:
        for cell in eq_class:
            if cell == eq_class[0]: continue
            df_character_matrix = df_character_matrix.drop(cell, axis=0)

    return eq_class_dict, df_character_matrix

def prune_tree(tree, eq_class_leaders):
    T = tree.copy()
    while True:
        leaves_to_prune = [v for v in T.nodes if is_leaf(T, v) and v not in eq_class_leaders]
        if not leaves_to_prune: break
        for leaf in leaves_to_prune:
            T.remove_node(leaf)

    return T

def parse_args():
    p = argparse.ArgumentParser(
        description=("Partitions a character matrix into sets of cells"
                     "belonging to different equivalence classes.")
    )

    p.add_argument("character_matrix", help="Character matrix (TSV)")
    p.add_argument("-n", "--newick_tree", help="Newick tree to prune")
    p.add_argument("output", help="Output files prefix")

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    character_matrix = pd.read_csv(args.character_matrix, index_col=[0], sep='\t', dtype=str)
    character_matrix = character_matrix.replace('-', '-1')

    eq_class_dict, character_matrix = compute_equivalence_classes(character_matrix)
    with open(f"{args.output}_eq_classes.json", "w") as f:
        json.dump(eq_class_dict, f)

    # if newick
    if args.newick_tree:
        tree = from_newick_get_nx_tree(args.newick_tree)
        tree = prune_tree(tree, set(eq_class_dict.keys()))

        newick_tree = tree_to_newick(tree)
        with open(f"{args.output}_pruned_tree.newick", 'w') as f:
            f.write(f"{newick_tree};")

    character_matrix.to_csv(f"{args.output}_pruned_character_matrix.txt", sep='\t')
