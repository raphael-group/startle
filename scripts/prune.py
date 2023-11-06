import pandas as pd
import numpy as np

import json
import argparse

from utilities import *

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
