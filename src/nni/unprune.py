import pandas as pd
import numpy as np

import json
import argparse

from collections import deque
from utilities import *

def tree_to_newick_eq_classes(T, eq_class_dict, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            subgs.append(tree_to_newick_eq_classes(T, eq_class_dict, root=child))
        else:
            eq_class = eq_class_dict[child]
            if len(eq_class) > 1:
                subgs.append("(" + ','.join(map(str, eq_class)) + ")")
            else:
                subgs.append(eq_class[0])
    return "(" + ','.join(map(str, subgs)) + ")"

def parse_args():
    p = argparse.ArgumentParser(
        description="Reattaches equivalence classes to generated Newick tree."
    )

    p.add_argument("eq_classes", help="Equivalence classes dictionary")
    p.add_argument("newick_tree", help="Newick tree generated from a pruned character matrix")

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    with open(args.eq_classes, "r") as f:
        eq_class_dict = json.load(f)

    newick_tree = from_newick_get_nx_tree(args.newick_tree)
    newick_tree = tree_to_newick_eq_classes(newick_tree, eq_class_dict)
    print(f"{newick_tree};")
