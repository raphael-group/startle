import pandas as pd
import numpy as np

import json
import argparse

from collections import deque
from utilities import *

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
