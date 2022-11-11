import argparse
import shutil
import os
import re
import tempfile
import io
import sys

import subprocess as sp
import pandas as pd
import networkx as nx

from collections import defaultdict, deque
from loguru import logger

sys.path.append("src/nni")

from clucknni import small_parsimony
from utilities import *

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Placement algorithm for lentivirally barcoded " +
            "lineage tracing data."
        )
    )

    p.add_argument(
        "character_matrix",
        help="Character matrix in CSV format."
    )

    p.add_argument(
        "mutation_priors",
        help="Mutation priors in CSV format."
    )

    return p.parse_args()

"""
Runs the Startle-ILP and returns the inferred
tree. This routine makes heavy use of the subprocess
module and has many implicit dependencies (perl/C++/etc).
"""
def startle_ilp(character_matrix, mutation_priors):
    subprocess_dir = tempfile.mkdtemp(dir=".")

    character_matrix_file = f"{subprocess_dir}/character_matrix.csv"
    character_matrix.to_csv(character_matrix_file)

    mutation_priors_file = f"{subprocess_dir}/mutation_priors.csv"
    mutation_priors.to_csv(mutation_priors_file)
    
    cmd = [
        "perl", "scripts/startle_ilp.pl", 
        "-m", mutation_priors_file,
        "-c", character_matrix_file,
        "-o", subprocess_dir
    ]

    logger.info("Started ILP")
    proc = sp.Popen(cmd, stderr=sp.PIPE)
    for line in iter(proc.stderr.readline, b""):
        line = line.decode("utf-8")
        m = re.match(r'Step (\d+) -- number of violations: (\d+)', line)
        if m is not None:
            logger.info(f"ILP Step {m.group(1)} has {m.group(2)} violations.")
    logger.info("ILP complete")

    inferred_tree_file = f"{subprocess_dir}/startle_tree.newick"

    inferred_tree = from_newick_get_nx_tree(inferred_tree_file)
    shutil.rmtree(subprocess_dir)
    return inferred_tree

"""
Cell is a row of the character matrix.
"""
def star_homoplasy_distance(tree, labeling, cell):
    # TODO: use weighted distance instead
    def distance(node, cell):
        d, visit_child = 0, True
        for (character, character_labeling) in labeling.items():
            node_state = character_labeling[node]
            cell_state = cell[character]

            if node_state != cell_state and cell_state != '-1':
                d += 1

            if node_state != '0' and cell_state == '0':
                visit_child = False

        return d, visit_child

    min_distance = np.inf

    stack = deque(['root'])
    while stack:
        current_node = stack.pop()
        d, visit_child = distance(current_node, cell)
        if d < min_distance:
            min_distance = d
        if not visit_child: continue
        for child in tree[current_node]:
            stack.append(child)

    return min_distance

if __name__ == "__main__":
    args = parse_args()
    
    input_character_matrix = pd.read_csv(args.character_matrix, index_col=0)
    input_mutation_priors = pd.read_csv(args.mutation_priors, index_col=0)

    mutation_prior_dict = defaultdict(dict)
    for (character, row) in input_mutation_priors.iterrows():
        character_idx = int(character[1:])
        mutation_prior_dict[character_idx][str(int(row['state']))] = row['probability']

    lenti_groups = input_character_matrix['c0'].unique()
    lenti_trees  = {}
    for lenti_group in lenti_groups:
        if lenti_group == -1: continue

        lenti_cm = input_character_matrix[input_character_matrix['c0'] == lenti_group].copy().drop(columns=['c0'])
        lenti_tree = startle_ilp(lenti_cm, input_mutation_priors)
        lenti_cm = lenti_cm.astype(str)

        # labeling is a function from (character, nodes) -> state
        parsimony, parsimony_scores, labeling = small_parsimony(mutation_prior_dict, lenti_tree, lenti_cm)
        lenti_trees[lenti_group] = (labeling, lenti_tree)

        # test_cell = input_character_matrix[input_character_matrix['c0'] == -1].copy().iloc[0].astype(str)
        # test_cell = {c: s for c, s in test_cell.items()}
        # cell_distance = star_homoplasy_distance(lenti_tree, labeling, test_cell)
        # print(cell_distance)
