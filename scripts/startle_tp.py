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

from collections import defaultdict
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
        print(line)
        m = re.match(r'Step (\d+) -- number of violations: (\d+)', line)
        if m is not None:
            logger.info(f"ILP Step {m.group(1)} has {m.group(2)} violations.")
    logger.info("ILP complete")

    inferred_tree_file = f"{subprocess_dir}/startle_tree.newick"

    inferred_tree = from_newick_get_nx_tree(inferred_tree_file)
    shutil.rmtree(subprocess_dir)
    return inferred_tree

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

        parsimony, parsimony_scores, labeling = small_parsimony(mutation_prior_dict, lenti_tree, lenti_cm)
