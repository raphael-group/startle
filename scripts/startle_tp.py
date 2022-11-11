import argparse

import shutil
import os
import re
import tempfile
import io
import subprocess as sp
import sys
from loguru import logger

import pandas as pd
import networkx as nx

sys.path.append("src/nni")

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
        m = re.match(r'Step (\d+) -- number of violations: (\d+)', line.decode("utf-8"))
        if m is not None:
            logger.info(f"ILP Step {m.group(1)} has {m.group(2)} violations.")
    logger.info("ILP complete")

    inferred_tree_file = f"{subprocess_dir}/startle_tree.newick"

    inferred_tree = from_newick_get_nx_tree(inferred_tree_file)
    shutil.rmtree(subprocess_dir)
    return inferred_tree

if __name__ == "__main__":
    args = parse_args()

    # split into lenti groups
    
    input_character_matrix = pd.read_csv(args.character_matrix, index_col=0)
    input_mutation_priors = pd.read_csv(args.mutation_priors, index_col=0)

    lenti_groups = input_character_matrix['c0'].unique()
    lenti_trees  = {}
    for lenti_group in lenti_groups:
        if lenti_group == -1: continue

        lenti_cm = input_character_matrix[input_character_matrix['c0'] == lenti_group].copy().drop(columns=['c0'])
        lenti_tree = startle_ilp(lenti_cm, input_mutation_priors)
