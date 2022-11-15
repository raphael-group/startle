import argparse
import math
import random
import pickle
import re
import glob
import sys

import pandas as pd
import networkx as nx
import numpy as np
import multiprocessing as mp

from collections import deque
from loguru import logger
from tqdm import tqdm
from funcy import *
from utilities import *

"""
Computes small parsimony under the Star Camin-Sokal model
for a fixed character.
"""
def small_parsimony_fixed_character(T, character_series, start_node, weight_function, root_state='0', missing_state='-1'):
    stack = deque([start_node]) # stack to simulate DFS search

    # these are vertex level values
    labeling = {}
    parsimony_score = {}

    while stack:
        node = stack.pop()

        if is_leaf(T, node):
            labeling[node] = character_series[node]
            parsimony_score[node] = 0
            continue

        # i.e. if all children have been visited
        if all(child in labeling for child in T[node]):
            parsimony_score[node] = 0
            child_states = set()
            for child in T[node]:
                child_states.add(labeling[child])
                parsimony_score[node] += parsimony_score[child]

            if len(child_states - {missing_state}) > 1:
                labeling[node] = root_state
                for state in child_states:
                    if state != root_state and state != missing_state:
                        parsimony_score[node] += weight_function(state)
            else:
                if len(child_states - {missing_state}) == 0:
                    labeling[node] = missing_state
                else:
                    labeling[node] = list(child_states - {missing_state})[0]

            continue

        # if all children have not been visited, push current node on stack and then
        # push all children on stack so that next time we visit the current node,
        # we will have visited all children
        stack.append(node)
        for child in T[node]:
            stack.append(child)

    return labeling, parsimony_score

"""
Computes small parsimony score under Star Camin-Sokal model,
assumes a 0 rooted vertex.
"""
def small_parsimony(mutation_prior_dict, T, df, weighted=True):
    parsimony_scores, labelings = {}, {}

    parsimony = 0
    for character in df:
        character_idx = int(character[1:]) # i.e. 'r7' -> 7

        if weighted:
            def weight_function(s):
                return -np.log(mutation_prior_dict[character_idx][s])
        else:
            weight_function = lambda s: 1

        labeling, parsimony_score = small_parsimony_fixed_character(T, df[character], 'root', weight_function)

        root_state = labeling['root']
        if root_state == '-1':
            root_state = '0'
            labeling['root'] = root_state

        if root_state != '0':
            parsimony_score['root'] += weight_function(root_state)

        parsimony += parsimony_score['root']

        parsimony_scores[character] = parsimony_score
        labelings[character] = labeling

    return parsimony, parsimony_scores, labelings

"""
Recomputes the parsimony score efficiently for a particular
NNI move on T without having to traverse the entire tree. Time
complexity is O(d) where d is the depth of the vertex v in the NNI
move.
"""
def recompute_small_parsimony(mutation_prior_dict, T, df, parsimony_scores, labelings, nni_move, missing_state='-1', root_state='0'):
    (u, w), (v, z) = nni_move

    # short circuit if no change when swapping
    if all(labelings[char][w] == labelings[char][z] for char in df):
        return sum(parsimony_scores[char]['root'] for char in df)

    parent = v
    path_to_root = deque()
    while True:
        path_to_root.append(parent)

        preds = list(T.predecessors(parent))
        assert len(preds) <= 1
        if len(preds) == 0:
            break
        else:
            parent = preds[0]

    # swap e1 and e2
    T.remove_edge(u, w) 
    T.remove_edge(v, z)
    T.add_edge(v, w)
    T.add_edge(u, z) 

    # save labelings and parsimony scores
    saved_labelings = {}
    saved_parsimony_scores = {}
    for character in df:
        character_idx = int(character[1:]) # i.e. 'r7' -> 7
        weight_function = lambda s: -np.log(mutation_prior_dict[character_idx][s])

        labeling = labelings[character]
        parsimony_score = parsimony_scores[character]

        saved_parsimony_scores[character] = {}
        saved_labelings[character] = {}

        # perform saving
        for node in path_to_root:
            saved_parsimony_scores[character][node] = parsimony_score[node]
            saved_labelings[character][node] = labeling[node]

        # update labeling and parsimony scores by walking
        # up path to root
        for node in path_to_root:
            parsimony_score[node] = 0

            child_states = set()
            for child in T[node]:
                child_states.add(labeling[child])
                parsimony_score[node] += parsimony_score[child]

            if len(child_states - {missing_state}) > 1:
                labeling[node] = root_state
                for state in child_states:
                    if state != root_state and state != missing_state:
                        parsimony_score[node] += weight_function(state)
            else:
                if len(child_states - {missing_state}) == 0:
                    labeling[node] = missing_state
                else:
                    labeling[node] = list(child_states - {missing_state})[0]

    parsimony = 0
    for character in df:
        character_idx = int(character[1:]) # i.e. 'r7' -> 7
        weight_function = lambda s: -np.log(mutation_prior_dict[character_idx][s])

        root_state = labelings[character]['root']
        if root_state != '0':
            parsimony_scores[character]['root'] += weight_function(root_state)

        parsimony += parsimony_scores[character]['root']

    # restore labelings and parsimony_scores
    for character in df:
        for node in path_to_root:
            parsimony_scores[character][node] = saved_parsimony_scores[character][node]
            labelings[character][node] = saved_labelings[character][node]

    # return T to original state
    T.add_edge(u, w) 
    T.add_edge(v, z)
    T.remove_edge(v, w)
    T.remove_edge(u, z) 

    return parsimony

"""
Collapses mutationless edges as defined by an
optimal Star Camin-Sokal labeling of
the internal vertices.
"""
def collapse_mutationless_edges(T, df):
    T = T.copy()

    # property of labeling:
    #   if a vertex v is labeled by a missing state,
    #   then all children of v are labeled by a missing
    #   state
    labeling_dict, _ = compute_labeling_count(T, df)

    while True:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))

        no_contractions = True
        for (u, v) in T.edges:
            if is_leaf(T, v):
                continue

            if not all((labeling[u] == labeling[v] or labeling[v] == '-1') for labeling in labeling_dict.values()):
                continue

            T.remove_edge(u, v)
            for w in list(T[v].keys()):
                T.remove_edge(v, w)
                T.add_edge(u, w)

            T.remove_node(v)

            no_contractions = False
            break

        if no_contractions:
            break

    return T

"""
Utility function for computing two types of scores
"""
def compute_parsimony(mutation_prior_dict, T, character_matrix):
    unweighted_score, _, _ = small_parsimony(mutation_prior_dict, T, character_matrix, weighted=False)
    weighted_score, _, _ = small_parsimony(mutation_prior_dict, T, character_matrix, weighted=True)

    return {
        "parsimony": unweighted_score,
        "weighted_parsimony": weighted_score
    }

"""
Stochastically pertubs T using NNI operations.
"""
def stochastic_nni(T, aggression=0.30):
    T = T.copy()

    internal_edges = [(u, v) for (u, v) in T.edges if not is_leaf(T, v)]
    num_perturbations = math.floor(len(internal_edges) * aggression)
    count = 0
    while count < num_perturbations:
        internal_edges = [(u, v) for (u, v) in T.edges if not is_leaf(T, v)]
        u, v = random.sample(list(internal_edges), 1)[0]

        if is_leaf(T, v):
            continue

        u_children = list(set(T[u].keys()) - set([v]))
        v_children = list(T[v].keys())

        u_edges = [(u, w) for w in u_children]
        v_edges = [(v, w) for w in v_children]
        if not u_edges or not v_edges:
            continue

        u, w = random.sample(u_edges, 1)[0]
        v, z = random.sample(v_edges, 1)[0]

        # swap e1 and e2
        T.remove_edge(u, w) 
        T.remove_edge(v, z)
        T.add_edge(v, w)
        T.add_edge(u, z) 

        count += 1

    return T

"""
Performs all NNI operation(s) using the edge e of T,
returning only those that improve the score.
"""
def greedy_nni(current_score, scoring_function, T, e):
    u, v = e[0], e[1]

    if is_leaf(T, v):
        return []

    u_children = list(set(T[u].keys()) - set([v]))
    v_children = list(T[v].keys())

    u_edges = [(u, w) for w in u_children]
    v_edges = [(v, w) for w in v_children]

    nni_moves = []
    for (u, w) in u_edges:
        for (v, z) in v_edges:
            score = scoring_function.rescore(T, ((u, w), (v, z)))
            if score < current_score:
                nni_moves.append({
                    "move": ((u, w), (v, z)),
                    "score": score
                })

    return nni_moves

def hill_climb_on_edges(scoring_function, current_tree, edge_set):
    current_score = scoring_function.most_recent_score
    nni_best_score, nni_best_move = current_score, None

    positive_moves = []
    for e in edge_set:
        # loop invariant: current_tree remains identical
        nni_moves = greedy_nni(nni_best_score, scoring_function, current_tree, e)

        for nni_move in nni_moves:
            if nni_move["score"] > current_score: continue
            positive_moves.append(nni_move)
            if nni_move["score"] < nni_best_score:
                nni_best_move = nni_move["move"]
                nni_best_score = nni_move["score"]

    if nni_best_move is None:
        return current_tree, nni_best_score

    # # perform all non-conflicting moves
    # conflicting_vertices = set()
    # non_conflicting_positive_moves = []
    # for nni_move in positive_moves:
    #     (u, w), (v, z) = nni_move["move"]
    #     if u in conflicting_vertices or v in conflicting_vertices:
    #         continue

    #     non_conflicting_positive_moves.append(nni_move)
    #     conflicting_vertices.update([u, w, v, z])

    #     current_tree.remove_edge(u, w) 
    #     current_tree.remove_edge(v, z)
    #     current_tree.add_edge(v, w)
    #     current_tree.add_edge(u, z) 

    # # check if all moves is better than best move
    # non_conflicting_score = scoring_function.score(current_tree)
    # if non_conflicting_score > nni_best_score:
    #     logger.info(f"Many simultaneous NNIs applied: {non_conflicting_score} > {nni_best_score}")
    #     return current_tree, non_conflicting_score

    # # undo all moves
    # for nni_move in non_conflicting_positive_moves:
    #     (u, w), (v, z) = nni_move["move"]

    #     current_tree.add_edge(u, w) 
    #     current_tree.add_edge(v, z)
    #     current_tree.remove_edge(v, w)
    #     current_tree.remove_edge(u, z) 

    # else use only one move
    (u, w), (v, z) = nni_best_move 

    current_tree.remove_edge(u, w) 
    current_tree.remove_edge(v, z)
    current_tree.add_edge(v, w)
    current_tree.add_edge(u, z) 

    return current_tree, nni_best_score

"""
Performs hill climbing until a local
optima has been reached.
"""
def hill_climb(T, scoring_function, threads=8):
    overall_best_tree = T
    processor_pool = mp.Pool(threads)
    
    chunk_size = math.ceil(len(list(overall_best_tree.edges)) / threads)

    while True:
        overall_best_score = scoring_function.score(overall_best_tree)
        all_edges = list(overall_best_tree.edges)
        random.shuffle(all_edges)
        
        edge_sets = chunks(chunk_size, all_edges)

        proc_arguments = [(scoring_function, overall_best_tree, edge_set) for edge_set in edge_sets]
        hill_climb_trees = processor_pool.starmap(hill_climb_on_edges, proc_arguments)

        nni_best_tree, nni_best_score = overall_best_tree, overall_best_score
        for hill_climb_tree, hill_climb_score in hill_climb_trees:
            if hill_climb_score < nni_best_score:
                nni_best_tree, nni_best_score = hill_climb_tree, hill_climb_score

        if nni_best_score == overall_best_score:
            break

        logger.info(f"Score improved from {overall_best_score} to {nni_best_score}")
        
        overall_best_tree = nni_best_tree

    processor_pool.close()
    return overall_best_tree

"""
Converts a non-binary phylogenetic tree to a binary
phylogenetic tree on the same leaf set. 
"""
def arbitrarily_resolve_polytomies(T):
    T = T.copy()

    largest_clade_idx = max(int(node[6:]) for node in T.nodes if 'clade' in node)
    clade_idx = largest_clade_idx + 1

    for u in list(T.nodes):
        children = list(T[u].keys())

        if len(children) <= 2:
            continue

        for v in children:
            T.remove_edge(u, v)

        us = []
        for i in range(len(children) - 2):
            new_u = f'clade_{clade_idx}'
            T.add_node(new_u)
            us.append(new_u)
            clade_idx += 1

            if i == 0:
                T.add_edge(u, new_u)
            else:
                T.add_edge(us[i - 1], new_u)

        us = [u] + us 
        assert len(us) + 1 == len(children)

        T.add_edge(us[-1], children[-1])
        for w, v in zip(us, children[:-1]):
            T.add_edge(w, v)

    return T

def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("seed_tree", help="Seed tree in Newick format.")
    p.add_argument("character_matrix", help="Character matrix.")
    p.add_argument("-d", help="Mutation prior pickle file.")
    p.add_argument("-m", help="Mutation prior CSV table.")

    p.add_argument(
        "--iterations",
        help="Number of iterations to run stochastic hill climbing before giving up.",
        type=int,
        default=250
    )

    p.add_argument(
        "--mode",
        choices=['collapse', 'score', 'infer'],
        default='infer',
        help="Different modes have different functionality."
    )

    p.add_argument("--output", help="Output file for newick tree.", required=True)

    return p.parse_args()

def load_mutation_prior_dict(args):
    if args.m:
        mutation_prior_dict = defaultdict(dict)
        input_mutation_priors = pd.read_csv(args.m, index_col=0)
        for (character, row) in input_mutation_priors.iterrows():
            character_idx = int(character[1:])
            mutation_prior_dict[character_idx][str(int(row['state']))] = row['probability']

        return mutation_prior_dict
    elif args.e:
        with open(args.e, 'rb') as f:
            return pickle.load(f)

def score_mode(args):
    seed_tree = from_newick_get_nx_tree(args.seed_tree)

    character_matrix = pd.read_csv(args.character_matrix, index_col=[0], sep='\t', dtype=str)
    character_matrix = character_matrix.replace('-', '-1')

    mutation_prior_dict = load_mutation_prior_dict(args)

    seed_parsimony = compute_parsimony(mutation_prior_dict, seed_tree, character_matrix) 
    logger.info(f"Seed tree parsimony: {seed_parsimony}")

def collapse_mode(args):
    input_tree = from_newick_get_nx_tree(args.seed_tree)

    character_matrix = pd.read_csv(args.character_matrix, index_col=[0], sep='\t', dtype=str)
    character_matrix = character_matrix.replace('-', '-1')

    mutation_prior_dict = load_mutation_prior_dict(args)

    input_tree_contracted = collapse_mutationless_edges(input_tree, character_matrix)
    newick_tree = tree_to_newick(input_tree_contracted)
    with open(args.output, 'w') as f:
        f.write(f"{newick_tree};")

"""
Need to make ScoringFunction a class instead of a function object 
so that it can be pickled and can mantain internal state
for efficient rescoring.
"""
class ScoringFunction:
    def __init__(self, mutation_prior_dict, character_matrix):
        self.mutation_prior_dict = mutation_prior_dict
        self.character_matrix = character_matrix
    
    def score(self, T):
        score, parsimony_scores, labelings = small_parsimony(self.mutation_prior_dict, T, self.character_matrix)
        self.parsimony_scores = parsimony_scores
        self.labelings = labelings
        self.most_recent_score = score
        return score

    def rescore(self, T, nni_move):
        score = recompute_small_parsimony(
            self.mutation_prior_dict, T, self.character_matrix,
            self.parsimony_scores, self.labelings, nni_move
        )

        return score

def infer_mode(args):
    random.seed(73)

    seed_tree = from_newick_get_nx_tree(args.seed_tree)
    seed_tree = arbitrarily_resolve_polytomies(seed_tree)

    character_matrix = pd.read_csv(args.character_matrix, index_col=[0], dtype=str)
    character_matrix = character_matrix.replace('-', '-1')
    character_matrix.index = character_matrix.index.map(str)

    mutation_prior_dict = load_mutation_prior_dict(args)

    scoring_function = ScoringFunction(mutation_prior_dict, character_matrix)

    seed_parsimony = compute_parsimony(mutation_prior_dict, seed_tree, character_matrix) 
    logger.info(f"Seed tree parsimony: {seed_parsimony}")

    candidate_trees = []
    for i in range(1):
        candidate_tree = seed_tree.copy() 
        candidate_tree = stochastic_nni(candidate_tree, aggression=1.3)
        candidate_parsimony = scoring_function.score(candidate_tree) 
        candidate_trees.append((candidate_tree, candidate_parsimony))

    candidate_trees = sorted(candidate_trees, key=lambda x: x[1])[:5]

    count = 0
    while count < args.iterations:
        candidate_trees = sorted(candidate_trees, key=lambda x: x[1])

        for (i, (T, _)) in enumerate(candidate_trees):
            candidate_parsimony = compute_parsimony(mutation_prior_dict, T, character_matrix)
            logger.info(f"Candidate tree {i + 1} has parsimony {candidate_parsimony}")

        candidate_tree, _ = random.sample(candidate_trees, 1)[0]
        if len(candidate_trees) != 1:
            candidate_tree = stochastic_nni(candidate_tree)
        candidate_tree_optimized = hill_climb(candidate_tree, scoring_function)
        candidate_tree_optimized_parsimony = scoring_function.score(candidate_tree_optimized)

        if len(candidate_trees) < 5:
            candidate_trees.append((candidate_tree_optimized, candidate_tree_optimized_parsimony))
            continue

        if candidate_trees[0][1] <= candidate_tree_optimized_parsimony:
            if candidate_tree_optimized_parsimony < candidate_trees[-1][1]:
                candidate_trees = candidate_trees[:-1] + [(candidate_tree_optimized, candidate_tree_optimized_parsimony)]

            logger.info(f"Did not update candidate trees.")
            count += 1
            continue

        count = 0
        logger.info(f"Updated candidate trees w/ parsimony: {candidate_tree_optimized_parsimony}")
        candidate_trees = [(candidate_tree_optimized, candidate_tree_optimized_parsimony)] + candidate_trees[:-1]

        best_tree_contracted = collapse_mutationless_edges(candidate_trees[0][0], character_matrix)
        newick_tree = tree_to_newick(best_tree_contracted)
        with open(args.output, 'w') as f:
            f.write(f"{newick_tree};")

if __name__ == "__main__":
    args = parse_args()

    if not args.m and not args.e:
        print("Please specify either a mutation prior dictionary or CSV.")
        sys.exit(1)

    if args.mode == 'collapse':
        collapse_mode(args)
    elif args.mode == 'infer':
        infer_mode(args)
    elif args.mode == 'score':
        score_mode(args)
