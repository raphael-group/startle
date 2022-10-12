#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 2022

@author: Palash Sashittal
"""

import pandas as pd
import seaborn as sns
import numpy as np
import math
import itertools
from collections import Counter
from collections import defaultdict
from Bio import Phylo
import networkx as nx
# import cassiopeia as cas
# from cassiopeia.critique import critique_utilities

def get_binary_dataframe(df, root_state = '0', missing_state='-1'):
    character_state_dict = {c:list(df[c].unique()) for c in df.columns}

    bin_char_state_list = []
    bin_char_list = []
    for character, state_list in character_state_dict.items():
        for state in state_list:
            if state != root_state and state != missing_state:
                binary_array = np.array(df[character] == state, dtype=int)
                binary_array[np.where(df[character] == missing_state)] = -1
                bin_char_state_list.append(binary_array)
                bin_char_list.append(f'{character}_{state}')

    return pd.DataFrame(np.vstack(bin_char_state_list).T, columns=bin_char_list, index=df.index)

def from_newick_get_castree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    # new_net_tree = net_tree.copy()
    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
    #     print(node)    
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'

    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))

    return cas.data.CassiopeiaTree(tree=directed_tree)


def get_correct_resolved_triplets(T1, T2, number_of_trials = 100, min_triplets_at_depth = 1):

    all_triplets_correct = defaultdict(int)
    unresolved_triplets_correct = defaultdict(int)
    resolvable_triplets_correct = defaultdict(int)
    proportion_unresolvable = defaultdict(int)
    
    # create copies of the trees and collapse process
#     T1 = copy.copy(tree1)
#     T2 = copy.copy(tree2)

    T1.collapse_unifurcations()
    T2.collapse_unifurcations()

    # set depths in T1 and compute number of triplets that are rooted at
    # ancestors at each depth
    depth_to_nodes = critique_utilities.annotate_tree_depths(T1)

    max_depth = np.max([T1.get_attribute(n, "depth") for n in T1.nodes])
    for depth in range(max_depth):

        score = 0
        number_unresolvable_triplets = 0

        # check that there are enough triplets at this depth
        candidate_nodes = depth_to_nodes[depth]
        total_triplets = sum([T1.get_attribute(v, "number_of_triplets") for v in candidate_nodes])
        if total_triplets < min_triplets_at_depth:
            continue

        for _ in range(number_of_trials):

            (i, j, k), out_group = critique_utilities.sample_triplet_at_depth(T1, depth, depth_to_nodes)

            reconstructed_outgroup = critique_utilities.get_outgroup(T2, (i, j, k))
            
            is_resolvable = True
            if reconstructed_outgroup == 'None':
#             if out_group == "None":
                number_unresolvable_triplets += 1
                is_resolvable = False

            # increment score if the reconstructed outgroup is the same as the
            # ground truth
            score = int(reconstructed_outgroup == out_group)

            all_triplets_correct[depth] += score
            if is_resolvable:
                resolvable_triplets_correct[depth] += score
            else:
                unresolved_triplets_correct[depth] += score

        all_triplets_correct[depth] /= number_of_trials

        if number_unresolvable_triplets == 0:
            unresolved_triplets_correct[depth] = 1.0
        else:
            unresolved_triplets_correct[depth] /= number_unresolvable_triplets

        proportion_unresolvable[depth] = (
            number_unresolvable_triplets / number_of_trials
        )

        if proportion_unresolvable[depth] < 1:
            resolvable_triplets_correct[depth] /= (
                number_of_trials - number_unresolvable_triplets
            )
        else:
            resolvable_triplets_correct[depth] = 1.0

    return (
        all_triplets_correct,
        resolvable_triplets_correct,
#         unresolved_triplets_correct,
        proportion_unresolvable,
    )

def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    # new_net_tree = net_tree.copy()
    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = str(node)
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'
    
    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))

    return directed_tree

def is_leaf(T, node):
    if len(T[node]) == 0:
        return True
    else:
        return False

def update_labeling_count(T, character_series, labeling, count_dict, node, root_state = '0', missing_state = '-1'):
    if is_leaf(T, node):
        labeling[node] = character_series[node]
    else:
        child_states = []
        for child in T[node]:
            update_labeling_count(T, character_series, labeling, count_dict, child)
            child_states.append(labeling[child])
        
        if len(set(child_states) - {missing_state}) > 1:
            labeling[node] = root_state
            for state in child_states:
                if state != root_state and state != missing_state:
                    count_dict[state] += 1
        else:
            if len(set(child_states) - {missing_state}) == 0:
                labeling[node] = missing_state
            else:
                labeling[node] = list(set(child_states) - {missing_state})[0]

def compute_labeling_count(T, df):
    count_dict = {x:Counter() for x in df}
    labeling_dicts = {x:{node:'' for node in T.nodes} for x in df}
    for character in df:
        update_labeling_count(T, df[character], labeling_dicts[character], count_dict[character], 'root')
    
    return labeling_dicts, count_dict

def compute_likelihood(df, count_dict):
    likelihood = 0
    for character, state_counter in count_dict.items():
        for state, count in state_counter.items():
            likelihood += count * (-math.log(float(df.loc[character, state]['probability'])))
    return likelihood

def tree_to_newick(T, root=None):
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
            subgs.append(tree_to_newick(T, root=child))
        else:
            subgs.append(child)
    return "(" + ','.join(map(str, subgs)) + ")"

