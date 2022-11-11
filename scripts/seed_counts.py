import pandas as pd
import sys
import argparse
import numpy as np
import random

import networkx as nx
import cassiopeia as cas

from collections import defaultdict

def get_character_taxon_dataframe(fname):
    cell_list = []
    state_list = []
    with open(fname) as inp:
        flag = 0
        for line in inp:
            if flag == 0:
                flag = 1
                continue
            data = line.rstrip('\n')
            cellname, state = data.split('\t')
            cell_list.append(cellname)
            state_list.append(np.array([int(x) for x in state]))
    
    return pd.DataFrame(np.vstack(state_list), index=cell_list, columns = [f'c{i}' for i in range(len(state_list[0]))])

def display_counts(character_matrix, tree):
    tree.collapse_mutationless_edges(True)

    count_dict = defaultdict(dict)
    for (u, v) in tree.edges:
        for (c, s) in tree.get_mutations_along_edge(u, v):
            if s in count_dict[c]:
                count_dict[c][s] += 1
            else:
                count_dict[c][s] = 1

    root_states = tree.get_character_states(tree.root)
    for (c, s) in enumerate(root_states):
        if s in count_dict[c]:
            count_dict[c][s] += 1
        else:
            count_dict[c][s] = 1

    characters = list(character_matrix.columns)

    print('character,state,mutation,count')
    for idx, count_data in count_dict.items():
        for state, count in count_data.items():
            if state == 0: continue
            print(f'{characters[idx]},{state},{characters[idx]}_{state},{count}')

def main(args):
    random.seed(0)
    np.random.seed(0)
    
    if args.f:
        df = get_character_taxon_dataframe(args.f)
    elif args.d:
        df = pd.read_csv(args.d, index_col = 0)
    
    df.index = df.index.astype(str, copy = False)
    
    if args.p:
        df_mutation_prior = pd.read_csv(args.p)
        df_mutation_prior['state'] = df_mutation_prior['state'].map(str)
        df_mutation_prior = df_mutation_prior.set_index(['character', 'state'])
        character_state_prior = {}
        for idx, character in enumerate(list(df.columns)):
            state_list = list(df_mutation_prior.loc[character].index)
            character_state_prior[idx] = {}
            for state in state_list:
                character_state_prior[idx][int(state)] = df_mutation_prior.loc[character, state]['probability']
    else:
        character_state_prior = None

    reconstructed_tree = cas.data.CassiopeiaTree(character_matrix = df,
                                                 missing_state_indicator = int(args.m),
                                                 priors = character_state_prior)
    
    solver_dict = {'greedy': cas.solver.VanillaGreedySolver(),
                   'upgma': cas.solver.UPGMASolver(dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance),
                   'nj': cas.solver.NeighborJoiningSolver(dissimilarity_function=cas.solver.dissimilarity.weighted_hamming_distance, add_root=True),
                   'ILP': cas.solver.ILPSolver(convergence_time_limit=3600),
                   'hybrid': cas.solver.HybridSolver(cas.solver.VanillaGreedySolver(), cas.solver.ILPSolver(), cell_cutoff=35)}
    
    solver = solver_dict[args.s]
    solver.solve(reconstructed_tree)

    display_counts(df, reconstructed_tree)

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=str, help='solver name [greedy, upgma, nj, maxcut]', default='greedy')
    parser.add_argument('-f', type=str, help='cell-state file in tsv format')
    parser.add_argument('-d', type=str, help='csv file with character-taxon dataframe')
    parser.add_argument('-m', type=int, help='missing data state for the characters [-1]', default=-1)
    parser.add_argument('-p', type=str, help='csv file with probability distribution of occurence of the character-state pairs [uniform]')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    if not args.d and not args.f:
        raise Exception('need either cell-state file or character-taxon dataframe as input')
    elif args.s not in ['greedy', 'upgma', 'nj', 'maxcut', 'ILP', 'hybrid']:
            raise Exception('solver must be one of the following -- greedy, upgma, nj, maxcut, ILP, hybrid')
    else:
        main(args)
