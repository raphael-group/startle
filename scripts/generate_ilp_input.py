#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import math
import itertools
import json

from collections import deque

def get_binary_matrix(df_character_matrix):
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

    pruned_mutations = []
    pruned_cells = []
    while False:
        delta_mutations = []
        for mutation in df_binary.columns:
            ncounts_one = np.sum(df_binary[mutation] == 1)
            ncounts_zero = np.sum(df_binary[mutation] == 0)
            if ncounts_one <= 1 or ncounts_zero == 0:
                delta_mutations.append(mutation)

        if len(delta_mutations) == 0:
            break

        pruned_mutations += delta_mutations

        df_binary = df_binary.drop(delta_mutations, axis = 1)

        print(f'removed {len(delta_mutations)} mutations')
      
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
            df_binary = df_binary.drop(cell, axis=0)

    return df_binary, pruned_mutations, eq_class_dict

def main(args):
    df_character_matrix = pd.read_csv(args.d, index_col = 0)
    df_binary, pruned_mutations, eq_class_dict = get_binary_matrix(df_character_matrix)

    with open(f'{args.o}_equivalence_classes.json', 'w') as out:
        out.write(json.dumps(eq_class_dict))

    with open(f'{args.o}_pruned_mutations.txt', 'w') as out:
        for mutation in pruned_mutations:
            out.write(f'{mutation}\n')

    df_count = pd.read_csv(args.c, index_col = 2)
    max_allowed_homoplasy = df_count['count']
    
    cell_list = list(df_binary.index)
    mutation_list = list(df_binary.columns)
    mutation_to_index = {x: idx for idx, x in enumerate(mutation_list)}
    
    with open(f'{args.o}_binary_character_matrix.txt', 'w') as out:
        for cell in cell_list:
            cell_entries = []
            for mutation in mutation_list:
                cell_entries.append(df_binary.loc[cell][mutation])
            out.write(' '.join(map(str, cell_entries)))
            out.write('\n')
    
    with open(f'{args.o}_counts.txt', 'w') as out:
        for mutation in mutation_list:
            out.write(f'{max_allowed_homoplasy[mutation]}\n')
        
    with open(f'{args.o}_mutations.txt', 'w') as out:
        for mutation in mutation_list:
            out.write(f'{mutation}\n')
            
    with open(f'{args.o}_cells.txt', 'w') as out:
        for cell in cell_list:
            out.write(f'{cell}\n')
        
    character_list = list(set([x.split('_')[0] for x in df_binary.columns]))
            
    one_cell_mut_list = []
    var_index_list = []
    for cell_idx, cell in enumerate(cell_list):
        for mut_idx, mut in enumerate(mutation_list):
            if df_binary.loc[cell][mut] == 1:
                one_cell_mut_list.append((cell_idx, mut_idx))

    missing_cell_character_list = []
    for character_idx, character in enumerate(character_list):
        for cell_idx, cell in enumerate(cell_list):
            if df_character_matrix.loc[cell][character] == -1:
                missing_cell_character_list.append((cell_idx, character_idx))            
                
    with open(f'{args.o}_one_indices.txt', 'w') as out:
        for cell_idx, mut_idx in one_cell_mut_list:
            out.write(f'{cell_idx} {mut_idx}\n')
    
    with open(f'{args.o}_missing_indices.txt', 'w') as out:
        for cell_idx, character_idx in missing_cell_character_list:
            out.write(f'{cell_idx} {character_idx}\n')
    
    with open(f'{args.o}_character_mutation_mapping.txt', 'w') as out:
        for _, character in enumerate(character_list):
            character_mutation_list = [mutation_to_index[x] for x in mutation_list if x.startswith(f'{character}_')]
            out.write(' '.join(map(str, character_mutation_list)) + '\n')
    
    df_weights = pd.read_csv(args.w)
    df_weights['mutation'] = df_weights.apply(lambda x: f"{x['character']}_{x['state']}", axis = 1)
    if 'probability' in df_weights.columns:
        df_weights['weight'] = -np.log(df_weights['probability'])
    else:
        df_weights['weight'] = -np.log(df_weights['freq'])
    df_weights = df_weights.set_index('mutation')    
    
    with open(f'{args.o}_weights.txt', 'w') as out:
        for mutation in mutation_list:
            out.write(f"{df_weights.loc[mutation]['weight']}\n")    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='character-state matrix')
    parser.add_argument('-c', type=str, help='count matrix')
    parser.add_argument('-w', type=str, help='weights')
    parser.add_argument('-o', type=str, help='output prefix', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
