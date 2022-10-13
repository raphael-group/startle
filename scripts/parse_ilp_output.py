import argparse
import sys
import pandas as pd
import seaborn as sns
import numpy as np
import networkx as nx
import math
import itertools
import json

from collections import Counter
from collections import defaultdict
from Bio import Phylo

def generate_perfect_phylogeny(df_binary):
    solT_mut = nx.DiGraph()
    solT_mut.add_node('root')

    solT_cell = nx.DiGraph()
    solT_cell.add_node('root')

    df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

    for cell_id, row in df_binary.iterrows():
        if cell_id == 'root':
            continue

        curr_node = 'root'
        for column in df_binary.columns[row.values == 1]:
            if column in solT_mut[curr_node]:
                curr_node = column
            else:
                if column in solT_mut.nodes:
                    raise NameError(f'{column} is being repeated')
                solT_mut.add_edge(curr_node, column)
                solT_cell.add_edge(curr_node, column)
                curr_node = column

        solT_cell.add_edge(curr_node, cell_id)   

    return solT_mut, solT_cell

def tree_to_newick(T, eq_class_dict, root=None):
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
            subgs.append(tree_to_newick(T, eq_class_dict, root=child))
        else:
            eq_class = eq_class_dict[child]
            if len(eq_class) > 1:
                subgs.append("(" + ','.join(map(str, eq_class)) + ")")
            else:
                subgs.append(eq_class[0])
    return "(" + ','.join(map(str, subgs)) + ")"

def main(args):
    df_cpp_output = pd.read_csv(args.i, header=None, sep=' ')
    df_cpp_output = df_cpp_output.rename(columns={0:'cell_idx', 1:'mut_idx', 2:'state_idx', 3:'entry'})

    cell_list = []
    mutation_list = []
    with open(args.c, 'r') as inp:
        for line in inp:
            cell_name = line.rstrip('\n')
            if len(cell_name) > 0:
                cell_list.append(cell_name)

    with open(args.m, 'r') as inp:
        for line in inp:
            mutation_name = line.rstrip('\n')
            if len(mutation_name) > 0:
                mutation_list.append(mutation_name)
                
    with open(args.e, 'r') as inp:
        eq_class_dict = json.load(inp)
    
    ncells = len(cell_list)
    nmutations = len(mutation_list)
    
    df_cpp_output['name'] = df_cpp_output.apply(lambda x: f"{mutation_list[x['mut_idx']]}_{x['state_idx']}", axis =1)
    sol_columns = list(df_cpp_output['name'].unique())
    nsol_columns = len(sol_columns)
    
    sol_entries = np.zeros((ncells, nsol_columns), dtype=int)
    for mut_idx, mut in enumerate(sol_columns):
        for cell_idx in df_cpp_output[(df_cpp_output['entry'] == 1) & (df_cpp_output['name'] == mut)]['cell_idx']:
            sol_entries[cell_idx][mut_idx] = 1    

    df_sol_binary = pd.DataFrame(sol_entries, columns=sol_columns, index=cell_list)
    df_sol_binary.to_csv(f'{args.o}_completed_binary_character_matrix.csv')
    
    solT_mut, solT_cell = generate_perfect_phylogeny(df_sol_binary)

    node_list = list(solT_mut.nodes)
    node_to_index = {node:node_idx for node_idx, node in enumerate(node_list)}
    with open(f'{args.o}_edge_labels.csv', 'w') as out:
        out.write('src,dst,character,state\n')
        for edge in solT_mut.edges:
            src_node = edge[0]
            dst_node = edge[1]

            src_node_idx = node_to_index[src_node]
            dst_node_idx = node_to_index[dst_node]

            edge_character = dst_node.split('_')[0]
            edge_state = dst_node.split('_')[1]

            out.write(f'{src_node_idx},{dst_node_idx},{edge_character},{edge_state}\n')

    with open(f'{args.o}_tree.newick', 'w') as out:
        out.write(f"{tree_to_newick(solT_cell, eq_class_dict)};")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', type=str, help='cpp code output')
    parser.add_argument('-e', type=str, help='equivalence classes')
    parser.add_argument('-m', type=str, help='list of mutations')
    parser.add_argument('-c', type=str, help='list of cells')
    parser.add_argument('-o', type=str, help='output prefix', required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
