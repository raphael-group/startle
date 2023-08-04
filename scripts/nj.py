from skbio import DistanceMatrix, TreeNode
from skbio.tree import nj

import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Runs neighbor joining on a character-state matrix using the Hamming distance."
    )

    parser.add_argument(
        "character_matrix", help="Distance matrix to do NJ on"
    )

    parser.add_argument(
        "--output", help="Output tree.", required=True
    )

    return parser.parse_args()

def hamming_distance(x, y):
    return np.sum(x != y)

if __name__ == "__main__":
    args = parse_arguments()
    
    character_matrix = pd.read_csv(args.character_matrix, index_col=0)
    pairwise_distances = np.zeros((len(character_matrix), len(character_matrix)))
    for i, row in enumerate(character_matrix.iterrows()):
        for j, row2 in enumerate(character_matrix.iterrows()):
            pairwise_distances[i][j] = hamming_distance(row[1], row2[1])

    names = character_matrix.index.tolist()
    dm = DistanceMatrix(pairwise_distances, list(map(lambda n: str(n).replace(" ", "_"), names)))

    tree = nj(dm)

    tree.write(f"{args.output}")
