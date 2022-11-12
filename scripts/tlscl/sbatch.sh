#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 32G
#SBATCH -t 1-00:00:00 

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=hs2435@cs.princeton.edu

PYTHON=/n/fs/ragr-data/users/schmidt/miniconda3/envs/lineage-tracing/bin/python

character_matrix=../lineage-tracing/tlscl-data/bar"$experiment"_character_matrix.csv
mutation_prior=../lineage-tracing/tlscl-data/bar"$experiment"_mutation_prior.csv

export HOME=/n/fs/grad/hs2435
export PATH=/u/hs2435/.local/bin:/n/fs/ragr-research/users/schmidth/julia-1.7.3/bin:/n/fs/ragr-data/users/schmidt/texlive/2022/bin/x86_64-linux:/n/fs/ragr-data/users/schmidt/miniconda3/envs/lineage-tracing/bin:/n/fs/ragr-data/users/schmidt/miniconda3/condabin:/u/hs2435/bin:/usr/share/Modules/bin:/u/hs2435/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/puppetlabs/bin:/u/hs2435/bin:/n/fs/ragr-research/projects/network-mutations/08-25-2020/gurobi903/linux64/bin:/u/hs2435/.fzf/bin

$PYTHON scripts/startle_tp.py $character_matrix $mutation_prior --output results/bar$experiment
