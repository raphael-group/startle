# startle

startle is a method for lineage tree reconstruction that uses the
*star homoplasy* evolutionary model. 

If you found the tool useful in your research, please cite us at:
```
```

## startle-ILP

startle-ILP is a command-line tool for inferring lineage trees from 
cell-by-state matrices. It takes as input a character-matrix and
a set of estimated mutation probabilities for each character-state pair
and outputs a inferred tree in Newick format.

### Installation

To install startle-ILP you will need the following set of libraries
and programming languages available on your device.

Requirements:
- [LEMON](https://lemon.cs.elte.hu/trac/lemon) Graph Library
- IBM ILOG [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) Optimizer
- C++11 Compatible Compiler
- [perl](https://www.perl.org/) v5.10 or higher
- [python3](https://www.python.org/downloads/)

After ensuring the aforementioned libraries are properly installed,
run the following command to build the `startle` integer linear programming
solver - the output file will be named `startle` and will be located under 
the `build/` directory.

```
$ mkdir build; cd build
$ cmake -DLIBLEMON_ROOT=/n/fs/ragr-data/users/palash/libraries/lemon\
        -DCPLEX_INC_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/cplex/include\
        -DCPLEX_LIB_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic\
        -DCONCERT_INC_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/concert/include\
        -DCONCERT_LIB_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic\
        ..
$ make
```

To install the required python dependencies, we recommend
you create a new conda environment for this project. 

```
$ conda create -n startle
$ conda install python=3.8
$ conda install -c conda-forge biopython
$ conda install -c conda-forge funcy loguru tqdm networkx pandas numpy seaborn
```

### Running

```
$ perl scripts/startle-ilp.pl --help
usage: Infer lineage trees using Startle-ILP solver. --output|-o
--mutation-priors|-m --character-matrix|-c

required named arguments:
  --output, -o OUTPUT                          Output directory name.
  --mutation-priors, -m MUTATION-PRIORS        Input mutation prior probabilities as CSV.
  --character-matrix, -c CHARACTER-MATRIX      Input character matrix as CSV.
```
