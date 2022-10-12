# startle

## start-ILP

Requirements:
    - [LEMON](https://lemon.cs.elte.hu/trac/lemon) Graph Library
    - IBM ILOG [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) Optimizer

build command:
```
    $ cmake -DLIBLEMON_ROOT=/n/fs/ragr-data/users/palash/libraries/lemon -DCPLEX_INC_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/cplex/include -DCPLEX_LIB_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic -DCONCERT_INC_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/concert/include -DCONCERT_LIB_DIR=/n/fs/ragr-code/general/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic ..
```
