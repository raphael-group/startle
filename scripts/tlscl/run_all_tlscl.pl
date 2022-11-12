#!/bin/perl

use v5.10;

use strict;
use warnings;

for my $i (1..24) {
    `sbatch -e results/bar${i}.err.out -o results/bar${i}.std.out --export=experiment=$i scripts/tlscl/sbatch.sh`;
}
