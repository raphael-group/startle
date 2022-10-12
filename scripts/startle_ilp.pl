#!/bin/perl

use v5.10;

use strict;
use warnings;

use Text::CSV qw( csv );
use Getopt::ArgParse;

my $STARTLE = "./build/startle";
 
my $ap = Getopt::ArgParse->new_parser(
    prog => 'Infer lineage trees using Startle-ILP solver.',
);

$ap->add_args(
    ['--character-matrix', '-c', required => 1, help => 'Input character matrix as CSV.'],
    ['--mutation-priors', '-m', required => 1, help => 'Input mutation prior probabilities as CSV.'],
    ['--output', '-o', required => 1, help => 'Output directory name.']
);

my $ns = $ap->parse_args();

my $character_matrix = $ns->character_matrix;
my $mutation_priors = $ns->mutation_priors;
my $output = $ns->output;

`mkdir -p ${output}`;

my $CPP_INPUTS = "${output}/cpp_inputs";

# SEED COUNTS AND GENERATE ILP INPUT
my $counts = "${CPP_INPUTS}_counts.csv";
`python scripts/seed_counts.py -d ${character_matrix} -p ${mutation_priors} -s greedy > ${counts}`;
`python scripts/generate_ilp_input.py -d ${character_matrix} -c ${counts} -w ${mutation_priors} -o ${CPP_INPUTS}`;

# RUN ILP
$character_matrix = "${CPP_INPUTS}_binary_character_matrix.txt";
$counts = "${CPP_INPUTS}_counts.txt";

my $one_indices = "${CPP_INPUTS}_one_indices.txt";
my $missing_indices = "${CPP_INPUTS}_missing_indices.txt";
my $char_mut_mapping = "${CPP_INPUTS}_character_mutation_mapping.txt";
my $weights = "${CPP_INPUTS}_weights.txt";
my $results = "${output}/ilp_results.txt";

my $aoa = csv (in => $character_matrix, sep_char => " ");  
my $n = scalar @$aoa;
my $m = scalar @{$aoa->[0]};

my @startle_args = (
    $STARTLE, "-n", "$n", "-m", "$m",
    "-t", "48", "-T", "7200", "-g", "0.05",
    $one_indices, $missing_indices, $counts,
    $char_mut_mapping, $weights, $results,
);

my $cmd = join (' ' , @startle_args);
`$cmd`;

# PARSE ILP OUTPUT
my @parse_ilp_output_args = (
    "python", "scripts/parse_ilp_output.py",
    "-i", $results,
    "-m", "${CPP_INPUTS}_mutations.txt",
    "-c", "${CPP_INPUTS}_cells.txt",
    "-e", "${CPP_INPUTS}_equivalence_classes.json",
    "-o", "${output}/startle"
);

$cmd = join (' ', @parse_ilp_output_args);
`$cmd`;
