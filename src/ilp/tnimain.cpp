/*
 * tnimain.cpp
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#include "utils.h"
#include "ilpsolver.h"
#include "matrix.h"
#include <fstream>
#include <lemon/arg_parser.h>

//#ifdef CPLEX
//#include "ilpsolvercplex.h"
//// #else
//// #include "ilpsolvergurobi.h"
//#endif

int main(int argc, char** argv)
{
	std::string one_indices_filename, missing_indices_filename, count_filename, character_mutation_filename, weight_filename, sol_filename;
		
	int ncells, nmutations;
	int nrThreads = 1;
    int timeLimit = -1;
    int memLimit = -1;
    double EPgap = 0.01;
    bool allViolations = false;
    
	lemon::ArgParser ap(argc, argv);
	ap.refOption("n", "number of cells", ncells);
	ap.refOption("m", "number of mutations", nmutations);
    ap.refOption("t", "number of threads (default: 1)", nrThreads);
    ap.refOption("T", "time limit in seconds (default: unlimited)", timeLimit);
    ap.refOption("M", "memory limit in MB (default: unlimited)", memLimit);
    ap.refOption("g", "gap in optimality (default: 0.01)", EPgap);
    ap.refOption("f", "forces all violations to be created at beginning (default: false)", allViolations);

	ap.other("<one indices file>");
	ap.other("<missing indices file>");
	ap.other("<count file>");
	ap.other("<character-mutation mapping>");
	ap.other("<weight file>");
	ap.other("<output file>");
	ap.parse();

    std::cerr << allViolations << std::endl;
		
    if (ap.files().size() != 6)
    {
        std::cerr << "Error1: expected <one indices file> <missing indices file> <count file> <character mutation file> <weight file> <output file>" << std::endl;
        return 1;
    }

	one_indices_filename = ap.files()[0];
	missing_indices_filename = ap.files()[1];
	count_filename = ap.files()[2];
	character_mutation_filename = ap.files()[3];
	weight_filename = ap.files()[4];
	sol_filename = ap.files()[5];
    
	std::ifstream one_indices_file(one_indices_filename.c_str());
	if (!one_indices_file.good())
	{
			std::cerr << "Error2: failed opening '" << one_indices_filename << "' for reading" << std::endl;
			return 1;
	}
	std::ifstream missing_indices_file(missing_indices_filename.c_str());
	if (!missing_indices_file.good())
	{
			std::cerr << "Error3: failed opening '" << missing_indices_filename << "' for reading" << std::endl;
			return 1;
	}
	std::ifstream count_file(count_filename.c_str());
	if (!count_file.good())
	{
			std::cerr << "Error4: failed opening '" << count_filename << "' for reading" << std::endl;
			return 1;
	}
	std::ifstream character_mutation_file(character_mutation_filename.c_str());
	if (!character_mutation_file.good())
	{
			std::cerr << "Error5: failed opening '" << character_mutation_filename << "' for reading" << std::endl;
			return 1;
	}
	std::ifstream weight_file(weight_filename.c_str());
	if (!weight_file.good())
	{
			std::cerr << "Error6: failed opening '" << weight_filename << "' for reading" << std::endl;
			return 1;
	}
  
	Matrix input(ncells, nmutations);
	
	// read one index file
	input.readOneIndices(one_indices_file);
	one_indices_file.close();
	// read missing index file
	input.readMissingIndices(missing_indices_file);
	missing_indices_file.close();
	// read count file
	input.readMaxAllowedHomoplasy(count_file);
	count_file.close();
	// read character mutation mapping file
	input.readCharacterMutationMapping(character_mutation_file);
	character_mutation_file.close();
	// read weights
	input.readWeights(weight_file);
	weight_file.close();
	// initialize free indices
	input.initFreeIndices();
	
	IlpSolver solver(input);
	solver.init();
	
	std::ofstream output(sol_filename.c_str());
	if (solver.solve(timeLimit, memLimit, nrThreads, EPgap, true, allViolations))
	{		
		solver.writeSolution(output);
	}
	else
	{
		output << "no solution";
	}
	output.close();
	
	return 0;
}
