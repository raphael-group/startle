/*
 * ilpsolver.cpp
 *
 *  Created on: 1-aug-2022
 *      Author: P. Sashittal
 */

#include "ilpsolver.h"

IlpSolver::IlpSolver(const Matrix& input)
: _input(input)
, _ncells(input.getNumCells())
, _nmutations(input.getNumMutations())
, _ncharacters(input.getNumCharacters())
, _violations()
, _violatedMutations()
, _nActivated()
, _nConstraints()
, _nviolations()
, _homoplasyVector()
, _env()
, _model(_env)
, _cplex(_model)
, _obj(_env)
, _x(_env)
, _y1(_env)
, _y2(_env)
, _z(_env)
{
	std::cout << "number of cells is " << _ncells << std::endl;
	std::cout << "number of mutations is " << _nmutations << std::endl;
	std::cout << "number of characters is " << _ncharacters << std::endl;
}

void IlpSolver::init()
{
//	IntPairVector one_indices = _input.getOneIndices();
//	std::cout << "size of one indices is: " << one_indices.size() << std::endl;
    initVariables();
    initConstraints();
	initObjective();
	
	std::cout << "finished initialization" << std::endl;
	std::cout << "number of activated variables: " << _nActivated << std::endl;
}

void IlpSolver::initVariables()
{
	char buf[1024];
	
	// initialize variables
	
	// x
	_x = IloBoolVar3Matrix(_env, _ncells);
//	_activeVariables = Bool3Matrix(_ncells);
	for (int i = 0; i < _ncells; ++i)
	{
		_x[i] = IloBoolVarMatrix(_env, _nmutations);
//		_activeVariables[i] = BoolMatrix(_nmutations);
		for (int j = 0; j < _nmutations; ++j)
		{
			int nstates = _input.getMaxCharacterHomoplasy(j);
			_x[i][j] = IloBoolVarArray(_env, nstates);
//			_activeVariables[i][j] = BoolVector(nstates);
			for (int s = 0; s < nstates; ++s)
			{
				snprintf(buf, 1024, "x_%d_%d_%d", i, j, s);
				_x[i][j][s] = IloBoolVar(_env, buf);
//				_activeVariables[i][j][s] = false;
			}
		}
	}
	
	// w
	_w = IloNumVarMatrix(_env, _nmutations);
	for (int j = 0; j < _nmutations; ++j)
	{
		int nstates = _input.getMaxCharacterHomoplasy(j);
		_w[j] = IloNumVarArray(_env, nstates);
		for (int s = 0; s < nstates; ++s)
		{
			snprintf(buf, 1024, "w_%d_%d", j, s);
			_w[j][s] = IloNumVar(_env, 0, 1, buf);
		}
	}
		
	// init num of conflict constraints
	_nConflictConstraints = 0;
	
	// init homoplasy vector
	_homoplasyVector = IntVector(_nmutations, 1);

//	// activate 1 variables
//	for (IntPair cell_mut : _input.getOneIndices())
//	{
//		int cell_idx = cell_mut.first;
//		int mut_idx = cell_mut.second;
//		_activeVariables[cell_idx][mut_idx][0] = true;
//	}
//
//	// activate missing variables
//	for (IntPair cell_character : _input.getMissingIndices())
//	{
//		int cell_idx = cell_character.first;
//		int char_idx = cell_character.second;
//		for (int mut_idx : _input.getCharacterMutations(char_idx))
//		{
//			_activeVariables[cell_idx][mut_idx][0] = true;
//		}
//	}
	
	// update variable activation
//	updateVariableActivation();
	
	// update variable bounds for x
	updateVariableBounds();
}

//void IlpSolver::updateVariableActivation()
//{
//	_nActivated = 0;
//
//	// activate 1 variables
//	for (IntPair cell_mut: _input.getOneIndices())
//	{
//		int cell_idx = cell_mut.first;
//		int mut_idx = cell_mut.second;
//		for (int s = 0; s < _homoplasyVector[mut_idx]; ++s)
//		{
//			_activeVariables[cell_idx][mut_idx][s] = true;
//			_nActivated++;
//		}
//	}
//
//	// activate missing variables
//	for (IntPair cell_character : _input.getMissingIndices())
//	{
//		int cell_idx = cell_character.first;
//		int char_idx = cell_character.second;
//		for (int mut_idx : _input.getCharacterMutations(char_idx))
//		{
//			for (int s = 0; s < _homoplasyVector[mut_idx]; ++s)
//			{
//				_activeVariables[cell_idx][mut_idx][0] = true;
//				_nActivated++;
//			}
//		}
//	}
//}

void IlpSolver::updateVariableBounds()
{
	_nActivated = 0;

	for (IntPair cell_mut : _input.getFreeIndices())
	{
		int cell_idx = cell_mut.first;
		int mut_idx = cell_mut.second;
		int nstates = _input.getMaxCharacterHomoplasy(mut_idx);
		
//		std::cout << "here1 " << mut_idx << std::endl;
		
		for (int s = 0; s < nstates; ++s)
		{
			if (s < _homoplasyVector[mut_idx])
			{
//				std::cout << "here2 " << mut_idx << std::endl;
				_x[cell_idx][mut_idx][s].setUB(1);
				_nActivated++;
			}
			else
			{
				_x[cell_idx][mut_idx][s].setUB(0);
			}
		}
	}

//	std::cout << "number of active variables: " << _nActivated << std::endl;
//	for (int i = 0; i < _ncells; ++i)
//	{
//		for (int j = 0; j < _nmutations; ++j)
//		{
//			int nstates = _input.getMaxCharacterHomoplasy(j);
//			for (int s = 0; s < nstates; ++s)
//			{
//				if (_activeVariables[i][j][s])
//				{
//					_x[i][j][s].setUB(1);
//				}
//				else
//				{
//					_x[i][j][s].setUB(0);
//				}
//			}
//		}
//	}
	
	for (int j = 0; j < _nmutations; ++j)
	{
		int nstates = _input.getMaxCharacterHomoplasy(j);
		for (int s = 0; s < nstates; ++s)
		{
			if (s < _homoplasyVector[j])
			{
				_w[j][s].setUB(1);
			}
			else
			{
				_w[j][s].setUB(0);
			}
		}
	}
}

void IlpSolver::initConstraints()
{
	IloExpr xsum(_env);
	
	// 1 variable constraints
	for (IntPair cell_mut : _input.getOneIndices())
	{
		int cell_idx = cell_mut.first;
		int mut_idx = cell_mut.second;
		int nstates = _input.getMaxCharacterHomoplasy(mut_idx);
		for (int s = 0; s < nstates; s++)
		{
			xsum += _x[cell_idx][mut_idx][s];
		}
		_model.add(xsum == 1);
		xsum.clear();
	}
	
	// missing vairable constraints
	for (IntPair cell_character : _input.getMissingIndices())
	{
		int cell_idx = cell_character.first;
		int char_idx = cell_character.second;
		for (int mut_idx : _input.getCharacterMutations(char_idx))
		{
			int nstates = _input.getMaxCharacterHomoplasy(mut_idx);
			for (int s = 0; s < nstates; s++)
			{
				xsum += _x[cell_idx][mut_idx][s];
			}
		}
		_model.add(xsum <= 1);
		xsum.clear();
	}
	
	// zero indices
	for (int i = 0; i < _ncells; ++i)
	{
		for (int j = 0; j < _nmutations; ++j)
		{
			if (!_input.isFree(i, j))
			{
				int nstates = _input.getMaxCharacterHomoplasy(j);
				for (int s = 0; s < nstates; ++s)
				{
					_model.add(_x[i][j][s] == 0);
				}
			}
		}
	}
	
	// symmetry breaking
	IloExpr prev_xsum(_env);
	int maxHomoplasy = _input.getMaxHomoplasy();
	for (int s = 1; s < maxHomoplasy; ++s)
	{
		for (IntPair cell_mut : _input.getFreeIndices())
		{
			int cell_idx = cell_mut.first;
			int mut_idx = cell_mut.second;
			int nstates = _input.getMaxCharacterHomoplasy(mut_idx);
//			if (s - 1 < _homoplasyVector[mut_idx])
			if (s - 1 < nstates)
			{
				prev_xsum += _x[cell_idx][mut_idx][s-1];
			}
//			if (s  < _homoplasyVector[mut_idx])
			if (s < nstates)
			{
				xsum += _x[cell_idx][mut_idx][s];
			}
		}
		
		_model.add(prev_xsum >= xsum);
		prev_xsum.clear();
		xsum.clear();
	}
	
	prev_xsum.end();
	xsum.end();
	
	// w constraints
	for (IntPair cell_mut : _input.getFreeIndices())
	{
		int cell_idx = cell_mut.first;
		int mut_idx = cell_mut.second;
		int nstates = _input.getMaxCharacterHomoplasy(mut_idx);
		for (int s = 0; s < nstates; ++s)
		{
			_model.add(_w[mut_idx][s] >= _x[cell_idx][mut_idx][s]);
		}
	}
}

void IlpSolver::setAllConstraints()
{
	_violations.clear();
	_violatedMutations.clear();

	for (int c1 = 0; c1 < _ncharacters; c1++)
	{
		for (int c2 = 0; c2 < c1; c2++)
		{
			IntVector c1_mutations = _input.getCharacterMutations(c1);
			IntVector c2_mutations = _input.getCharacterMutations(c2);
			
			for (int j1 : c1_mutations)
			{
				for (int j2 : c2_mutations)
				{
					for (int s1 = 0; s1 < _homoplasyVector[j1]; s1++)
					{
						for (int s2 = 0; s2 < _homoplasyVector[j2]; s2++)
						{
                            IntPair mutation_state1 = std::make_pair(j1, s1);
                            IntPair mutation_state2 = std::make_pair(j2, s2);
                            _violations.push_back(IntQuad(mutation_state1, mutation_state2));
                            _violatedMutations.insert(mutation_state1);
                            _violatedMutations.insert(mutation_state2);
                            _nviolations++;
						}
					}
				}
			}
		}
	}
}


int IlpSolver::checkViolation()
{
	_violations.clear();
//	IntPairSet violated_mutations;
	_violatedMutations.clear();
	_nviolations = 0;
	// if the one-set of two mutations is not related by containment
	// or disjoint, then they are in violation
	for (int c1 = 0; c1 < _ncharacters; c1++)
	{
		for (int c2 = 0; c2 < c1; c2++)
		{
			IntVector c1_mutations = _input.getCharacterMutations(c1);
			IntVector c2_mutations = _input.getCharacterMutations(c2);
			
			for (int j1 : c1_mutations)
			{
				for (int j2 : c2_mutations)
				{
					for (int s1 = 0; s1 < _homoplasyVector[j1]; s1++)
					{
						for (int s2 = 0; s2 < _homoplasyVector[j2]; s2++)
						{
							// compare (j1,s1) with (j1,s2)
							// gamete1 indicates (0,1)
							// gamete2 indicates (1,0)
							// gamete3 indicates (1,1)
							bool gamete1 = false, gamete2 = false, gamete3 = false;
							bool violated = false;
							
							for (int cell_idx = 0; cell_idx < _ncells; cell_idx++)
							{
								int entry1 = _cplex.getValue(_x[cell_idx][j1][s1]);
								int entry2 = _cplex.getValue(_x[cell_idx][j2][s2]);
//								int entry1, entry2;
//
//								if (_input.isFree(cell_idx, j1))
//								{
//									entry1 = _cplex.getValue(_x[cell_idx][j1][s1]);
//								}
//								else
//								{
//									entry1 = 0;
//								}
//
//								if (_input.isFree(cell_idx, j2))
//								{
//									entry2 = _cplex.getValue(_x[cell_idx][j2][s2]);
//								}
//								else
//								{
//									entry2 = 0;
//								}
								
								if (entry2 > 0.5)
								{
									if (entry1 > 0.5)
									{
										gamete3 = true;
									}
									else
									{
										gamete1 = true;
									}
								}
								else
								{
									if (entry1 > 0.5)
									{
										gamete2 = true;
									}
								}
							}
							if (gamete1 && gamete2 && gamete3)
							{
								IntPair mutation_state1 = std::make_pair(j1, s1);
								IntPair mutation_state2 = std::make_pair(j2, s2);
								_violations.push_back(IntQuad(mutation_state1, mutation_state2));
//								violated_mutations.insert(mutation_state1);
//								violated_mutations.insert(mutation_state2);
								_violatedMutations.insert(mutation_state1);
								_violatedMutations.insert(mutation_state2);
                                _nviolations++;
							}
						}
					}
				}
			}
		}
	}
	
	std::cout << "size of _violations is " << _violations.size() << " and _nviolations is " << _nviolations << std::endl;
		
	return _nviolations;
}

void IlpSolver::updateHomoplasyVector()
{
	// update homoplasy vector for violated mutations
	int idx = 0, nExpanded = 0;
	for (IntPair violated_mutation : _violatedMutations)
	{
		int mutation_idx = violated_mutation.first;
		int state_idx  = violated_mutation.second;
		int curr_homoplasy = _homoplasyVector[mutation_idx];
		
		if ((state_idx == curr_homoplasy - 1) && (curr_homoplasy < _input.getMaxCharacterHomoplasy(mutation_idx)))
		{
			_homoplasyVector[mutation_idx]++;
			nExpanded++;
		}
		++idx;
	}
	
	std::cout << "expanded " << nExpanded << " mutations out of " << idx << " violated mutations" << std::endl;

}

void IlpSolver::updateVariables()
{
	char buf[1024];
	
	for (int idx = 0; idx < _nviolations; ++idx)
	{
		int j1 = _violations[idx].first.first;
		int s1 = _violations[idx].first.second;
		int j2 = _violations[idx].second.first;
		int s2 = _violations[idx].second.second;
		
		snprintf(buf, 1024, "y1_%d_%d_%d_%d", j1, s1, j2, s2);
		_y1.add(IloNumVar(_env, 0, 1, buf));
		snprintf(buf, 1024, "y2_%d_%d_%d_%d", j1, s1, j2, s2);
		_y2.add(IloNumVar(_env, 0, 1, buf));
		snprintf(buf, 1024, "z_%d_%d_%d_%d", j1, s1, j2, s2);
		_z.add(IloNumVar(_env, 0, 1, buf));
	}
}

void IlpSolver::updateConstraints()
{
	int nPrevConflictConstraints = _nConflictConstraints;
	
	for (int idx = 0; idx < _nviolations; ++idx)
	{
		int j1 = _violations[idx].first.first;
		int s1 = _violations[idx].first.second;
		int j2 = _violations[idx].second.first;
		int s2 = _violations[idx].second.second;
		int curr_idx = nPrevConflictConstraints + idx;

		for (int cell_idx = 0; cell_idx < _ncells; ++cell_idx)
		{
			// y1 constraints
			_model.add(_y1[curr_idx] <= 1 - _x[cell_idx][j1][s1] + _x[cell_idx][j2][s2]);
			
			// y2 constraints
			_model.add(_y2[curr_idx] <= 1 - _x[cell_idx][j2][s2] + _x[cell_idx][j1][s1]);
			
			// z constraints
			_model.add(_z[curr_idx] <= 2 - _x[cell_idx][j1][s1] - _x[cell_idx][j2][s2]);
		}
		
		_model.add(_y1[curr_idx] + _y2[curr_idx] + _z[curr_idx] >= 1);
		_nConflictConstraints++;
	}
}

void IlpSolver::initObjective()
{
	IloExpr obj(_env);
	
	for (int j = 0; j < _nmutations; ++j)
	{
		for (int s = 0; s < _input.getMaxCharacterHomoplasy(j); ++s)
		{
//			obj += _w[j][s];
			obj += _input.getWeight(j) * _w[j][s];
		}
//		std::cout << "added cost of mutation " << j << std::endl;
	}
	
	_obj = IloMinimize(_env, obj);
	_model.add(_obj);
}

bool IlpSolver::solve(int timeLimit,
                      int memoryLimit,
                      int nrThreads,
                      double EPgap,
                      bool verbose,
                      bool allViolations)
{
  if (!verbose)
  {
    _env.setOut(_env.getNullStream());
    _env.setError(_env.getNullStream());
    _env.setWarning(_env.getNullStream());
    _cplex.setOut(_env.getNullStream());
    _cplex.setError(_env.getNullStream());
    _cplex.setWarning(_env.getNullStream());
  }
  else
  {
    _env.setOut(std::cerr);
    _env.setError(std::cerr);
    _env.setWarning(std::cerr);
    _cplex.setOut(std::cerr);
    _cplex.setError(std::cerr);
    _cplex.setWarning(std::cerr);
  }
  
  _cplex.setParam(IloCplex::ParallelMode, -1);

  if (nrThreads > 0)
  {
      std::cout << "setting nthreads " << nrThreads << std::endl;
    _cplex.setParam(IloCplex::Threads, nrThreads);
  }
  if (timeLimit > 0)
  {
//    int t = timeLimit - g_timer.realTime();
//    if (t < 0)
//      return false;
    _cplex.setParam(IloCplex::TiLim, timeLimit);
  }
  if (memoryLimit > 0)
  {
    _cplex.setParam(IloCplex::WorkMem, memoryLimit);
  }
    
  std::cout << "setting EPgap " << EPgap << std::endl; 
  _cplex.setParam(IloCplex::EpGap, EPgap);
	
  _cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);

  _cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);
	
  _nConstraints = _cplex.getNrows();

  bool res = false;
  if (allViolations) {
    _cplex.solve();

    for (int mutation_idx = 0; mutation_idx < _homoplasyVector.size(); mutation_idx++) {
      _homoplasyVector[mutation_idx] = _input.getMaxCharacterHomoplasy(mutation_idx);
    }

    setAllConstraints();
    updateHomoplasyVector();
    updateVariables();
    updateConstraints();
    updateVariableBounds();

    _cplex.solve();

	std::cerr << "CPLEX: [" << _cplex.getObjValue() << " , " << _cplex.getBestObjValue() << "]" << std::endl;
	res = (checkViolation() == 0);

    return res;
  }
  
  int iteration = 1;
  while (true)
  {
    std::cerr << "Step " << iteration << " -- number of constraints: " << _nConstraints << std::endl;
    std::cerr << "Step " << iteration << " -- number of active variables: " << _nActivated << std::endl;
		std::cerr << "Step " << iteration << " -- number of conflict constraints " << _nConflictConstraints << std::endl;


//		std::string export_filename = "export_it_" + std::to_string(iteration) + ".lp";
//		exportModel(export_filename);
		
    _cplex.solve();
    // if (_cplex.getStatus() != IloAlgorithm::Optimal || _cplex.getCplexStatus() == IloCplex::AbortTimeLim)
    // {
    //   res = false;
    //   break;
    // }
    
		std::cerr << "Step " << iteration << " -- solve complete" << std::endl;
		
		int additionalViolations = checkViolation();
//    int separatedConstraints = separate();
//    _nrConstraints += separatedConstraints;
    std::cerr << "Step " << iteration << " -- detected " << additionalViolations << " conflicts" << std::endl;
		std::cerr << "Step " << iteration << " -- number of violations: " << _nviolations << std::endl;
		std::cerr << "Step " << iteration << " -- number of violated mutations: " << _violatedMutations.size() << std::endl;

    if (additionalViolations == 0)
    {
      res = true;
      break;
    }
		else
		{
			updateHomoplasyVector();
			updateVariables();
			updateConstraints();
			updateVariableBounds();
		}
    
    ++iteration;
		
//		if (iteration == 10)
//		{
//			return false;
//		}
				
  }
  
	if (res)
	{
		std::cerr << "CPLEX: [" << _cplex.getObjValue() << " , " << _cplex.getBestObjValue() << "]" << std::endl;
	}
	
//	if (res)
//	{
//		writeSolution(output);
//	}
//	else
//	{
//		output << "no solution";
//	}
	
//  if (res)
//  {
//    processSolution();
//    std::cerr << "CPLEX: [" << _cplex.getObjValue() << " , " << _cplex.getBestObjValue() << "]" << std::endl;
//  }
//  std::cerr << "Elapsed time: " << g_timer.realTime() << std::endl;

  return res;
}

void IlpSolver::writeSolution(std::ostream& output)
{
	for (IntPair cell_mut : _input.getFreeIndices())
	{
		int cell_idx = cell_mut.first;
		int mut_idx = cell_mut.second;
		int nstates = _homoplasyVector[mut_idx];
		for (int s = 0; s < nstates; ++s)
		{
			output << cell_idx << " " << mut_idx << " " << s << " " << _cplex.getValue(_x[cell_idx][mut_idx][s]) << std::endl;
		}
	}
}
