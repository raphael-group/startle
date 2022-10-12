/*
 * ilpsolver.h
 *
 *  Created on: 1-aug-2022
 *      Author: P. Sashittal
 */

#ifndef ILPSOLVER_H
#define ILPSOLVER_H

#include "utils.h"
#include "matrix.h"
#include <fstream>
#include <ilcplex/ilocplex.h>

class IlpSolver
{
public:
	/// constructor
	IlpSolver(const Matrix& input);
	
	/// initialize the solver without repeat variables
	void init();
	/// solve
//	virtual bool solve(std::ostream& output, int timeLimit, int memoryLimit, int nThreads, bool verbose);
	virtual bool solve(int timeLimit, int memoryLimit, int nThreads, double EPgap, bool verbose, bool allViolations);
	/// write solution
	virtual void writeSolution(std::ostream& output);

protected:

	/// initialize variables and indices
	virtual void initVariables();
	/// update virable activation
//	virtual void updateVariableActivation() = 0;
	///  update vairable bounds
	virtual void updateVariableBounds();
	/// initialize constraints
	virtual void initConstraints();
	/// update variables
	virtual void updateVariables();
	/// update constraints
	virtual void updateConstraints();
    // set all constrains
    virtual void setAllConstraints();
	/// initialize objective
	virtual void initObjective();
	/// check violation
	virtual int checkViolation();
	/// update homoplasy vector
	virtual void updateHomoplasyVector();
	
	/// clear the current model
	virtual void clearModel()
	{
			_cplex.clearModel();
			_cplex.end();
			_model.end();
			_model = IloModel(_env);
			_cplex = IloCplex(_model);
	}
	/// export the current model
	void exportModel(const std::string& filename) const
	{
			_cplex.exportModel(filename.c_str());
	}
	
  typedef IloArray<IloBoolVarArray> IloBoolVarMatrix;
  typedef IloArray<IloBoolVarMatrix> IloBoolVar3Matrix;
	typedef IloArray<IloNumVarArray> IloNumVarMatrix;

protected:
	/// input data
	const Matrix& _input;
	/// number of cells
	const int _ncells;
	/// number of mutations
	const int _nmutations;
	/// number of characters;
	const int _ncharacters;
	
	/// violated mutation-state pairs
	IntQuadVector _violations;
	/// number of violations
	int _nviolations;
	/// violated mutations
	IntPairSet _violatedMutations;

	/// bool 3d matrix indicating active variables
//	Bool3Matrix _activeVariables;
	/// number of active binary variables;
	int _nActivated;
	/// number of constraints
	int _nConstraints;
	/// number of conflict constraints
	int _nConflictConstraints;
	///  homoplasy vector for each mutation
	IntVector _homoplasyVector;
	
	/// cplex environment
	IloEnv _env;
	/// cplex model
	IloModel _model;
	/// cplex solver
	IloCplex _cplex;
	/// Objective
	IloObjective _obj;
	/// x[i,j,s] =1 if and only if a[i,j]=s
	IloBoolVar3Matrix _x;
	///  y1[j,s,j',s'] =1 if and only if (j,s) is within (j',s')
	IloNumVarArray _y1;
	///  y2[j,s,j',s'] =1 if and only if (j',s') is within (j,s)
	IloNumVarArray _y2;
	/// z[j,s,j',s'] = 1 if and only if (j,s) and (j',s') are disjoint
	IloNumVarArray _z;
	/// w[j,s] = 1 if and only if at least one cell has mutation (j,s)
	IloNumVarMatrix _w;
};

#endif // ILPSOLVER_H
