/*
 * matrix.h
 *
 *  Created on: 31-Jul-2022
 *      Author: P. Sashittal
 */

#ifndef MATRIX_H
#define MATRIX_H

#include "utils.h"

/// This class holds the input to the lineage tracing problem
class Matrix
{
public:
    /// Default constructor
    Matrix();
  
		Matrix(int n, int m);
	
    /// Destructor
    virtual ~Matrix()
    {
    }
    
    /// read character matrix
//	virtual bool readCharacterMatrix(std::istream& in);
    
	virtual bool readMaxAllowedHomoplasy(std::istream& in);
    
//	virtual bool readPrior(std::istream& in);
  
	virtual bool readCharacterMutationMapping(std::istream& in);
	
	virtual bool readWeights(std::istream& in);
	
	virtual bool readMissingIndices(std::istream& in);

	virtual bool readOneIndices(std::istream& in);
	
	virtual void initFreeIndices();
	
//     struct Violation
//     {
//     public:
//         Violation(int j1, int j2, int s1, int s2)
//         : _j1(j1)
//         , _j2(j2)
//         , _s1(s1)
//         , _s2(s2)
//         {
//         }
        
//         # mutation 1
//         const int j1;
//         # mutation 2
//         const int j2;
//         # state 1
//         const int s1;
//         # state 2
//         const int s2;
//     }
    
    
	void writeMatrix(std::ostream& out);
	
	int getNumCells() const
	{
			return _ncells;
	}
	
	int getNumMutations() const
	{
			return _nmutations;
	}
	
	int getNumCharacters() const
	{
		return _ncharacters;
	}
	
	int getMaxCharacterHomoplasy(int j) const
	{
		return _maxAllowedHomoplasy[j];
	}
	
	IntPairVector getOneIndices() const
	{
		return _oneIndices;
	}
	
	IntPairVector getMissingIndices() const
	{
		return _missingIndices;
	}
	
	IntPairVector getFreeIndices() const
	{
		return _freeIndices;
	}
	
	IntVector getCharacterMutations(int p) const
	{
		return _characterMutationMap[p];
	}
	
	double getWeight(int j) const
	{
		return _weights[j];
	}
	
	int getMaxHomoplasy() const
	{
		return *std::max_element(_maxAllowedHomoplasy.begin(), _maxAllowedHomoplasy.end());
	}
	
	bool isFree(int i, int j) const
	{
		return _freeIndexIndicator[i][j];
	}
	
protected:
	/// number of cells
	int _ncells;
	
	/// number of mutations
	int _nmutations;
	
	/// number of characters;
	int _ncharacters;

	/// character-state matrix
	IntMatrix _characterMatrix;
	
	/// max allowed homoplasy
	IntVector _maxAllowedHomoplasy;

	/// 1 cell mutation index pairs
	IntPairVector _oneIndices;
	
	/// missing cell character index pairs
	IntPairVector _missingIndices;

	/// free cell mutation index pairs
	IntPairVector _freeIndices;
	
	/// free cell mutation index indicator
	BoolMatrix _freeIndexIndicator;
	
	/// character mutation mapping
	IntMatrix _characterMutationMap;
	
	/// weights vector
	DoubleVector _weights;
};

#endif // MATRIX_H
