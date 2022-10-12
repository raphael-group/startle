/*
 * matrix.cpp
 *
 *  Created on: 31-Jul-2022
 *      Author: P. Sashittal
 */

#include "matrix.h"

Matrix::Matrix()
: _ncells(0)
, _nmutations(0)
, _characterMatrix()
, _maxAllowedHomoplasy()
, _oneIndices()
, _missingIndices()
, _characterMutationMap()
, _weights()
{
}

Matrix::Matrix(int n, int m)
: _ncells(n)
, _nmutations(m)
, _characterMatrix()
, _maxAllowedHomoplasy()
, _oneIndices()
, _missingIndices()
, _characterMutationMap()
, _weights()
{
}

//bool Matrix::readCharacterMatrix(std::istream& in)
//{
//	int idx = 0;
//	_characterMatrix = IntMatrix(_ncells);
//	while (in.good())
//	{
//		std::string line;
//		getline(in, line);
//
//		if (line.empty())
//			break;
//
//		StringVector s;
//		boost::split(s, line, boost::is_any_of("\t "));
//
//		if (s.size() != _nmutations)
//		{
//			std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
//			return false;
//		}
//
//		_characterMatrix[idx] = IntVector(_nmutations);
//
//		for (int jdx = 0; jdx < _nmutations; ++jdx)
//		{
//			_characterMatrix[idx][jdx] = atoi(s[jdx].c_str());
//		}
//
//		++idx;
//	}
//
//	return true;
//}

bool Matrix::readMaxAllowedHomoplasy(std::istream &in)
{
	int idx = 0;
	_maxAllowedHomoplasy = IntVector(_nmutations);
	while(in.good())
	{
		std::string line;
		getline(in, line);
		
		if (line.empty())
			break;
		
		StringVector s;
		boost::split(s, line, boost::is_any_of("\t "));

		if (s.size() != 1)
		{
			std::cerr << "Error: line '" << line << "' incorrect number of entries for maximum allowed homoplasy for mutation" << std::endl;
			return false;
		}
		
		_maxAllowedHomoplasy[idx] = atoi(s[0].c_str());
		
		++idx;
	}
	
	return true;
}

bool Matrix::readOneIndices(std::istream &in)
{
	int idx = 0;
	_oneIndices = IntPairVector();
	while(in.good())
	{
		std::string line;
		getline(in, line);
		
		if (line.empty())
			break;
		
		StringVector s;
		boost::split(s, line, boost::is_any_of("\t "));

		if (s.size() != 2)
		{
			std::cerr << "Error: line '" << line << "' incorrect number of entries for one indices" << std::endl;
			return false;
		}
		
		int cell_idx = atoi(s[0].c_str());
		int mutation_idx = atoi(s[1].c_str());
		
		_oneIndices.push_back(IntPair(cell_idx, mutation_idx));
		++idx;
	}
	
	std::cout << "read " << _oneIndices.size() << " index pairs with 1 entry in the character-state matrix" << std::endl;
	
//	for (IntPair cell_mut : _oneIndices)
//	{
//		std::cout << cell_mut.first << " " << cell_mut.second << std::endl;
//	}
	
	return true;
}

bool Matrix::readMissingIndices(std::istream &in)
{
	int idx = 0;
	_missingIndices = IntPairVector();
	while(in.good())
	{
		std::string line;
		getline(in, line);
		
		if (line.empty())
			break;
		
		StringVector s;
		boost::split(s, line, boost::is_any_of("\t "));

		if (s.size() != 2)
		{
			std::cerr << "Error: line '" << line << "' incorrect number of entries for one indices" << std::endl;
			return false;
		}
		
		int cell_idx = atoi(s[0].c_str());
		int character_idx = atoi(s[1].c_str());
		
		_missingIndices.push_back(IntPair(cell_idx, character_idx));
		++idx;
	}
	
	std::cout << "read " << idx << " index pairs with missing entry in the character-state matrix" << std::endl;
	
	return true;
}

bool Matrix::readCharacterMutationMapping(std::istream &in)
{
	int idx = 0;
	_characterMutationMap = IntMatrix();
	while(in.good())
	{
		std::string line;
		getline(in, line);
		
		if (line.empty())
			break;
		
		StringVector s;
		boost::split(s, line, boost::is_any_of("\t "));
		
		int numMutations = s.size();
		if (numMutations > _nmutations)
		{
			std::cerr << "character " << idx << " has more mutations than nmutations which is " << _nmutations << std::endl;
			return false;
		}
		
		IntVector mutationVector(numMutations);
		for (int jdx = 0; jdx < numMutations; ++jdx)
		{
			int mutation_idx = atoi(s[jdx].c_str());
			if (mutation_idx > _nmutations)
			{
				std::cerr << "for character " << idx << " mutation " << mutation_idx << " is out of bounds since nmutations is " << _nmutations << std::endl;
				return false;
			}
			mutationVector[jdx] = mutation_idx;
		}
		
		_characterMutationMap.push_back(mutationVector);
		
		++idx;
	}
	
	std::cout << "read mutation map for " << idx << " characters" << std::endl;
	
	_ncharacters = idx;
	
	return true;
}

bool Matrix::readWeights(std::istream &in)
{
	int idx = 0;
	_weights = DoubleVector(_nmutations);
	while (in.good())
	{
		std::string line;
		getline(in, line);
		
		if (line.empty())
			break;
		
		StringVector s;
		boost::split(s, line, boost::is_any_of("\t "));
		
		double weight = std::stod(s[0].c_str());
//		std::cout << "mutation " << idx << " weight is " << weight << std::endl;
		
		_weights[idx] = weight;
		idx ++;
	}
	
	return true;
}

void Matrix::initFreeIndices()
{
	int idx = 0;
	
	_freeIndexIndicator = BoolMatrix(_ncells);
	for (int i = 0; i < _ncells; ++i)
	{
		_freeIndexIndicator[i] = BoolVector(_nmutations, false);
	}
	
	for (IntPair cell_mut : _oneIndices)
	{
		int cell_idx = cell_mut.first;
		int mut_idx = cell_mut.second;
		_freeIndexIndicator[cell_idx][mut_idx] = true;
		_freeIndices.push_back(cell_mut);
		idx++;
	}
	
	for (IntPair cell_character : _missingIndices)
	{
		int cell_idx  = cell_character.first;
		int character_idx = cell_character.second;
		for (int mut_idx : getCharacterMutations(character_idx))
		{
			_freeIndexIndicator[cell_idx][mut_idx] = true;
			_freeIndices.push_back(std::make_pair(cell_idx, mut_idx));
			idx ++;
		}
	}
	
	std::cout << "number of free indices: " << idx << std::endl;
}
