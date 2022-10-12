/*
 * utils.h
 *
 *  Created on: 9-nov-2018
 *      Author: P. Sashittal
 */

#ifndef UTILS_H
#define UTILS_H

#include <lemon/list_graph.h>
#include <lemon/adaptors.h>
#include <lemon/tolerance.h>
#include <cassert>
#include <iostream>
#include <set>
#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <random>

typedef lemon::ListDigraph Digraph;
typedef lemon::ListGraph Graph;
DIGRAPH_TYPEDEFS(Digraph);
typedef Digraph::NodeMap<Digraph::Node> NodeNodeMap;
typedef Digraph::NodeMap<std::string> StringNodeMap;
typedef Digraph::NodeMap<double> DoubleNodeMap;
typedef Digraph::NodeMap<uint64_t> UInt64NodeMap;
typedef Digraph::ArcMap<uint64_t> UInt64ArcMap;
typedef std::vector<uint64_t> UInt64Vector;
typedef Digraph::ArcMap<UInt64Vector> UInt64VectorArcMap;
typedef Digraph::NodeMap<UInt64Vector> UInt64VectorNodeMap;
typedef std::map<std::string, Node> StringToNodeMap;
typedef std::map<std::string, int> StringToIntMap;
typedef std::set<Node> NodeSet;
typedef std::set<Arc> ArcSet;
typedef std::vector<Arc> ArcVector;
typedef std::vector<ArcVector> ArcMatrix;
typedef Digraph::NodeMap<NodeSet> NodeNodeSetMap;
typedef std::pair<int, Node> IntNodePair;
typedef std::list<IntNodePair> IntNodePairList;
typedef std::vector<Node> NodeVector;
typedef std::vector<NodeVector> NodeMatrix;
typedef std::list<Node> NodeList;
typedef NodeList::const_iterator NodeListIt;
typedef std::vector<NodeList> NodeListVector;
typedef NodeListVector::const_iterator NodeListVectorIt;
typedef std::list<NodeListIt> NodeListItList;
typedef std::vector<NodeListIt> NodeListItVector;
typedef std::list<NodeList> NodeListList;
typedef std::vector<std::string> StringVector;
typedef std::pair<Node, Node> NodePair;
typedef std::list<NodePair> NodePairList;
typedef Digraph::NodeMap<NodePair> NodePairMap;
typedef std::set<std::string> StringSet;
typedef std::multiset<StringSet> Split;
typedef std::set<Split> SplitSet;
typedef std::set<StringSet> StringSetSet;
typedef std::vector<int> IntVector;
typedef std::vector<IntVector> IntMatrix;
typedef std::vector<IntMatrix> IntTensor;
typedef std::vector<bool> BoolVector;
typedef std::vector<BoolVector> BoolMatrix;
typedef std::vector<BoolMatrix> Bool3Matrix;
typedef std::vector<double> DoubleVector;
typedef std::vector<DoubleVector> DoubleMatrix;
typedef std::vector<DoubleMatrix> DoubleTensor;
typedef std::pair<double, double> DoublePair;
typedef std::vector<DoublePair> DoublePairVector;
typedef std::set<int> IntSet;
typedef std::list<int> IntList;
typedef IntSet::const_iterator IntSetIt;
typedef std::list<std::string> StringList;
typedef std::set<std::string> StringSet;
typedef std::pair<std::string, std::string> StringPair;
typedef std::list<StringPair> StringPairList;
typedef Digraph::NodeMap<StringList> StringListNodeMap;
typedef Digraph::NodeMap<BoolVector> BoolVectorNodeMap;
typedef Digraph::NodeMap<IntVector> IntVectorNodeMap;
typedef Digraph::NodeMap<DoubleVector> DoubleVectorNodeMap;
typedef Digraph::NodeMap<IntSet> IntSetNodeMap;
typedef Digraph::ArcMap<IntSet> IntSetArcMap;
typedef std::pair<int, int> IntPair;
typedef std::set<IntPair> IntPairSet;
typedef std::vector<IntSet> IntSetVector;
typedef Digraph::NodeMap<IntPair> IntPairNodeMap;
typedef Digraph::NodeMap<IntPairSet> IntPairSetNodeMap;
typedef std::map<IntPair, Node> IntPairToNodeMap;
typedef std::map<IntPair, NodeSet> IntPairToNodeSetMap;
typedef std::pair<int, IntPair> IntTriple;
typedef std::pair<IntPair, IntPair> IntQuad;
typedef std::vector<IntQuad> IntQuadVector;
typedef std::pair<int, std::string> IntStringPair;
typedef std::vector<IntStringPair> IntStringPairVector;
typedef std::vector<IntStringPairVector> IntStringPairMatrix;
typedef std::pair<uint64_t, std::string> UInt64StringPair;
typedef std::vector<UInt64StringPair> UInt64StringPairVector;
typedef std::vector<UInt64StringPairVector> UInt64StringPairMatrix;
typedef std::vector<IntPair> IntPairVector;
typedef Digraph::ArcMap<IntVector> IntVectorArcMap;
typedef std::pair<double, Arc> DoubleArcPair;
typedef std::vector<DoubleArcPair> DoubleArcPairVector;
typedef std::vector<IntTriple> IntTripleVector;
typedef std::pair<uint64_t, int> UInt64IntPair;
typedef std::vector<UInt64IntPair> UInt64IntPairVector;
typedef std::map<IntPairSet, int> IntPairSetIntMap;
typedef Digraph::ArcMap<DoubleVector> DoubleVectorArcMap;
typedef Digraph::NodeMap<int> IntNodeMap;
typedef Digraph::NodeMap<double> DoubleNodeMap;
typedef Digraph::ArcMap<int> IntArcMap;

/// Get line from stream in a platform independent manner
std::istream& getline(std::istream& is, std::string& t);

/// Random number generator
extern std::mt19937 g_rng;

/// Tolerance for floating point comparisons
extern lemon::Tolerance<double> g_tol;

#endif // UTILS_H
