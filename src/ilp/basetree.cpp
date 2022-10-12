/*
 * basetree.cpp
 *
 *  Created on: 08-Jan-2019
 *      Author: P. Sashittal
 */

#include "basetree.h"
#include <lemon/bfs.h>

BaseTree::BaseTree()
: _tree()
, _nodeToId(_tree)
, _idToNode()
, _nodeToIndex(_tree)
, _indexToNode()
, _arcToIndex(_tree)
, _timestamp(_tree)
, _label(_tree)
, _hostLabel()
, _nhosts()
, _unshosts()
, _nTimePoints()
, _enttime()
, _remtime()
, _lowerTime(_tree)
, _upperTime(_tree)
, _timePoints()
, _maxInf()
, _totalTime()
{
    _unshosts = 0;
}

BaseTree::BaseTree(int unhosts)
: _tree()
, _nodeToId(_tree)
, _idToNode()
, _nodeToIndex(_tree)
, _indexToNode()
, _arcToIndex(_tree)
, _timestamp(_tree)
, _label(_tree)
, _hostLabel()
, _nhosts()
, _unshosts(unhosts)
, _nTimePoints()
, _enttime()
, _remtime()
, _lowerTime(_tree)
, _upperTime(_tree)
, _timePoints()
, _maxInf()
, _totalTime()
{
}

void BaseTree::writePtree(std::ostream& out) const
{
    out << "no solution";
}


void BaseTree::writePtree(std::ostream& out, const std::string& msg) const
{
    out << msg;
}

void BaseTree::writePtree(std::ostream& out, const IntNodeMap& sol_label) const
{
    int ntotal = lemon::countNodes(_tree);
    
    for (int i = 1; i <= ntotal; ++i)
    {
        assert(_idToNode.count(std::to_string(i)));
        Node v = _idToNode.find(std::to_string(i))->second;
        if (lemon::countOutArcs(_tree, v) == 0)
        {
            // is leaf
            out << _timestamp[v] << "\t0\t0\t" << _label[v] + 1 << std::endl;
        }
        else
        {
            // is internal node
            out << _timestamp[v] << "\t";
            for (OutArcIt e(_tree, v); e != lemon::INVALID; ++e) {
                out << _nodeToId[_tree.target(e)] << "\t";
            }
            //out << _label[_idToNode[std::to_string(i)]] + 1 << std::endl;
            out <<  sol_label[v] + 1 << std::endl;
        }
    }
}


void BaseTree::writePtree(std::ostream& out,
                          const IntNodeMap& sol_label,
                          const IntArcMap& xi_label,
                          const DoubleArcMap& xi_time) const
{
    int ntotal = lemon::countNodes(_tree);
    int leaf_counter = 0; // count the number of leaves written
    //int internal_counter = 0; // count the number of internal nodes written
    int row_counter = 0; // count the number of rows written
    
    for (int i = 1; i <= ntotal; ++i)
    {
        assert(_idToNode.count(std::to_string(i)));
        Node v = _idToNode.find(std::to_string(i))->second;
        if (lemon::countOutArcs(_tree, v) == 0)
        {
            // is leaf
            out << _timestamp[v] << "\t0\t0\t" << _label[v] + 1 << std::endl;
            
            ++leaf_counter;
            ++row_counter;
        }
        else
        {
            // is internal node
            // first write xi of each outgoing arc of the node
            for (OutArcIt e(_tree, v); e != lemon::INVALID; ++e)
            {
                out << xi_time[e] << "\t";
                if (std::stoi(_nodeToId[_tree.target(e)]) <= leaf_counter)
                {
                    out << _nodeToId[_tree.target(e)] << "\t0\t";
                }
                else
                {
                    out << 3*(std::stoi(_nodeToId[_tree.target(e)])) - 2*leaf_counter << "\t0\t";
                }
                
                out << xi_label[e] + 1 << std::endl;
                
                ++row_counter;
            }
            
            out << _timestamp[v] << "\t" << row_counter - 1 << "\t" << row_counter << "\t" << sol_label[v] + 1 << std::endl;

            ++row_counter;
        }
    }
}


void BaseTree::writeTtree(std::ostream& out) const
{
    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
        Node u = _tree.source(a);
        Node v = _tree.target(a);
        
        out << _nodeToId[u] << " " << _nodeToId[v] << std::endl;
    }
}

bool BaseTree::readTtree(std::istream& in)
{    
    int idx = 0;
    
    while (in.good())
    {
        std::string line;
        getline(in, line);
        
        if (line.empty())
            break;
        
        StringVector s;
        boost::split(s, line, boost::is_any_of("\t "));
        
        if (s.size() != 4)
        {
            std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
            return false;
        }
        
        double time = std::stod(s[1]);

        _enttime.push_back(0);
        _remtime.push_back(time);
        
        idx++;
    }
    
    /// initialize nhosts
    _nhosts = idx;
    
    /// initiaize totalTime
    _totalTime = *std::max_element(_remtime.begin(), _remtime.end());
    
    /*
    std::cout << "number of hosts is: " << _nhosts << std::endl;
    std::cout << "removal times are: " << std::endl;
    for (int s = 0; s < _nhosts; ++s)
    {
        std::cout << "host " << s+1 << " has " << _remtime[s] << std::endl;
    }
    
    std::cout<< "total time is " << _totalTime << std::endl;
    */
    
    return true;
}

void BaseTree::writeDOT(const ArcMatrix& N,
                        const IntNodeMap& l,
                        std::ostream& out) const
{
  StringVector colorMap({
    "#3243BA",
    "#0363E1",
    "#0D75DC",
    "#1485D4",
    "#0998D1",
    "#06A7C6",
    "#15B1B4",
    "#38B99E",
    "#65BE86",
    "#92BF73",
    "#B7BD64",
    "#D9BA56",
    "#F8BB44",
    "#FCCE2E",
    "#F5E41D",
    "#F9FB0E"
  });
  
  out << "digraph N {" << std::endl;
  out << "\toverlap=\"false\"" << std::endl;
  out << "\trankdir=\"LR\"" << std::endl;
  for (int s = 0; s < _nhosts; ++s)
  {
    out << "\t" << s << " [label=\"" << _hostLabel[s] << "\",width=1.2,height=1.2,style=\"\",penwidth=3,color=\"" << colorMap[s % colorMap.size()] << "\"]" << std::endl;
  }
  
  for (const ArcVector& event : N)
  {
    int s = l[_tree.source(event.front())];
    int t = l[_tree.target(event.front())];
    
    out << "\t" << s << " -> " << t << " [penwidth=1,color=black,label=\"" << event.size() << "\"]" << std::endl;
  }
  
  out << "}" << std::endl;
}

bool BaseTree::readHost(std::istream& in, int nTimePoints)
{
    // TODO: Clear everything before starting
    int idx = 0;
    _hostLabel.clear();
    while (in.good())
    {
        std::string line;
        getline(in, line);
        
        if (line.empty())
            break;
        
        StringVector s;
        boost::split(s, line, boost::is_any_of("\t "));
        
        if (s.size() != 3)
        {
            std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
            return false;
        }
        
        double etime = std::stod(s[1]);
        double rtime = std::stod(s[2]);
      
        _hostLabel.push_back(s[0]);
        _enttime.push_back(etime);
        _remtime.push_back(rtime);
        
        idx++;
    }
    
    /// initiaize totalTime
    _totalTime = *std::max_element(_remtime.begin(), _remtime.end());
    
    /// add unsampled hosts
    for (int s = 0; s < _unshosts; ++s)
    {
        _enttime.push_back(0);
        _remtime.push_back(_totalTime);
        idx++;
    }
    
    // initialize nhosts
    _nhosts = idx;

		// initialize the timepoints
		_nTimePoints = nTimePoints;
		_timePoints = DoubleMatrix(_nhosts);
		
		for (int s = 0; s < _nhosts; ++s)
		{
				_timePoints[s] = DoubleVector(nTimePoints);
				double timeStep = (_remtime[s] - _enttime[s])/double(nTimePoints - 1);
				
				for (int i = 0; i < nTimePoints; ++i)
				{
						_timePoints[s][i] = _enttime[s] + timeStep*double(i);
				}
		}
    
    std::cout << "number of hosts is: " << _nhosts << std::endl;
    /*
    std::cout << "removal times are: " << std::endl;
    for (int s = 0; s < _nhosts; ++s)
    {
        std::cout << "host " << s+1 << " has " << _enttime[s] << " and " << _remtime[s] << std::endl;
    }
    */
    std::cout<< "total time is " << _totalTime << std::endl;
		
		
		// output timepoint for each host
		std::cout << "time points for each host are: " << std::endl;
		for (int s = 0; s < _nhosts; ++s)
		{
				std::cout << "host " << s + 1;
				for (int i = 0; i < nTimePoints; ++i)
				{
						std::cout << " " << _timePoints[s][i];
				}
				std::cout << std::endl;
		}
		
    return true;
}

bool BaseTree::readPtree(std::istream& in, bool binaryTree)
{
    _idToNode.clear();
    
    int idx = 1;
    
    while (in.good())
    {
        std::string line;
        getline(in, line);
        
        if (line.empty())
            break;
        
        StringVector s;
        boost::split(s, line, boost::is_any_of("\t "));
        
        if (binaryTree)
        {
            if (s.size() != 4)
            {
                std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
                return false;
            }
            
            double time = std::stod(s[0]);
            std::string child1 = s[1];
            std::string child2 = s[2];
            int vlabel = std::stoi(s[3]) - 1;
            
            if((child1 == "0") && (child2 == "0"))
            {
                // node is a leaf
                Node u = _tree.addNode();
                _label[u] = vlabel;
                _timestamp[u] = time;
                _nodeToId[u] = std::to_string(idx);
                _nodeToIndex[u] = idx - 1;
								_indexToNode.push_back(u);
                _idToNode[std::to_string(idx)] = u;
            }
            else
            {
                // node is internal
                Node u = _tree.addNode();
                _label[u] = vlabel;
                _timestamp[u] = time;
                _tree.addArc(u, _idToNode[child1]);
                _tree.addArc(u, _idToNode[child2]);
                _nodeToId[u] = std::to_string(idx);
                _nodeToIndex[u] = idx - 1;
								_indexToNode.push_back(u);
                _idToNode[std::to_string(idx)] = u;
            }
        }
        else
        {
            if (s.size() < 2)
            {
                std::cerr << "Error: line '" << line << "' incorrect number of columns" << std::endl;
                return false;
            }

            double time = std::stod(s[0]);
            int vlabel = std::stoi(s[1]) - 1;

            
            if (s.size() == 2)
            {
                // is leaf
                Node u  = _tree.addNode();
                _label[u] = vlabel;
                _timestamp[u] = time;
                _nodeToId[u] = std::to_string(idx);
                _nodeToIndex[u] = idx - 1;
								_indexToNode.push_back(u);
                _idToNode[std::to_string(idx)] = u;
            }
            else
            {
                // is internal node
                Node u = _tree.addNode();
                _label[u] = vlabel;
                _timestamp[u] = time;
                
                for (int i = 2; i < s.size(); ++i)
                {
                    _tree.addArc(u, _idToNode[s[i]]);
                }
                _nodeToId[u] = std::to_string(idx);
                _nodeToIndex[u] = idx - 1;
								_indexToNode.push_back(u);
                _idToNode[std::to_string(idx)] = u;
            }
        }
        ++idx;
    }
    
    // check that tree has a single root
    _root = lemon::INVALID;
    for (NodeIt node(_tree); node != lemon::INVALID; ++node)
    {
        if (InArcIt(_tree, node) == lemon::INVALID)
        {
            if (_root != lemon::INVALID)
            {
                std::cerr << "Error: multiple root node '" << _nodeToId[node]
                << "' and '" << _nodeToId[_root] << "'" << std::endl;
                return false;
            }
            _root = node;
        }
    }
 
    // check if tree is empty
    if (_idToNode.empty())
    {
        std::cerr << "Error: empty tree" << std::endl;
        return false;
    }
    
    // initialize arc-to-index
    idx = 0;
    for (ArcIt a(_tree); a!= lemon::INVALID; ++a)
    {
        _arcToIndex[a] = idx;
        ++idx;
    }
		
		// initialize the lower and upper time variables
		for (ArcIt a(_tree); a != lemon::INVALID; ++a)
		{
				_lowerTime[a] = DoubleVector(_nhosts);
				
				for (int s = 0; s < _nhosts; ++s)
				{
						double sourceTime = getTime(_tree.source(a));
						if (sourceTime <= _enttime[s])
						{
								_lowerTime[a][s] = 0;
						}
						else
						{
								_lowerTime[a][s] = sourceTime - _enttime[s];
						}
				}
		}

		for (ArcIt a(_tree); a != lemon::INVALID; ++a)
		{
				_upperTime[a] = DoubleVector(_nhosts);
				
				for (int s = 0; s < _nhosts; ++s)
				{
						double targetTime = getTime(_tree.target(a));
						if (targetTime >= _remtime[s])
						{
								_upperTime[a][s] = 0;
						}
						else
						{
								_upperTime[a][s] = targetTime - _remtime[s];
						}
				}
		}
		
		/*
    // some print statements for debugging
    for (NodeIt v(_tree); v != lemon::INVALID; ++v)
    {
        std::cout << _nodeToId[v] << " has name: " << _nodeToId[v] <<
        " and time: " << _timestamp[v] << std::endl;
    }

    for (ArcIt a(_tree); a != lemon::INVALID; ++a)
    {
        Node s = _tree.source(a);
        Node t = _tree.target(a);

        std::cout << "Arc: (" << _nodeToId[s] << ","
        << _nodeToId[t] << ")" << std::endl;
    }
		*/
		
		std::cout << "root is " << _nodeToId[_root] << std::endl;

		/*
		// output the upper and lower limits for each (arc,host)
		for (ArcIt a(_tree); a != lemon::INVALID; ++a)
		{
				for (int s = 0; s < _nhosts; ++s)
				{
						Node ai = _tree.source(a);
						Node aj = _tree.target(a);
						
						std::cout << "Arc: (" << _nodeToId[ai] << "," << _nodeToId[aj] << ")" << " for host " << s+1 << " has [" << _lowerTime[a][s] << "," << _upperTime[a][s] << "]" <<std::endl;
				}
		}
		*/
    
    return true;
}

int BaseTree::getMu(const IntNodeMap& l) const
{
  int mu = 0;
  for (ArcIt a(_tree); a != lemon::INVALID; ++a)
  {
    Node u = _tree.source(a);
    Node v = _tree.target(a);
    
    if (l[u] != l[v])
      ++mu;
  }
  
  return mu;
}

ArcMatrix BaseTree::getN(const IntNodeMap& l) const
{
  ArcMatrix N;
  const int nrInfectedHosts = getNrHost();
  
  int gamma = 0;
  for (int s = 0; s < nrInfectedHosts; ++s)
  {
    for (int t = 0; t < nrInfectedHosts; ++t)
    {
      if (s != t)
      {
        // timeIntervals - ((endtime, starttime), arc)
        std::vector<std::pair<std::pair<double, double>, Arc> > timeIntervals;
        
        for (ArcIt eij(_tree); eij != lemon::INVALID; ++eij)
        {
          Node vi = _tree.source(eij);
          Node vj = _tree.target(eij);
          
          if (l[vi] == s && l[vj] == t) {
            timeIntervals.push_back(std::make_pair(std::make_pair(getTime(vj), getTime(vi)), eij));
          }
        }
        
        while (timeIntervals.size() > 0)
        {
          auto minTime = *std::min_element(timeIntervals.cbegin(), timeIntervals.cend());
          std::vector<std::pair<std::pair<double, double>, Arc> > timeIntervalsToRemove;

          N.push_back(ArcVector());
          for (const auto& x : timeIntervals)
          {
            if (x.first.second <= minTime.first.first)
            {
              N.back().push_back(x.second);
            }
          }

          timeIntervals.erase(std::remove_if(timeIntervals.begin(), timeIntervals.end(),
                                             [minTime](const std::pair<std::pair<double, double>, Arc>& x){
                                               return x.first.second <= minTime.first.first;
                                             }), timeIntervals.end());
          
          ++gamma;
        }
      }
    }
  }
  
  return N;
}

/*
void BaseTree::readInputWithHosts(std::istream& hfile,
                                  std::istream& pfile,
                                  int unsampled_hosts,
                                  bool binaryTree,
                                  BaseTree& B)
{
  // construct the basetree
  B.setUnHosts(unsampled_hosts);
  
  // read the host file
  B.readHost(hfile);
  
  for (int i = 0; i < unsampled_hosts; ++i)
  {
    B._hostLabel.push_back("Unsampled");
  }
  
  // read the pfile
  B.readPtree(pfile, binaryTree);
}
*/

void BaseTree::readInputWithTransmission(std::istream& tfile,
                                         std::istream& pfile,
                                         BaseTree& B)
{
  // read the tfile
  B.readTtree(tfile);
  
  // read the pfile
  B.readPtree(pfile, false);
}


int BaseTree::getGamma(const IntNodeMap& l) const
{
  return getN(l).size();
}
