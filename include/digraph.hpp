#ifndef _DIGRAPH_H
#define _DIGRAPH_H

#include <map>
#include <unordered_map>
#include <random>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

/*
  Simple adjacency list representation of a DiGraph where the 
  vertices are represented as integers from 0 to N - 1. And
  the vertices have data associated with them.
*/

template <class T>
class vertex {
public:
    int id;
    T data;
    vertex() {};
    vertex(int id, T data) : id(id), data(data) {};
};

template <class T>
class digraph {
private:
    int id_counter = 0;

    std::unordered_map<int, std::set<int>> succ;
    std::unordered_map<int, std::set<int>> pred;

    std::unordered_map<int, vertex<T>> vertices;
public:
    // returns id of created vertex
    int add_vertex(T data) {
        vertex<T> v(id_counter, data);
        vertices[v.id] = v;
        succ[v.id] = std::set<int>();
        pred[v.id] = std::set<int>();
        id_counter++;
        return v.id;
    }

    void add_edge(int u, int v) {
        succ[u].insert(v);
        pred[v].insert(u);
    }

    void remove_edge(int u, int v) {
        succ[u].erase(v);
        pred[v].erase(u);
    }

    std::vector<int> nodes() const {
        std::vector<int> vertices;
        for (int i = 0; i < id_counter; i++) {
            vertices.push_back(i);
        }
        return vertices;
    }

    std::vector<std::pair<int, int>> edges() const {
        std::vector<std::pair<int, int>> edges;
        for (const auto &[u, vs] : succ) {
            for (const auto &v : vs) {
                edges.push_back(std::make_pair(u, v));
            }
        }
        return edges;
    }

    /* TODO
       WARNING: does not maintain invariant
       that all edges are between 0 and N. We 
       will need to update the code to fix this.
       Probably should return a map from old to new
       vertex labels.
     */
    void delete_vertex(int u) {
        // removes u and all (v, u) edges
        vertices.erase(u);
        for (const vertex<T>& v : pred[u]) {
            succ[v].erase(u);
        }

        succ.erase(u);
        pred.erase(u);
    }

    vertex<T>& operator[](int u) {
        return vertices.at(u);
    }

    const vertex<T>& operator[](int u) const {
        return vertices.at(u);
    }

    const std::set<int>& predecessors(int u) const {
        return pred.at(u);
    }

    const std::set<int>& successors(int u) const {
        return succ.at(u);
    }

    bool contains(int u) const {
        return vertices.find(u) != vertices.end();
    }

    size_t out_degree(int u) const {
        return succ.at(u).size();
    }

    /*
     * Returns true if u is an ancestor of v in the given tree.
     */
    friend bool ancestor(const digraph<T>& tree, int u, int v) {
        if (u == v) {
            return true;
        }

        for (int w : tree.succ.at(u)) {
            if (ancestor(tree, w, v)) {
                return true;
            }
        }

        return false;
    }
};

template <class T>
std::string to_adjacency_list(const digraph<T>& G, const std::unordered_map<int, int>& vertex_map) {
    std::stringstream ss;
    for (const auto& u : G.nodes()) {
        ss << vertex_map.at(u) << " ";
        for (const auto& v : G.successors(u)) {
            ss << vertex_map.at(v) << " ";
        }
        ss << std::endl;
    }
    return ss.str();
}

#endif
