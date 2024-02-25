#include "starhomoplasy.hpp"
#include <spdlog/spdlog.h>
#include <stack>
#include <algorithm>
#include <random>

namespace starhomoplasy {
    namespace {
        int rand_int(std::ranlux48_base& gen, int a, int b) {
            std::uniform_int_distribution<int> distrib(a, b);
            return distrib(gen);
        }
    }

    void arbitrarily_resolve_polytomies(digraph<star_homoplasy_data> &tree) {
        int internal_node_counter = 0;
        for (const auto &u : tree.nodes()) {
            if (tree[u].data.internal) {
                internal_node_counter++;
            }
        }

        for (const auto &u : tree.nodes()) {
            if (!tree[u].data.internal || tree.out_degree(u) <= 2) continue;

            std::vector<int> children;
            for (const auto &v : tree.successors(u)) {
                children.push_back(v);
            }

            for (const auto &v : children) {
                tree.remove_edge(u, v);
            }

            std::vector<int> us;
            for (size_t i = 0; i < children.size() - 2; i++) {
                int new_u = tree.add_vertex(star_homoplasy_data{
                        .character_states = std::vector<int>(tree[u].data.character_states.size(), -1),
                        .internal = true,
                        .cell = fmt::format("internal_{}", internal_node_counter++),
                        .valid = false,
                        .parsimony_score = -1
                        });

                us.push_back(new_u);
                if (i == 0) {
                    tree.add_edge(u, new_u);
                } else {
                    tree.add_edge(us[i - 1], new_u);
                }
            }

            tree.add_edge(u, children[children.size() - 1]);
            tree.add_edge(us[us.size() - 1], children[children.size() - 2]);
            for (size_t i = 0; i < children.size() - 2; i++) {
                tree.add_edge(us[i], children[i]);
            }
        }
    }

    void small_parsimony(digraph<star_homoplasy_data> &tree, 
                         const std::map<std::string, std::map<int, double>>& mutation_priors,
                         const character_state_matrix& M,
                         int vertex, int root) {
        if (tree[vertex].data.valid) {
            return;
        }

        if (!tree[vertex].data.internal) {
            tree[vertex].data.valid = true;
            tree[vertex].data.parsimony_score = 0;
            return;
        }

        const auto successors = tree.successors(vertex);

        for (const auto &u : successors) {
            small_parsimony(tree, mutation_priors, M, u, root);
        }

        std::vector<int> character_states(tree[vertex].data.character_states.size(), 0);

        if (vertex != root) {
            for (size_t j = 0; j < tree[vertex].data.character_states.size(); j++) {
                int parent_state = -1;
                bool zero_parent_state = false;
                for (const auto &u : successors) {
                    if (tree[u].data.character_states[j] == 0) {
                        zero_parent_state = true;
                        break;
                    } 

                    if (tree[u].data.character_states[j] != -1) {
                        if (parent_state == -1) {
                            parent_state = tree[u].data.character_states[j];
                            continue;
                        } 

                        if (parent_state != tree[u].data.character_states[j]) {
                            zero_parent_state = true;
                            break;
                        }
                    }
                }

                if (!zero_parent_state) {
                    character_states[j] = parent_state;
                }
            }
        }

        double parsimony_score = 0;
        for (const auto &u : successors) {
            for (size_t j = 0; j < tree[vertex].data.character_states.size(); j++) {
                if (tree[u].data.character_states[j] != -1 && tree[u].data.character_states[j] != character_states[j]) {
                    parsimony_score += mutation_priors.at(M.characters[j]).at(tree[u].data.character_states[j]);
                }
            }
        }

        for (const auto &u : successors) {
            parsimony_score += tree[u].data.parsimony_score;
        }

        tree[vertex].data.character_states = character_states;
        tree[vertex].data.valid = true;
        tree[vertex].data.parsimony_score = parsimony_score;
    }

     void nni(digraph<star_homoplasy_data>& t, int u, int w, int v, int z) {
        t.remove_edge(u, w);
        t.remove_edge(v, z);
        t.add_edge(v, w);
        t.add_edge(u, z);
    }

    void undo_nni(digraph<star_homoplasy_data>& t, int u, int w, int v, int z) {
        t.add_edge(u, w);
        t.add_edge(v, z);
        t.remove_edge(v, w);
        t.remove_edge(u, z);
    }

    void invalidate(digraph<star_homoplasy_data> &t, int root, int u) {
        t[root].data.valid = false;
        int current_node = u;
        do {
            t[current_node].data.valid = false;
            current_node = *t.predecessors(current_node).begin(); // requires only one parent exists, i.e. t is a tree
        } while (current_node != root);
    }

    void invalidate(digraph<star_homoplasy_data> &t, int root) {
        std::stack<int> callstack;
        callstack.push(root);
        while (!callstack.empty()) {
            int node = callstack.top();
            callstack.pop();

            t[node].data.valid = false;
            for (const auto& child : t.successors(node)) {
                callstack.push(child);
            }
        }
    }

    /*
      Performs all NNIs in the immediate neighborhood of the passed in
      tree and returns the best move. Does not modify the input tree.
     */
    std::optional<std::tuple<int, int, int, int>> greedy_nni(
            digraph<star_homoplasy_data> &t,
            const std::map<std::string, std::map<int, double>>& mutation_priors,
            const character_state_matrix& M,
            const std::map<int, std::pair<int, int>> &indexed_edges,
            const std::vector<int> &edge_indices,
            bool greedy
    ) {
        double best_score = t[0].data.parsimony_score; // i.e. best_score = \infty
        std::optional<std::tuple<int, int, int, int>> best_move;
        for (int idx : edge_indices) {
            const auto& [u, v] = indexed_edges.at(idx);

            if (t.successors(v).empty()) continue; // i.e if not internal

            std::vector<int> u_children;
            std::vector<int> v_children;

            for (int w : t.successors(u)) {
                if (w != v) {
                    u_children.push_back(w);
                }
            }

            for (int w : t.successors(v)) {
                v_children.push_back(w);
            }

            for (auto w : u_children) {
                for (auto z : v_children) {
                    nni(t, u, w, v, z);
                    invalidate(t, 0, v);
                    small_parsimony(t, mutation_priors, M, 0, 0);

                    double score = t[0].data.parsimony_score;
                    if (score < best_score) {
                        best_score = score;
                        best_move = std::make_tuple(u, w, v, z);

                        if (greedy) {
                            undo_nni(t, u, w, v, z);
                            invalidate(t, 0, v);
                            return best_move;
                        }
                    }

                    undo_nni(t, u, w, v, z);
                    invalidate(t, 0, v);

                }
            }
        }

        return best_move;
    }

    digraph<star_homoplasy_data> hill_climb(
        digraph<star_homoplasy_data> t, 
        const std::map<std::string, std::map<int, double>>& mutation_priors,
        const character_state_matrix& M,
        std::ranlux48_base& gen, 
        bool greedy
    ) {
        std::map<int, std::pair<int, int>> index_to_edges;
        std::map<std::pair<int, int>, int> edges_to_index;
        std::vector<int> random_indices;
        
        int idx = 0;
        for (const auto &p : t.edges()) {
            index_to_edges[idx] = p;
            edges_to_index[p] = idx;
            random_indices.push_back(idx);
            idx++;
        }
            
        std::shuffle(random_indices.begin(), random_indices.end(), gen);

        double current_score = t[0].data.parsimony_score;
        int iterations = 0;
        for (; true; iterations++) {
            auto best_move = greedy_nni(t, mutation_priors, M, index_to_edges, random_indices, greedy);
            if (!best_move) break;

            auto [u, w, v, z] = *best_move;
            nni(t, u, w, v, z);

            // update edge map by deleting (u, w) and (v, z) and adding
            // (v, w) and (u, z)
            int i1 = edges_to_index[std::make_pair(u, w)];
            int i2 = edges_to_index[std::make_pair(v, z)];
            index_to_edges[i1] = std::make_pair(v, w);
            edges_to_index[std::make_pair(v, w)] = i1;
            index_to_edges[i2] = std::make_pair(u, z);
            edges_to_index[std::make_pair(u, z)] = i2;

            invalidate(t, 0, v);
            small_parsimony(t, mutation_priors, M, 0, 0);

            double new_score = t[0].data.parsimony_score;

            if (current_score <= new_score) break;
            current_score = new_score;
        }

        return t;
    }

    digraph<star_homoplasy_data> stochastic_nni(const digraph<star_homoplasy_data>& t, std::ranlux48_base& gen, float aggression) {
        digraph<star_homoplasy_data> perturbed_t = t;

        std::vector<std::pair<int, int>> internal_edges;
        for (auto [u, v] : perturbed_t.edges()) {
            if (perturbed_t.successors(v).empty()) continue;
            internal_edges.push_back(std::make_pair(u, v));
        }

        int num_perturbations = internal_edges.size() * aggression;
        for (int i = 0; i < num_perturbations; i++) {
            internal_edges.clear();
            for (auto [u, v] : perturbed_t.edges()) {
                if (perturbed_t.successors(v).empty()) continue;
                internal_edges.push_back(std::make_pair(u, v));
            }

            int index = rand_int(gen, 0, internal_edges.size() - 1);
            auto [u, v] = internal_edges[index];

            std::vector<int> u_children;
            std::vector<int> v_children;

            for (int w : perturbed_t.successors(u)) {
                if (w != v) {
                    u_children.push_back(w);
                }
            }

            for (int w : perturbed_t.successors(v)) {
                v_children.push_back(w);
            }

            int w = u_children[rand_int(gen, 0, u_children.size() - 1)];
            int z = v_children[rand_int(gen, 0, v_children.size() - 1)];

            nni(perturbed_t, u, w, v, z);
        }

        invalidate(perturbed_t, 0);
        return perturbed_t;
    }
};
