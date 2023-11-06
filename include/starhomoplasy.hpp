#ifndef STARHOMOPLASY_HPP
#define STARHOMOPLASY_HPP

#include <vector>
#include <optional>
#include <string>
#include "digraph.hpp"

namespace starhomoplasy {
    /* 
    * A character-state matrix $M$ is a $n \times m$ matrix where $n$ is the number of cells
    * and $m$ is the number of characters. Each cell $i$ is associated with a string $c_i$ and
    * each character $j$ is associated with a string $s_j$. The entry $M_{ij}$ is the character
    * state of character $j$ in cell $i$.
    */
    struct character_state_matrix {
        std::vector<std::vector<int>> matrix;
        std::vector<std::string> characters;
        std::vector<std::string> cells;
        size_t num_characters;
        size_t num_cells;
    };

    struct star_homoplasy_data {
        std::vector<int> character_states;
        bool internal = false;
        std::optional<std::string> cell;
        bool valid = false;
        double parsimony_score = 0;
    };

    /* 
    * Converts a non-binary phylogenetic tree to a binary
    * phylogenetic tree on the same leaf set. 
    */
    void arbitrarily_resolve_polytomies(digraph<star_homoplasy_data> &tree);

    /*
    * Fills in the optiomal star homoplasy small parsimony labeling for all subtrees of a tree $T$ 
    * rooted at $vertex$. Avoids recomputing the optimal labeling for subtrees that have already been
    * processed.
    */
    void small_parsimony(digraph<star_homoplasy_data> &tree, const std::map<std::string, std::map<int, double>>& mutation_priors,
                         const character_state_matrix& M, int vertex, int root);

    /*
      Performs (or undos) a NNI operation on edges (u, w) and (v, z) by
      swapping the edges.
      
      Does not necessarily maintain validity, but it can be
      re-maintained by calling invalidate(t, u) and invalidate(t, w).
    */
    void nni(digraph<star_homoplasy_data>& t, int u, int w, int v, int z);
    void undo_nni(digraph<star_homoplasy_data>& t, int u, int w, int v, int z);

    /*
      Invalidates all vertices on path from root to u.
     */
    void invalidate(digraph<star_homoplasy_data> &t, int root, int u);

    /*
      Invalidates all vertices in sub-tree rooted at root.
     */
    void invalidate(digraph<star_homoplasy_data> &t, int root);

    /*
      Assumes the root is the 0 vertex.
    */
    digraph<star_homoplasy_data> stochastic_nni(const digraph<star_homoplasy_data>& t, std::ranlux48_base& gen, float aggression);


    /*
    * Performs hill climbing on the small parsimony score of the input tree 
    * until no more improvement is found using NNI operations
    */
    digraph<star_homoplasy_data> hill_climb(
        digraph<star_homoplasy_data> t, 
        const std::map<std::string, std::map<int, double>>& mutation_priors,
        const character_state_matrix& M,
        std::ranlux48_base& gen, 
        bool greedy
    );

}

#endif // STARHOMOPLASY_HPP
