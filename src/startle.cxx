#include <pprint.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <argparse/argparse.hpp>
#include <csv.hpp>

#include "startle.hpp"
#include "digraph.hpp"
#include "treeio.hpp"
#include "starhomoplasy.hpp"

#include <cmath>
#include <unordered_map>
#include <random>
#include <chrono>
#include <functional>
#include <random>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <optional>
#include <stack>
#include <tuple>

using json = nlohmann::json;
using namespace starhomoplasy;

/*
 * Parses a CSV file describing the character-state 
 * matrix $M$ of the form:
 *      ,c0,c1,...,cn
 *    r0,0, -1,...,2
 *    r1,1,  0,...,1
 *    ...
 *    rm,2,  1,...,0
 * where $M_{ij}$ is the character state of character $j$ in cell $i$.
 */
character_state_matrix parse_character_state_matrix(const std::string& filename){
    character_state_matrix M;

    csv::CSVReader reader(filename);
    auto columns = reader.get_col_names();
    for (auto it = columns.begin() + 1; it != columns.end(); ++it) {
        M.characters.push_back(*it);
    }

    for (auto row: reader) {
        M.cells.push_back(row.begin()->get<std::string>());
        std::vector<int> row_values;
        for (auto it = row.begin() + 1; it != row.end(); ++it) {
            row_values.push_back(it->get<int>());
        }
        M.matrix.push_back(row_values);
    }

    if (M.matrix.size() == 0) {
        throw std::runtime_error("character state matrix is empty");
    }

    size_t num_characters = M.matrix[0].size();
    for (auto row: M.matrix) {
        if (row.size() != num_characters) {
            throw std::runtime_error("character state matrix is not rectangular");
        }
    }

    M.num_characters = num_characters;
    M.num_cells = M.matrix.size();
    return M;
}

/* 
 * Parses a CSV file describing the mutation priors of the form:
 *    character,state,probability
 *    c0,0,0.1
 *    c0,1,0.2
 *    ...
 *    cn,0,0.3
 *    cn,1,0.4
 * where $P_{ij}$ is the probability of a mutation from state $j$ to state $i$.
 */
std::map<std::string, std::map<int, double>> parse_mutation_priors(const std::string& filename) {
    std::map<std::string, std::map<int, double>> mutation_priors;
    csv::CSVReader reader(filename);
    for (auto row: reader) {
        mutation_priors[row["character"].get<std::string>()][row["state"].get<int>()] = row["probability"].get<double>();
    }
    return mutation_priors;
}

void solve_small_parsimony(argparse::ArgumentParser small) {
    std::stringstream ss;
    pprint::PrettyPrinter printer(ss);    
    auto console = spdlog::get("startle");

    console->info("reading character-state matrix from {}", small.get<std::string>("character_matrix"));
    const character_state_matrix M = parse_character_state_matrix(small.get<std::string>("character_matrix"));
    console->info("character matrix has {} cells and {} characters", M.num_cells, M.num_characters);

    console->info("reading mutation priors from {}", small.get<std::string>("mutation_priors"));
    std::map<std::string, std::map<int, double>> mutation_priors = parse_mutation_priors(small.get<std::string>("mutation_priors"));

    for (const auto &[character, priors] : mutation_priors) {
        for (const auto &[state, probability] : priors) {
            if (small.get<bool>("unweighted")) {
                mutation_priors[character][state] = 1;
            } else {
                mutation_priors[character][state] = -log(probability);
            }
        }
    }

    console->info("reading tree from {}", small.get<std::string>("tree"));
    std::ifstream in(small.get<std::string>("tree"));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string tree_newick = buffer.str();

    digraph<treeio::newick_vertex_data> raw_tree = treeio::read_newick_node(tree_newick);
    digraph<star_homoplasy_data> tree;

    std::map<std::string, int> cell_to_idx;
    for (size_t i = 0; i < M.num_cells; i++) {
        cell_to_idx[M.cells[i]] = i;
    }

    for (const auto &u : raw_tree.nodes()) {
        if (!raw_tree.out_degree(u)) {
            if (cell_to_idx.find(raw_tree[u].data.name) == cell_to_idx.end()) {
                throw std::runtime_error("cell " + raw_tree[u].data.name + " not found in character-state matrix");
            }

            std::vector<int> character_states(M.num_characters, -1);
            int row_idx = cell_to_idx[raw_tree[u].data.name];
            for (size_t j = 0; j < M.num_characters; j++) {
                character_states[j] = M.matrix[row_idx][j];
            }

            star_homoplasy_data u_data = {
                .character_states = character_states,
                .internal = false,
                .cell = raw_tree[u].data.name,
                .valid = true,
                .parsimony_score = 0
            };

            tree.add_vertex(u_data);
            continue;
        }

        star_homoplasy_data u_data = {
            .character_states = std::vector<int>(M.num_characters, -1),
            .internal = true,
            .cell = raw_tree[u].data.name,
            .valid = false,
            .parsimony_score = 0
        };

        tree.add_vertex(u_data);
    }

    for (const auto &u : raw_tree.nodes()) {
        for (const auto &v : raw_tree.successors(u)) {
            tree.add_edge(u, v);
        }
    }

    console->info("arbitrarily resolving polytomies");
    arbitrarily_resolve_polytomies(tree);

    console->info("solving small parsimony problem");
    small_parsimony(tree, mutation_priors, M, 0, 0);
    console->info("small parsimony score = {}", tree[0].data.parsimony_score);

    json output;
    output["objective_value"] = 0;
    std::ofstream json_output(small.get<std::string>("output") + "_results.json");
    json_output << output.dump(4) << std::endl;

    return;

}

void solve_large_parsimony(argparse::ArgumentParser large) {
    std::stringstream ss;
    pprint::PrettyPrinter printer(ss);    
    auto console = spdlog::get("startle");

    console->info("reading character-state matrix from {}", large.get<std::string>("character_matrix"));
    const character_state_matrix M = parse_character_state_matrix(large.get<std::string>("character_matrix"));
    console->info("character matrix has {} cells and {} characters", M.num_cells, M.num_characters);

    console->info("reading mutation priors from {}", large.get<std::string>("mutation_priors"));
    std::map<std::string, std::map<int, double>> mutation_priors = parse_mutation_priors(large.get<std::string>("mutation_priors"));

    for (const auto &[character, priors] : mutation_priors) {
        for (const auto &[state, probability] : priors) {
            if (large.get<bool>("unweighted")) {
                mutation_priors[character][state] = 1;
            } else {
                mutation_priors[character][state] = -log(probability);
            }
        }
    }

    console->info("reading tree from {}", large.get<std::string>("tree"));
    std::ifstream in(large.get<std::string>("tree"));
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string tree_newick = buffer.str();

    digraph<treeio::newick_vertex_data> raw_tree = treeio::read_newick_node(tree_newick);
    digraph<star_homoplasy_data> tree;

    std::map<std::string, int> cell_to_idx;
    for (size_t i = 0; i < M.num_cells; i++) {
        cell_to_idx[M.cells[i]] = i;
    }

    for (const auto &u : raw_tree.nodes()) {
        if (!raw_tree.out_degree(u)) {
            if (cell_to_idx.find(raw_tree[u].data.name) == cell_to_idx.end()) {
                throw std::runtime_error("cell " + raw_tree[u].data.name + " not found in character-state matrix");
            }

            std::vector<int> character_states(M.num_characters, -1);
            int row_idx = cell_to_idx[raw_tree[u].data.name];
            for (size_t j = 0; j < M.num_characters; j++) {
                character_states[j] = M.matrix[row_idx][j];
            }

            star_homoplasy_data u_data = {
                .character_states = character_states,
                .internal = false,
                .cell = raw_tree[u].data.name,
                .valid = true,
                .parsimony_score = 0
            };

            tree.add_vertex(u_data);
            continue;
        }

        star_homoplasy_data u_data = {
            .character_states = std::vector<int>(M.num_characters, -1),
            .internal = true,
            .cell = raw_tree[u].data.name,
            .valid = false,
            .parsimony_score = 0
        };

        tree.add_vertex(u_data);
    }

    for (const auto &u : raw_tree.nodes()) {
        for (const auto &v : raw_tree.successors(u)) {
            tree.add_edge(u, v);
        }
    }

    console->info("arbitrarily resolving polytomies");
    arbitrarily_resolve_polytomies(tree);

    console->info("solving large parsimony problem");
    std::random_device rd;
    std::ranlux48_base gen(rd());

    /*
      Candidate tree set is obtained by randomly
      perturbing candidate trees.
     */
    std::vector<digraph<star_homoplasy_data>> candidate_trees;
    unsigned int num_candidate_trees = large.get<unsigned int>("--num-candidates");
    for (int i = 0; i < num_candidate_trees; i++) {
        float aggression = 0.2 * i;
        digraph<star_homoplasy_data> t = stochastic_nni(tree, gen, aggression);
        small_parsimony(t, mutation_priors, M, 0, 0);
        candidate_trees.push_back(t);
        spdlog::info("candidate tree (aggression {}) score: {}", aggression, t[0].data.parsimony_score);
    }


    unsigned int num_threads = large.get<unsigned int>("-t");
    if (num_threads > 1) {
        spdlog::info("using {} threads", num_threads);
    }
    
    std::mutex progress_mutex;
    std::vector<std::thread> threads;

    json progress_information;
    std::atomic<size_t> counter = 0, iteration = 0;
    for (unsigned int i = 0; i < num_threads; i++) {
        threads.push_back(std::thread([&]() {
            int thread_id = i;
            std::random_device rd;
            std::ranlux48_base gen(rd());

            while (counter < large.get<size_t>("-i")) {
                int current_iteration = iteration.load();
                iteration++;

                std::vector<double> scores;
                {
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    for (auto& candidate_tree : candidate_trees) {
                        scores.push_back(candidate_tree[0].data.parsimony_score);
                    }

                    json progress_information_i;
                    progress_information_i["scores"] = scores;
                    progress_information_i["iteration"] = current_iteration;
                    progress_information.push_back(progress_information_i);
                }
                
                /* compute summary statistics of candidate tree scores */
                { 
                    double min = *std::min_element(scores.begin(), scores.end());
                    double max = *std::max_element(scores.begin(), scores.end());
                    double median = scores[scores.size() / 2];
                    double mean = std::accumulate(scores.begin(), scores.end(), 0.0) / scores.size();

                    spdlog::info("(thread {}, iteration {}, counter {}) parsimony scores min: {:.2f}, max: {:.2f}, median: {:.2f}, mean: {:.2f}", 
                                 thread_id, current_iteration, counter, min, max, median, mean );
                }

                // Select and perturb candidate tree.
                std::uniform_int_distribution<int> distrib(0, num_candidate_trees - 1);
                int candidate_tree_idx = distrib(gen);

                digraph<star_homoplasy_data> candidate_tree;

                {
                    // copy candidate tree
                    std::lock_guard<std::mutex> lock(progress_mutex);
                    candidate_tree = candidate_trees[candidate_tree_idx];
                }

                std::uniform_real_distribution<double> aggression_distrib(0, large.get<double>("-a"));
                candidate_tree = stochastic_nni(candidate_tree, gen, aggression_distrib(gen)); // seems to be a NO OP

                digraph<star_homoplasy_data> updated_tree = hill_climb(candidate_tree, mutation_priors, M, gen, large.get<bool>("-g"));
                spdlog::info("thread ID {}: updated tree score is {}", thread_id, updated_tree[0].data.parsimony_score);

                std::lock_guard<std::mutex> lock(progress_mutex);
                counter++;
                auto maximum_it = std::max_element(
                        candidate_trees.begin(), candidate_trees.end(), 
                        [](const digraph<star_homoplasy_data> &a, const digraph<star_homoplasy_data> &b) {
                            return a[0].data.parsimony_score < b[0].data.parsimony_score;
                });

                if (updated_tree[0].data.parsimony_score < (*maximum_it)[0].data.parsimony_score) {
                    *maximum_it = updated_tree;
                    spdlog::info("thread ID {}: updated candidate tree set.", thread_id);
                    counter = 0;
                    continue;
                } 
            }
        }));
    }

    for (auto& thread : threads) {
        thread.join();
    }

    for (auto& candidate_tree : candidate_trees) {
        small_parsimony(candidate_tree, mutation_priors, M, 0, 0);
    }

    std::sort(candidate_trees.begin(), candidate_trees.end(),
              [](const digraph<star_homoplasy_data> &a, const digraph<star_homoplasy_data> &b) {
                  return a[0].data.parsimony_score < b[0].data.parsimony_score;
              });

    std::string newick_string = treeio::print_newick_tree(candidate_trees[0]);
    newick_string += ";";

    console->info("writing output files");

    console->info("writing newick tree to {}", large.get<std::string>("-o") + "_tree.newick");
    std::ofstream newick_output(large.get<std::string>("-o") + "_tree.newick", std::ios::out);
    newick_output << newick_string;

    console->info("writing newick tree to {}", large.get<std::string>("-o") + "_tree.json");
    std::ofstream info_output(large.get<std::string>("-o") + "_info.json", std::ios::out);
    info_output << progress_information.dump();

    return;
}


int main(int argc, char *argv[])
{
    auto console_logger = spdlog::stdout_color_mt("startle");
    spdlog::set_default_logger(console_logger);

    auto error_logger = spdlog::stderr_color_mt("error");

    argparse::ArgumentParser program(
        "startle",
        std::to_string(STARTLE_VERSION_MAJOR) + "." + std::to_string(STARTLE_VERSION_MINOR)
    );

    argparse::ArgumentParser small(
        "small"
    );

    argparse::ArgumentParser large(
        "large"
    );

    small.add_description("Solves the star homoplasy SMALL parsimony problem.");
    large.add_description("Solves the star homoplasy LARGE parsimony problem.");

    small.add_argument("character_matrix")
           .help("CSV file containing the character-state matrix");

    small.add_argument("mutation_priors")
           .help("CSV file describing the mutation priors");

    small.add_argument("tree")
            .help("Newick file describing the phylogenetic tree");

    small.add_argument("-o", "--output")
           .help("prefix of the output files")
           .required();

    small.add_argument("--unweighted")
            .help("use unweighted parsimony")
            .default_value(false)
            .implicit_value(true);

    large.add_argument("character_matrix")
           .help("CSV file containing the character-state matrix");

    large.add_argument("mutation_priors")
           .help("CSV file describing the mutation priors");

    large.add_argument("tree")
            .help("Newick file describing the starting (seed) phylogenetic tree");

    large.add_argument("-o", "--output")
          .help("prefix of the output files")
          .required();

    large.add_argument("-i", "--iterations")
          .help("number of iterations to use in the stochastic hill climbing algorithm")
          .default_value((size_t) 400)
          .scan<'u', size_t>();

    large.add_argument("-a", "--aggression")
          .help("aggression of the stochastic hill climbing algorithm")
          .default_value(0.4)
          .scan<'f', double>();

    large.add_argument("-r", "--random-seed") // TODO: get working
          .help("random seed")
          .default_value(0)
          .scan<'d', int>();

    large.add_argument("-t", "--threads")
          .help("number of threads to use")
          .default_value(std::thread::hardware_concurrency())
          .scan<'u', unsigned int>();

    large.add_argument("-g", "--greedy")
          .help("use greedy hill climbing strategy as opposed to full NNI neighborhood exploration")
          .default_value(false)
          .implicit_value(true);

    large.add_argument("--unweighted")
          .help("use unweighted parsimony")
          .default_value(false)
          .implicit_value(true);

    large.add_argument("--num-candidates")
          .help("number of candidate trees to use in the stochastic hill climbing algorithm")
          .default_value(10U)
          .scan<'u', unsigned int>();

    program.add_subparser(small);
    program.add_subparser(large);
    
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;

        if (program.is_subcommand_used(small)) {
            std::cerr << small;
        } else if (program.is_subcommand_used(large)) {
            std::cerr << large;
        } else {
            std::cerr << program;
        }

        std::exit(1);
    }

    if (program.is_subcommand_used(small)) {
        solve_small_parsimony(small);
    } else if (program.is_subcommand_used(large)) {
        solve_large_parsimony(large);
    } else {
        std::cerr << program;
    }

    return 0;
}
