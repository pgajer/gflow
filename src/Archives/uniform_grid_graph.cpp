#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <unordered_set>
#include <queue>
#include <utility>
#include <set>
#include <unordered_set>

#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"


extern "C" {
    SEXP S_create_uniform_grid_graph(SEXP s_input_adj_list,
                                     SEXP s_input_weight_list,
                                     SEXP s_grid_size,
                                     SEXP s_start_vertex);
}


// Hash function for std::pair<int,int>
struct edge_hash_t {
    std::size_t operator()(const std::pair<int,int>& edge) const {
        // Combine hashes of both integers
        // Using bit shifting to mix the values while avoiding overflow
        std::size_t h1 = std::hash<int>{}(edge.first);
        std::size_t h2 = std::hash<int>{}(edge.second);
        return h1 ^ (h2 << 1);  // or use boost::hash_combine
    }
};


/**
 * @brief Prints the structure of a uniform grid graph in either split or combined format
 *
 * @details This function provides two different ways to display a uniform grid graph:
 *          1. Split format (default): Prints separate adjacency and edge weight lists
 *          2. Combined format: Prints vertex-neighbor pairs with weights in a single list
 *
 * In split format, the function generates and prints:
 *   - An adjacency list showing connections between vertices
 *   - A corresponding list of edge weights
 *
 * In combined format, it prints each vertex followed by its neighbors and their
 * associated edge weights in the format: (neighbor_id, weight)
 *
 * @param x        Reference to the uniform grid graph structure to be printed
 * @param split    If true (default), prints separate adjacency and weight lists
 *                 If false, prints combined vertex-neighbor-weight information
 * @param shift    Integer value added to all vertex indices in the output
 *                 Useful for converting between 0-based and 1-based indexing
 *                 Default value is 0
 *
 * @note The function uses Rprintf for output in combined format
 * @note When split is true, it uses helper functions print_vect_vect for output
 *
 * @see print_vect_vect
 * @see uniform_grid_graph_t
 */
void print_uniform_grid_graph(uniform_grid_graph_t& x, bool split = true, int shift = 0) {

    int n_vertices = x.ugg.size();

    if (split) {
        std::vector<std::vector<int>> adj_list(n_vertices);
        std::vector<std::vector<double>> elen_list(n_vertices);

        for (int i = 0; i < n_vertices; ++i) {
            const auto& neighbors = x.ugg[i];
            for (const auto& [neighbor, weight] : neighbors) {
                adj_list[i].push_back(neighbor + shift);
                elen_list[i].push_back(weight);
            }
        }

        print_vect_vect(adj_list, "adj_list", 0, true, 0);
        print_vect_vect(elen_list, "elen_list");

    } else {
        for (int i = 0; i < n_vertices; ++i) {
            const auto& neighbors = x.ugg[i];
            Rprintf("i=%d: ", i + shift);
            for (const auto& [neighbor, weight] : neighbors) {
                Rprintf("(%d, %.2f)", neighbor + shift,weight);
            }
            Rprintf("\n");
        }
    }
}

void remove_edge(uniform_grid_graph_t& x, int v1, int v2) {
    // Create a dummy edge_info_t with any weight (0.0) since our comparator only checks vertices
    edge_info_t dummy_edge{v2, 0.0};
    // Use set's efficient find operation
    auto it1 = x.ugg[v1].find(dummy_edge);
    if (it1 != x.ugg[v1].end()) {
        x.ugg[v1].erase(it1);
        // Do the same for the reverse edge
        edge_info_t reverse_dummy{v1, 0.0};
        auto it2 = x.ugg[v2].find(reverse_dummy);
        if (it2 != x.ugg[v2].end()) {
            x.ugg[v2].erase(it2);
        }
    }
}

void add_edge(uniform_grid_graph_t& x, int v1, int v2, double weight) {
    x.ugg[v1].insert({v2, weight});
    x.ugg[v2].insert({v1, weight});
}


/**
 * @brief Computes grid spacing (eps) to achieve a specified number of grid points
 *
 * Given a graph and desired number of grid points to add, this function calculates
 * the appropriate epsilon value that will result in approximately that many new
 * grid points being added to the graph. This is analogous to creating a uniform
 * grid of specified size in the one-dimensional case.
 *
 * The computation considers the total length of all edges in the graph and
 * distributes the requested number of grid points along these edges. The original
 * vertices of the graph are not counted in the grid_size parameter, just as the
 * endpoints of an interval aren't counted in the grid size for a 1D uniform grid.
 *
 * @param adj_list Adjacency list of the graph
 * @param weight_list Edge weights corresponding to adj_list
 * @param grid_size Number of grid points to add (not counting original vertices)
 * @return Computed epsilon value for grid spacing
 */
double compute_grid_spacing(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int grid_size) {

    if (grid_size <= 0) {
        error("grid_size must be positive");
    }

    // Compute total graph length
    double total_length = 0.0;
    std::unordered_set<std::pair<int,int>, edge_hash_t> counted_edges;

    for (int i = 0; i < adj_list.size(); ++i) {
        for (size_t j = 0; j < adj_list[i].size(); ++j) {
            int neighbor = adj_list[i][j];
            auto edge = std::make_pair(std::min(i, neighbor),
                                     std::max(i, neighbor));

            if (counted_edges.insert(edge).second) {
                total_length += weight_list[i][j];
            }
        }
    }

    // Compute eps to get approximately grid_size new points
    // Note: This is now more directly analogous to the 1D case
    double eps = total_length / (grid_size - 1);

    return eps;
}


/**
* @brief Creates a graph with uniformly spaced vertices along edges of input graph
*
* @param input_adj_list Adjacency list representation of input graph
* @param input_weight_list Edge weights corresponding to input_adj_list
* @param grid_size Target number of grid vertices to create
* @param start_vertex Starting vertex for grid point placement (default: 0)
*
* @return uniform_grid_graph_t Expanded graph with added grid vertices
*
* @details Performs breadth-first traversal from start_vertex, subdividing edges
* into segments of approximately equal length (eps). Grid points are placed at
* segment endpoints. Edge weights determine spacing - total graph length is divided
* by (grid_size-1) to compute eps.
*
* @note Original vertices may or may not be grid points. Grid points are marked in
* returned graph's grid_vertices set.
*
* @throws std::runtime_error if grid_size <= 0
*/
uniform_grid_graph_t create_uniform_grid_graph(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    int grid_size,
    int start_vertex) {

    double eps = compute_grid_spacing(input_adj_list, input_weight_list, grid_size);

    // debugging
    Rprintf("\n\n-------------\nIn create_uniform_grid_graph()\n");
    print_vect_vect(input_adj_list, "input_adj_list");
    print_vect_vect(input_weight_list, "input_weight_list");
    Rprintf("grid_size: %d\tstart_vertex: %d\teps: %f\n", grid_size, start_vertex, eps);

    uniform_grid_graph_t expanded_graph;
    expanded_graph.n_original_vertices = input_adj_list.size();
    expanded_graph.ugg.resize(input_adj_list.size() + grid_size - 2);

    // First, collect all edges from the input graph
    std::vector<std::tuple<int, int, double>> edges_to_process;
    std::unordered_set<std::pair<int,int>, edge_hash_t> processed_edges;

    // Initial BFS to collect edges (only on original vertices)
    std::queue<int> discovery_queue;
    std::vector<bool> visited(input_adj_list.size(), false);

    discovery_queue.push(start_vertex);
    visited[start_vertex] = true;

    // Convert input graph and collect edges
    for (int i = 0; i < input_adj_list.size(); ++i) {
        for (size_t j = 0; j < input_adj_list[i].size(); ++j) {
            expanded_graph.ugg[i].insert({input_adj_list[i][j], input_weight_list[i][j]});
        }
    }

    // debug
    Rprintf("\nexpanded_graph.ugg\n");
    print_uniform_grid_graph(expanded_graph, false);

    // Collect edges in BFS order
    while (!discovery_queue.empty()) {
        int current_vertex = discovery_queue.front();
        discovery_queue.pop();

        for (const auto& edge_info : expanded_graph.ugg[current_vertex]) {
            int neighbor = edge_info.first;
            double edge_length = edge_info.second;

            std::pair<int,int> edge{
                std::min(current_vertex, neighbor),
                std::max(current_vertex, neighbor)
            };

            if (processed_edges.insert(edge).second) {
                edges_to_process.emplace_back(current_vertex, neighbor, edge_length);

                if (!visited[neighbor]) {
                    discovery_queue.push(neighbor);
                    visited[neighbor] = true;
                }
            }
        }
    }

    // Now process the collected edges to add grid vertices
    std::vector<double> distance_to_next_grid(input_adj_list.size(), eps);
    int next_vertex_id = input_adj_list.size();

    expanded_graph.grid_vertices.insert(start_vertex);  // start_vertex is the first element of the uniform grid

    std::vector<int> number_grid_points(input_adj_list.size(), 0); // number_grid_points[v] is the number of grid points on the shortest path from start_vertex to v;  number_grid_points[v] * eps + distance_to_next_grid[v] is the distance from

    // Process edges in the order they were discovered
    const double EPSILON = 1e-10;
    size_t edge_count = 0; // to identify the last vertex of the graph
    size_t index_of_last_edges_to_process = edges_to_process.size() - 1;
    for (const auto& [v1, v2, edge_length] : edges_to_process) {

        bool is_last_edge = (edge_count == index_of_last_edges_to_process);

        double remaining_length = edge_length;

        Rprintf("\n---------\nv1: %d; v2: %d; edge_length: %f; distance_to_next_grid[v1]: %f; remaining_length: %f\n\n",
                v1, v2, edge_length, distance_to_next_grid[v1], remaining_length);


        if (edge_length < EPSILON) {
            Rprintf("Warning: Zero-length edge detected between vertices %d and %d\n", v1, v2);
            continue;
        }

        if (std::abs(remaining_length - distance_to_next_grid[v1]) < EPSILON) {
            // We've hit exactly where a grid point should be

            //debugging
            Rprintf("\n\nWe've hit exactly where a grid point should be\n\n");

            expanded_graph.grid_vertices.insert(v2);
            distance_to_next_grid[v2] = eps;
        } else if (remaining_length < distance_to_next_grid[v1] + EPSILON) {

            //debugging
            Rprintf("\n\nShort edge\n\n");

            distance_to_next_grid[v2] = distance_to_next_grid[v1] - remaining_length;
        } else {

            // debugging
            Rprintf("Remove original edge and add grid points\n");

            // Remove original edge and add grid points
            remove_edge(expanded_graph, v1, v2);

            int prev_vertex = v1;
            int grid_point_number = 0;

            while (remaining_length > eps) {
                int grid_vertex = next_vertex_id++;
                expanded_graph.grid_vertices.insert(grid_vertex);

                // debugging
                Rprintf("grid_point_number: %d\n", grid_point_number);
                Rprintf("adding to expanded_graph: prev_vertex: %d; grid_vertex: %d; edge len: %f\n", prev_vertex, grid_vertex, (grid_point_number == 0 ? distance_to_next_grid[v1] : eps));

                add_edge(expanded_graph, prev_vertex, grid_vertex, (grid_point_number == 0 ? distance_to_next_grid[v1] : eps));

                prev_vertex = grid_vertex;

                remaining_length = edge_length - (grid_point_number * eps + distance_to_next_grid[v1]);

                // debugging
                Rprintf("remaining_length: %f\n", remaining_length);

                if (remaining_length < 0) {

                    // debugging
                    Rprintf("Warning: Negative remaining length detected while processing edge (%d,%d)\n", v1, v2);

                    remaining_length = 0.0;
                }

                if (std::abs(remaining_length - eps) < EPSILON) {
                    // v2 is a grid vertex
                    expanded_graph.grid_vertices.insert(v2);
                    distance_to_next_grid[v2] = eps;
                    add_edge(expanded_graph, prev_vertex, v2, eps);
                    break;
                }

                grid_point_number++;
            }

            // debugging
            Rprintf("v2: %d; expanded_graph.grid_vertices.count(v2): %ld\n", v2, expanded_graph.grid_vertices.count(v2));


            if (expanded_graph.grid_vertices.count(v2) == 0) {
                // Connect to end vertex

                // debugging
                Rprintf("after while adding to expanded_graph: prev_vertex: %d; v2: %d; edge len: %f\n", prev_vertex, v2, remaining_length);

                add_edge(expanded_graph, prev_vertex, v2, remaining_length);
                distance_to_next_grid[v2] = std::max(0.0, eps - remaining_length);
            }
        }

        if (is_last_edge) {
            expanded_graph.grid_vertices.insert(v2);
            distance_to_next_grid[v2] = 0.0;
        }

        ++edge_count;
    }

    while (expanded_graph.ugg.back().empty()) {
        expanded_graph.ugg.pop_back();
    }

    print_uniform_grid_graph(expanded_graph, false);
    //print_uniform_grid_graph(expanded_graph, true);

    return expanded_graph;
}


/**
* @brief R interface wrapper for create_uniform_grid_graph()
*
* @param s_input_adj_list R list representing graph adjacency list (1-based indices)
* @param s_input_weight_list R list of edge weights corresponding to adjacency list
* @param s_grid_size R integer scalar specifying target number of grid vertices
* @param s_start_vertex R integer scalar specifying starting vertex (1-based index)
*
* @return An R list with three named components:
*   - adj_list: List of integer vectors representing expanded graph adjacency list (1-based indices)
*   - weight_list: List of numeric vectors containing edge weights
*   - grid_vertices: Integer vector indicating which vertices are grid points (1-based indices)
*
* @details This function serves as the interface between R and the C++ implementation.
* It handles:
*   1. Converting R data structures to C++ types (adjusting for 1-based indexing)
*   2. Calling the core C++ implementation
*   3. Converting results back to R data structures (restoring 1-based indexing)
*
* The function manages R's garbage collection through PROTECT/UNPROTECT mechanism
* for all allocated R objects.
*
* @seealso create_uniform_grid_graph() for the core implementation
*
* @note Input indices in R are 1-based and are converted to 0-based for C++
* processing. Output indices are converted back to 1-based for R compatibility.
*/
SEXP S_create_uniform_grid_graph(SEXP s_input_adj_list,
                                 SEXP s_input_weight_list,
                                 SEXP s_grid_size,
                                 SEXP s_start_vertex) {

    // Convert R inputs to C++ types
    std::vector<std::vector<int>> input_adj_list = Rgraph_to_vector(s_input_adj_list);
    std::vector<std::vector<double>> input_weight_list = Rweights_to_vector(s_input_weight_list);
    int grid_size = INTEGER(s_grid_size)[0];
    int start_vertex = INTEGER(s_start_vertex)[0];

    uniform_grid_graph_t result = create_uniform_grid_graph(
        input_adj_list,
        input_weight_list,
        grid_size,
        start_vertex);

    // Creating return list
    int n_protected = 0;
    const int N_COMPONENTS = 3;
    SEXP r_result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Extract adjacency and weight lists from result
    int n_total_vertices = result.ugg.size();
    SEXP r_adj_list = PROTECT(allocVector(VECSXP, n_total_vertices)); n_protected++;
    SEXP r_weight_list = PROTECT(allocVector(VECSXP, n_total_vertices)); n_protected++;

    // Convert the set-based representation back to R lists
    for (int i = 0; i < n_total_vertices; ++i) {
        const auto& neighbors = result.ugg[i];

        // Create vectors for this vertex's adjacency list and weights
        SEXP r_adj = PROTECT(allocVector(INTSXP, neighbors.size()));
        SEXP r_weights = PROTECT(allocVector(REALSXP, neighbors.size()));

        // Fill the vectors
        int idx = 0;
        for (const auto& [neighbor, weight] : neighbors) {
            // Convert to 1-based indices for R
            INTEGER(r_adj)[idx] = neighbor + 1;
            REAL(r_weights)[idx] = weight;
            ++idx;
        }

        SET_VECTOR_ELT(r_adj_list, i, r_adj);
        SET_VECTOR_ELT(r_weight_list, i, r_weights);
        UNPROTECT(2); // for r_adj and r_weights
    }

    // Create grid vertices vector (1-based indices)
    int n_grid_vertices = result.grid_vertices.size();
    SEXP r_grid_vertices = PROTECT(allocVector(INTSXP, n_grid_vertices)); n_protected++;

    int counter = 0;
    for (const auto& i : result.grid_vertices) {
        // Convert to 1-based indices for R
        INTEGER(r_grid_vertices)[counter++] = i + 1;
    }

    // Set components in the result list
    SET_VECTOR_ELT(r_result, 0, r_adj_list);
    SET_VECTOR_ELT(r_result, 1, r_weight_list);
    SET_VECTOR_ELT(r_result, 2, r_grid_vertices);

    // Set names for return list
    SEXP names = PROTECT(allocVector(STRSXP, N_COMPONENTS)); n_protected++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("weight_list"));
    SET_STRING_ELT(names, 2, mkChar("grid_vertices"));
    setAttrib(r_result, R_NamesSymbol, names);

    UNPROTECT(n_protected);
    return r_result;
}
