#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <set>
#include <unordered_set>
#include <cmath>
#include <utility>
#include <tuple>
#include <numeric>

#include "graph_utils.hpp" // for get_grid_diameter()
#include "uniform_grid_graph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "error_utils.h"
#include "kernels.h"

extern "C" {
    SEXP S_create_uniform_grid_graph(SEXP s_input_adj_list,
                                     SEXP s_input_weight_list,
                                     SEXP s_grid_size,
                                     SEXP s_start_vertex,
                                     SEXP s_snap_tolerance);
}

double uniform_grid_graph_t::get_parent_graph_edge_weight(
    const edge_weights_t& weights,
    size_t v1, size_t v2
    ) const {

    #if 1
    // Input validation
    if (v1 >= adjacency_list.size() || v2 >= adjacency_list.size()) {
        REPORT_ERROR("Vertex index out of bounds: v1=%zu, v2=%zu, graph size=%zu",
                 v1, v2, adjacency_list.size());
    }
    #endif

    auto v_min = std::min(v1, v2);
    auto v_max = std::max(v1, v2);

    auto it = weights.find({v_min, v_max});
    if (it != weights.end()) {
        return it->second;
    }

    REPORT_WARNING("Warning: Edge (%d,%d) not found in parent graph weights, computing grid distance\n",
                   v_min, v_max);

    // Fallback to computing distance if edge not found
    return find_grid_distance(v1, v2);
}

/**
 * @brief Computes the shortest path distance between two vertices in a uniform grid graph
 *
 * This function performs a localized Dijkstra's algorithm starting from the source vertex
 * and continuing until we reach the target vertex. It's more efficient than a full BFS
 * because it can stop as soon as it finds the target.
 *
 * @param graph The uniform grid graph containing connectivity information
 * @param source Starting vertex for the path
 * @param target Ending vertex for the path
 * @return The shortest path distance between source and target vertices
 */
double uniform_grid_graph_t::find_grid_distance(
    size_t source,
    size_t target
    ) const {

    #if 0
    // Input validation
    if (source >= adjacency_list.size() || target >= adjacency_list.size()) {
        REPORT_ERROR("Warning: Vertex index out of bounds: source=%d, target=%d, graph size=%zu",
                       source, target, adjacency_list.size());
    }
    #endif

    // Early exit if source and target are the same
    if (source == target) return 0.0;

    // Priority queue entries will be pairs of (distance, vertex)
    // Using negative distance because C++'s priority_queue is a max heap
    using queue_entry = std::pair<double, size_t>;
    std::priority_queue<queue_entry, std::vector<queue_entry>, std::greater<>> pq;

    // Track distances to each vertex
    std::vector<double> distances(adjacency_list.size(), INFINITY);
    distances[source] = 0.0;

    pq.push({0.0, source});

    while (!pq.empty()) {
        auto [current_dist, current_vertex] = pq.top();
        pq.pop();

        // If we've reached our target, we're done
        if (current_vertex == target) {
            return current_dist;
        }

        // If we've found a longer path to current_vertex, skip it
        if (current_dist > distances[current_vertex]) {
            continue;
        }

        // Explore all neighbors of the current vertex
        for (const auto& [neighbor, edge_weight] : adjacency_list[current_vertex]) {
            double new_distance = current_dist + edge_weight;

            // If we've found a shorter path to the neighbor
            if (new_distance < distances[neighbor]) {
                distances[neighbor] = new_distance;
                pq.push({new_distance, neighbor});
            }
        }
    }

    // If we get here, there's no path from source to target
    REPORT_WARNING("Warning: No path found between vertices %d and %d\n", source, target);

    return INFINITY;
}

/**
 * @brief Precomputes and caches edge weights for efficient graph operations
 *
 * This function processes the graph's adjacency list and corresponding weight list
 * to create a cache of edge weights optimized for fast lookup. It stores edge weights
 * in a memory-efficient sparse representation using an unordered_map with custom hashing
 * for integer pairs.
 *
 * To minimize memory usage while maintaining O(1) lookup time, each edge weight is stored
 * exactly once, taking advantage of the graph's undirected nature. For any edge between
 * vertices i and j, the weight is stored only when i < j. This halves the storage
 * requirement while maintaining fast access through consistent key ordering.
 *
 * Performance Characteristics:
 *   - Time Complexity: O(E) where E is the number of edges
 *   - Space Complexity: O(E) for storing weights of existing edges
 *   - Hash map is pre-reserved to avoid rehashing during construction
 *
 * Memory Usage Example:
 *   For a graph with V vertices and approximately 6V edges:
 *   - V = 10K vertices:    ~2.8 MB
 *   - V = 100K vertices:   ~28.8 MB
 *   - V = 1M vertices:     ~288 MB
 *
 * @param adj_list Vector of adjacency lists where adj_list[i] contains all vertices
 *                 connected to vertex i in the original graph
 * @param weight_list Vector of weight lists where weight_list[i][j] contains the weight
 *                    of the edge between vertex i and the j-th vertex in adj_list[i]
 *
 * @return An unordered_map containing edge weights, keyed by vertex pairs
 *
 * @note Edge weights are stored only for vertex pairs (i,j) where i < j
 * @note The function includes bounds checking to ensure weight_list indices are valid
 *
 * Example Usage:
 * @code
 *     std::vector<std::vector<size_t>> adj_list = {{1,2}, {0,2}, {0,1}};
 *     std::vector<std::vector<double>> weights = {{1.0,2.0}, {1.0,3.0}, {2.0,3.0}};
 *     auto edge_weights = precompute_edge_weights(adj_list, weights);
 *     // Access weight between vertices 0 and 1:
 *     double weight_0_1 = edge_weights[{0,1}];
 * @endcode
 *
 * @throw Rf_error if weight_list indices are out of bounds
 *
 * @see edge_weights_t For the type definition of the returned map
 * @see int_pair_hash_t For the hash function used by the unordered_map
 */
edge_weights_t precompute_edge_weights(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list
) {
    edge_weights_t weights;

    // Reserve space to avoid rehashing
    size_t estimated_edges = 0;
    for (const auto& adj : adj_list) {
        estimated_edges += adj.size();
    }
    estimated_edges /= 2;  // Divide by 2 since we store unique edges (i < j)
    weights.reserve(estimated_edges);

    // Go through adjacency list to compute weights
    for (size_t i = 0; i < adj_list.size(); ++i) {
        for (size_t j_idx = 0; j_idx < adj_list[i].size(); ++j_idx) {
            size_t j = static_cast<size_t>(adj_list[i][j_idx]);
            if (i < j) {
                // Bounds checking for extra safety
                //if (i < weight_list.size() && j_idx < weight_list[i].size()) {
                weights[{i, j}] = weight_list[i][j_idx];
                // } else {
                //     Rf_error("Weight list index out of bounds: i=%d, j_idx=%d",
                //             (int)i, (int)j_idx);
                //}
            }
        }
    }

    return weights;
}


/**
 * @struct vertex_data_t
 * @brief Tracks metadata for each vertex during graph construction
 *
 * Maintains information about vertex distances, visitation status, and
 * parent relationships in the shortest path tree from the start vertex.
 */
struct vertex_data_t {
    double distance_from_start;  ///< Distance from the start vertex along shortest path
    bool visited;               ///< Whether vertex has been visited in BFS
    int parent;                 ///< Parent vertex in the shortest path tree

    /// Default constructor initializing a vertex with no parent
    vertex_data_t() : distance_from_start(0.0), visited(false), parent(-1) {}
};

/**
 * @struct edge_hash_t
 * @brief Hash function for std::pair<int,int> representing edges
 *
 * Enables the use of unordered containers with edge pairs as keys.
 */
struct edge_hash_t {
    /**
     * @brief Hash operator for edge pairs
     * @param edge The edge pair to hash
     * @return Size_t hash value for the edge
     */
    std::size_t operator()(const std::pair<int, int>& edge) const {
        std::size_t h1 = std::hash<int>{}(edge.first);
        std::size_t h2 = std::hash<int>{}(edge.second);
        return h1 ^ (h2 << 1);
    }
};

// detailed distance tracking
void print_vertex_distances(const std::vector<vertex_data_t>& vertex_data) {
    for (size_t i = 0; i < vertex_data.size(); ++i) {
        if (vertex_data[i].visited) {
            Rprintf("Vertex %zu: distance=%.4f, parent=%d\n",
                    i, vertex_data[i].distance_from_start, vertex_data[i].parent);
        }
    }
}

// grid point placement verification
void verify_grid_point_placement(const uniform_grid_graph_t& graph,
                               const std::vector<vertex_data_t>& vertex_data,
                               double grid_spacing) {
    for (int v : graph.grid_vertices) {
        double expected_position = std::round(vertex_data[v].distance_from_start / grid_spacing) * grid_spacing;
        double actual_position = vertex_data[v].distance_from_start;
        Rprintf("Grid point %d: expected=%.4f, actual=%.4f, diff=%.4f\n",
                v, expected_position, actual_position,
                std::abs(expected_position - actual_position));
    }
}

// Track edge processing order
void log_edge_processing(
    int v1,
    int v2,
    double d_start,
    double d_end,
    double next_grid_pos
    ) {
    Rprintf("Processing edge (%d,%d): d_start=%.4f, d_end=%.4f, "
            "next_grid_pos=%.4f\n",
            v1, v2, d_start, d_end, next_grid_pos);
}


void uniform_grid_graph_t::print_grid_vertices(size_t shift) const {
    Rprintf("Grid vertices: ");
    if (!grid_vertices.empty()) {
        auto it = grid_vertices.begin();
        Rprintf("%zu", *it + shift);
        ++it;
        for (; it != grid_vertices.end(); ++it) {
            Rprintf(", %zu", *it + shift);
        }
    }
    Rprintf("\n\n");
}


/**
 * @brief Prints the graph structure for debugging
 *
 * @param split If true, prints adjacency and weight lists separately
 * @param shift Optional vertex ID offset for printing
 *
 * Provides detailed output including:
 * - Number of original and total vertices
 * - Set of grid vertices
 * - Complete edge list or separate adjacency/weight lists
 * - Edge weights with high precision
 *
 * The shift parameter is useful when vertex IDs need to be adjusted for display
 * (e.g., for 1-based indexing in R).
 */
void uniform_grid_graph_t::print(
    const std::string& name,
    bool split,
    size_t shift
    ) const {

    if (!name.empty()) {
        Rprintf("\n\nuniform_grid_graph_t: '%s'\n", name.c_str());
    } else {
        Rprintf("\n");
    }

    Rprintf("\nGraph Structure:\n");
    Rprintf("Original vertices: %zu\n", n_original_vertices);
    Rprintf("Total vertices: %zu\n", adjacency_list.size());
    Rprintf("Grid vertices: ");
    for (size_t v : grid_vertices) {
        Rprintf("%zu ", v + shift);
    }
    Rprintf("\n\n");

    if (split) {
        std::vector<std::vector<size_t>> adj_list(adjacency_list.size());
        std::vector<std::vector<double>> weight_list(adjacency_list.size());

        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                adj_list[i].push_back(edge.vertex + shift);
                weight_list[i].push_back(edge.weight);
            }
        }

        Rprintf("Adjacency Lists:\n");
        for (size_t i = 0; i < adj_list.size(); ++i) {
            Rprintf("%zu: ", i + shift);
            for (size_t j = 0; j < adj_list[i].size(); ++j) {
                Rprintf("(%zu, %.6f) ", adj_list[i][j], weight_list[i][j]);
            }
            Rprintf("\n");
        }
    } else {
        Rprintf("Edge List:\n");
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                if (i < edge.vertex) { // Print each edge only once
                    Rprintf("(%zu, %zu): %.6f\n", i + shift, edge.vertex + shift, edge.weight);
                }
            }
        }
    }
}

/**
 * @brief Adds an edge between two vertices in the graph
 *
 * @param graph The graph to modify
 * @param v1 First vertex of the edge
 * @param v2 Second vertex of the edge
 * @param weight Weight (length) of the edge
 *
 * Adds the edge in both directions (undirected graph).
 * Prints debug information about the added edge.
 */
void uniform_grid_graph_t::add_edge(size_t v1, size_t v2, double weight) {
    // Check bounds
    size_t max_vertex = std::max(v1, v2);
    if (max_vertex >= adjacency_list.size()) {
        // This should not happen with the fixed create_uniform_grid_graph,
        // but we keep it as a safety check
        REPORT_ERROR("Vertex index out of bounds (v1=%zu, v2=%zu, adjacency_list.size()=%zu)\n",
                     v1, v2, adjacency_list.size());
    }

    adjacency_list[v1].insert({v2, weight});
    adjacency_list[v2].insert({v1, weight});
}

/**
 * @brief Removes an edge between two vertices in the graph
 *
 * @param graph The graph to modify
 * @param v1 First vertex of the edge
 * @param v2 Second vertex of the edge
 *
 * Removes the edge in both directions (undirected graph).
 * Prints debug information about the removed edge.
 */
void uniform_grid_graph_t::remove_edge(size_t v1, size_t v2) {

    // Create a dummy edge_info_t with any weight (0.0) since our comparator only checks vertices
    edge_info_t dummy{v2, 0.0};
    // Use set's efficient find operation
    auto it1 = adjacency_list[v1].find(dummy);
    if (it1 != adjacency_list[v1].end()) {
        adjacency_list[v1].erase(it1);

        // Do the same for the reverse edge
        edge_info_t reverse_dummy{v1, 0.0};
        adjacency_list[v2].erase(reverse_dummy);

        //Rprintf("DEBUG: Removed edge (%d, %d) with weight %.6f\n", v1, v2, weight);
    }
}

/**
 * @brief Computes the spacing between grid points
 *
 * @param adj_list Input graph adjacency list
 * @param weight_list Input graph edge weights
 * @param grid_size Desired number of grid points
 * @return double The computed grid spacing
 * @throw Runtime error if grid_size <= 1
 *
 * Calculates the uniform spacing between grid points based on the total
 * graph length and desired number of grid points. The spacing is computed
 * as total_length / (grid_size - 1).
 */
double compute_grid_spacing(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int grid_size) {

    #define DEBUG__compute_grid_spacing 0

    if (grid_size <= 1) {
        error("grid_size must be greater than 1");
    }

    double total_length = 0.0;
    std::unordered_set<std::pair<int,int>, edge_hash_t> counted_edges;

    // I wonder if we could much faster record discovered edges by using two vectors
    // std::vector<bool> first_vertex(adj_list.size(), false);
    // std::vector<bool> second_vertex(adj_list.size(), false);
    // each edge is uniquely identified by two size_t elements (first_vertex_index, second_vertex_index), where first_vertex_index < second_vertex_index
    // and so when we discover an edge (first_vertex_index, second_vertex_index), where first_vertex_index < second_vertex_index
    // we set
    // first_vertex[first_vertex_index] = true;
    // second_vertex[second_vertex_index] = true;
    // thus a vertex has been seen/detected if and only if  (first_vertex[first_vertex_index] && second_vertex[first_vertex_index])

    for (size_t i = 0; i < adj_list.size(); ++i) {
        for (size_t j = 0; j < adj_list[i].size(); ++j) {
            int neighbor = adj_list[i][j];
            auto edge = std::make_pair(std::min(static_cast<int>(i), neighbor),
                                       std::max(static_cast<int>(i), neighbor));

            if (counted_edges.insert(edge).second) {
                total_length += weight_list[i][j];
            }
        }
    }

    double grid_spacing = total_length / (grid_size - 1);

    #if DEBUG__compute_grid_spacing
    Rprintf("Total graph length: %.6f\n", total_length);
    Rprintf("Computed grid spacing: %.6f\n", grid_spacing);
    #endif

    return grid_spacing;
}

/**
 * @brief Creates a uniform grid graph from an input graph
 *
 * @param input_adj_list Input graph adjacency list
 * @param input_weight_list Input graph edge weights
 * @param grid_size Desired number of grid points
 * @param start_vertex Starting vertex for grid point placement
 * @param snap_tolerance Tolerance factor for snapping vertices to grid positions
 *
 * @return uniform_grid_graph_t The resulting uniform grid graph
 *
 * @details
 * This function creates a new graph where vertices are placed at uniform intervals
 * along paths from the start vertex. The algorithm:
 * 1. Computes appropriate grid spacing
 * 2. Uses BFS to process edges in order of distance from start
 * 3. Places grid points at multiples of grid_spacing from start
 * 4. Subdivides edges as needed to maintain uniform spacing
 * 5. Snaps vertices to grid points if they are very close to grid positions
 *
 * Edge cases handled:
 * - Vertices very close to grid positions (within EPSILON)
 * - Edges shorter than grid spacing
 * - Final vertex placement
 *
 * The resulting graph maintains the following properties:
 * - All grid points are at exact multiples of grid_spacing from start
 * - Original vertices may become grid points if they are close to grid positions
 * - Edge weights represent geometric distances
 */
uniform_grid_graph_t create_uniform_grid_graph(
    const std::vector<std::vector<int>>& input_adj_list,
    const std::vector<std::vector<double>>& input_weight_list,
    size_t grid_size,
    size_t start_vertex,
    double snap_tolerance) {

    double graph_diameter = get_vertex_eccentricity(input_adj_list, input_weight_list, start_vertex);
    double grid_spacing = graph_diameter / grid_size;
    double epsilon = grid_spacing * snap_tolerance;

    #define DEBUG__create_uniform_grid_graph 0
    #define DEBUG__create_uniform_grid_graph2 0

    #if DEBUG__create_uniform_grid_graph2
    Rprintf("\nInitial parameters:\n");
    Rprintf("Grid spacing: %.4f\n", grid_spacing);
    Rprintf("Snap tolerance: %f\n", snap_tolerance);
    Rprintf("Epsilon: %f\n", epsilon);
    #endif

    // Initialize graph structure
    uniform_grid_graph_t expanded_graph;
    expanded_graph.n_original_vertices = input_adj_list.size();

    // Start with a reasonable initial size, but prepare for dynamic growth
    // We allocate space for original vertices plus twice the grid_size as a safety margin
    size_t initial_capacity = input_adj_list.size() + 2 * grid_size;
    expanded_graph.adjacency_list.resize(initial_capacity);

    // Initialize vertex data tracking with the same initial capacity
    std::vector<vertex_data_t> vertex_data(initial_capacity);

    // Convert input graph to the uniform_grid_graph_t format
    for (size_t i = 0; i < input_adj_list.size(); ++i) {
        for (size_t j = 0; j < input_adj_list[i].size(); ++j) {
            size_t vertex = static_cast<size_t>(input_adj_list[i][j]);
            expanded_graph.add_edge(i, vertex, input_weight_list[i][j]);
        }
    }

    // Start with the initial vertex
    expanded_graph.grid_vertices.insert(start_vertex);
    vertex_data[start_vertex].visited = true;

    // BFS queue for processing vertices
    std::queue<int> bfs_queue;
    bfs_queue.push(start_vertex);

    int next_vertex_id = input_adj_list.size();

    // Phase 1: Complete BFS Discovery
    std::vector<std::tuple<int, int, double>> all_edges_to_process;
    while (!bfs_queue.empty()) {
        int current_vertex = bfs_queue.front();
        bfs_queue.pop();

        for (const auto& edge_info : expanded_graph.adjacency_list[current_vertex]) {
            // Collect ALL edges before processing any
            int neighbor = edge_info.vertex;
            double edge_length = edge_info.weight;

            if (!vertex_data[neighbor].visited) {
                // ... collect edge information ...
                vertex_data[neighbor].visited = true;
                vertex_data[neighbor].parent = current_vertex;
                vertex_data[neighbor].distance_from_start =
                    vertex_data[current_vertex].distance_from_start + edge_length;

                bfs_queue.push(neighbor);

                all_edges_to_process.emplace_back(current_vertex, neighbor, edge_length);
            }
        }
    }

    #if DEBUG__create_uniform_grid_graph
    Rprintf("\nVertex distances after BFS discovery:\n");
    print_vertex_distances(vertex_data);  // First placement
    #endif

    // Phase 2: Grid Point Placement
    for (const auto& [v1, v2, edge_length] : all_edges_to_process) {
        double d_start = vertex_data[v1].distance_from_start;
        double d_end = d_start + edge_length;

        // Calculate next grid point position
        int num_complete_intervals = static_cast<int>(d_start / grid_spacing);
        double next_grid_pos = (num_complete_intervals + 1) * grid_spacing;

        #if DEBUG__create_uniform_grid_graph
        log_edge_processing(v1, v2, d_start, d_end, next_grid_pos);
        #endif

        // Check if current vertex should be a grid point
        if (std::abs(next_grid_pos - d_start) < epsilon) {
            expanded_graph.grid_vertices.insert(v1);
            next_grid_pos += grid_spacing;
        }

        if (next_grid_pos < d_end - epsilon) {
            // Remove original edge
            expanded_graph.remove_edge(v1, v2);

            // Add grid points and edges
            int prev_vertex = v1;
            double curr_pos = next_grid_pos;

            while (curr_pos < d_end - epsilon) {
                // Check if we need to expand the graph
                if (next_vertex_id >= expanded_graph.adjacency_list.size()) {
                    size_t new_size = expanded_graph.adjacency_list.size() * 2;
                    expanded_graph.adjacency_list.resize(new_size);
                    vertex_data.resize(new_size);
                }

                int new_vertex = next_vertex_id++;
                expanded_graph.grid_vertices.insert(new_vertex);

                double edge_weight = (prev_vertex == v1) ?
                    curr_pos - d_start : grid_spacing;

                expanded_graph.add_edge(prev_vertex, new_vertex, edge_weight);

                vertex_data[new_vertex].distance_from_start = curr_pos;
                vertex_data[new_vertex].parent = prev_vertex;

                prev_vertex = new_vertex;
                curr_pos += grid_spacing;
            }

            // Connect to final vertex
            double final_weight = (prev_vertex == v1) ?
                edge_length : (d_end - vertex_data[prev_vertex].distance_from_start);

            expanded_graph.add_edge(prev_vertex, v2, final_weight);

            // Check if v2 should be a grid point
            double closest_grid_pos = std::round(d_end / grid_spacing) * grid_spacing;
            if (std::abs(d_end - closest_grid_pos) < epsilon) {
                expanded_graph.grid_vertices.insert(v2);
            }

            #if DEBUG__create_uniform_grid_graph
            Rprintf("\nGrid point verification after processing edge (%d, %d):\n", v1, v2);
            verify_grid_point_placement(expanded_graph, vertex_data, grid_spacing);
            Rprintf("\n");
            #endif

        } else {
            // Edge too short, check if v2 should be a grid point
            double closest_grid_pos = std::round(d_end / grid_spacing) * grid_spacing;
            if (std::abs(d_end - closest_grid_pos) < epsilon) {
                expanded_graph.grid_vertices.insert(v2);
            }
        }
    }

    #if DEBUG__create_uniform_grid_graph
    Rprintf("\nVertex distances after grid point placement:\n");
    print_vertex_distances(vertex_data);  // Second placement
    #endif

    // Trim any unused space
    while (!expanded_graph.adjacency_list.empty() &&
           expanded_graph.adjacency_list.back().empty()) {
        expanded_graph.adjacency_list.pop_back();
    }

    #if DEBUG__create_uniform_grid_graph
    Rprintf("\n");
    expanded_graph.print("In create_uniform_grid_graph(): expanded_graph",false,0);
    #endif

    return expanded_graph;
}

/**
* @brief R interface wrapper for create_uniform_grid_graph()
*
* @param s_adj_list R list representing graph adjacency list (1-based indices)
* @param s_weight_list R list of edge weights corresponding to adjacency list
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
SEXP S_create_uniform_grid_graph(SEXP s_adj_list,
                                 SEXP s_weight_list,
                                 SEXP s_grid_size,
                                 SEXP s_start_vertex,
                                 SEXP s_snap_tolerance) {

    // Convert R inputs to C++ types
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    int grid_size = INTEGER(s_grid_size)[0];
    int start_vertex = INTEGER(s_start_vertex)[0];
    double snap_tolerance = REAL(s_snap_tolerance)[0];

    uniform_grid_graph_t result = create_uniform_grid_graph(
        adj_list,
        weight_list,
        grid_size,
        start_vertex,
        snap_tolerance);

    // Creating return list
    int n_protected = 0;
    const int N_COMPONENTS = 3;
    SEXP r_result = PROTECT(allocVector(VECSXP, N_COMPONENTS)); n_protected++;

    // Extract adjacency and weight lists from result
    int n_total_vertices = result.adjacency_list.size();
    SEXP r_adj_list = PROTECT(allocVector(VECSXP, n_total_vertices)); n_protected++;
    SEXP r_weight_list = PROTECT(allocVector(VECSXP, n_total_vertices)); n_protected++;

    // Convert the set-based representation back to R lists
    for (int i = 0; i < n_total_vertices; ++i) {
        const auto& neighbors = result.adjacency_list[i];

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



/**
 * @brief Computes shortest path distances from a start vertex to all grid vertices
 *
 * @param start_vertex Starting vertex for path computation
 * @return std::vector<double> Vector containing shortest path distances to all grid vertices
 *        (infinity if no path exists)
 *
 * @throws std::out_of_range if start_vertex is out of range
 * @note Time complexity: O((V' + E') log V') where V' and E' are the number of vertices and edges
 *       explored before finding all grid vertices
 */
std::vector<double> uniform_grid_graph_t::compute_all_shortest_path_distances(size_t start_vertex) const {
    if (start_vertex >= adjacency_list.size()) {
        REPORT_ERROR("Start vertex out of range");
    }

    // Initialize distances only for vertices we've seen
    std::unordered_map<size_t, double> dist;
    dist[start_vertex] = 0;

    // Track remaining grid vertices to find
    std::unordered_set<size_t> remaining_targets = grid_vertices;
    if (remaining_targets.empty()) {
        return {};
    }

    // If start vertex is one of the targets, remove it immediately
    if (remaining_targets.count(start_vertex)) {
        remaining_targets.erase(start_vertex);
    }

    // Min-heap priority queue: pairs of (distance, vertex)
    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, start_vertex});

    while (!pq.empty() && !remaining_targets.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        d = -d;  // Convert back to positive distance

        // Skip if we already found a better path
        if (dist.count(u) && d > dist[u]) continue;

        // Process all neighbors
        for (const auto& edge : adjacency_list[u]) {
            size_t v = edge.vertex;
            double w = edge.weight;

            double new_dist = dist[u] + w;
            if (!dist.count(v) || new_dist < dist[v]) {
                dist[v] = new_dist;
                pq.push({-new_dist, v});

                // If this is a grid vertex, remove it from remaining targets
                if (remaining_targets.count(v)) {
                    remaining_targets.erase(v);
                }
            }
        }
    }

    // Create result vector for all grid vertices
    std::vector<double> result;
    result.reserve(grid_vertices.size());
    for (size_t v : grid_vertices) {
        result.push_back(dist.count(v) ? dist[v] : std::numeric_limits<double>::infinity());
    }

    return result;
}

/**
 * @brief Finds shortest paths from a start vertex to all vertices within a specified radius
 *
 * Implements Dijkstra's algorithm with a distance constraint to find paths
 * from the start vertex to all reachable vertices within the given radius.
 * For each reachable vertex, the function stores both the distance from the
 * start vertex and the predecessor vertex in the shortest path.
 *
 * @param start Starting vertex index
 * @param radius Maximum distance from start vertex to consider
 * @return An unordered map where keys are destination vertices and values are pairs of
 *         (distance from start, predecessor vertex in the path)
 * @see compute_reachability_map
 */
std::unordered_map<size_t, std::pair<double, size_t>> uniform_grid_graph_t::find_grid_paths_within_radius(
    size_t start,
    double radius) const {
    std::unordered_map<size_t, std::pair<double, size_t>> result;
    std::vector<double> dist(adjacency_list.size(), INFINITY);
    std::vector<size_t> prev(adjacency_list.size(), INVALID_VERTEX);
    dist[start] = 0;

    std::priority_queue<std::pair<double, size_t>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        d = -d;
        pq.pop();

        if (d > radius) break;
        if (d > dist[u]) continue;

        result.emplace(u, std::pair<double, size_t>(d, prev[u]));

        for (const auto& [v, w] : adjacency_list[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                prev[v] = u;
                pq.push({-dist[v], v});
            }
        }
    }

    return result;
}


/**
 * @brief Computes a reachability map for all vertices within a specified radius
 *
 * Creates a reachability map structure containing distances and predecessors for
 * all vertices reachable from the reference vertex within the given radius.
 * Also populates a list of original vertices sorted by distance (furthest first),
 * excluding the reference vertex itself.
 *
 * The resulting map can be used for spatial analysis and path reconstruction.
 *
 * @param ref_vertex Reference vertex from which to compute reachability
 * @param radius Maximum distance from reference vertex to consider
 * @return A reachability_map_t structure containing sorted vertices, distances, predecessors,
 *         and the reference vertex
 */
reachability_map_t uniform_grid_graph_t::compute_reachability_map(
    size_t ref_vertex,
    double radius) const {

    reachability_map_t result;
    result.ref_vertex = ref_vertex;

    auto shortest_paths_map = find_grid_paths_within_radius(ref_vertex, radius);

    for (const auto& [vertex, dist_parent_info] : shortest_paths_map) {
        result.distances[vertex] = dist_parent_info.first;
        result.predecessors[vertex] = dist_parent_info.second;
        if (vertex != ref_vertex &&
            is_original_vertex(vertex)) {
            result.sorted_vertices.push_back({vertex, dist_parent_info.first});
        }
    }

    std::sort(result.sorted_vertices.begin(), result.sorted_vertices.end(),
              [](const vertex_info_t& a, const vertex_info_t& b) {
                  return a.distance > b.distance;
              });

    return result;
}

/**
 * @brief Finds endpoints of all shortest paths within a radius from a reference vertex
 *
 * An endpoint is defined as the furthest vertex from the reference vertex along
 * a non-overlapping path. This function leverages the path reconstruction logic to
 * identify distinct paths and then extracts their endpoints.
 *
 * @param reachability_map The reachability map computed from a reference vertex
 * @return A vector of endpoint vertices with their distances from the reference vertex
 */
std::vector<vertex_info_t> uniform_grid_graph_t::find_path_endpoints(
    const reachability_map_t& reachability_map) const {

    // Reconstruct non-overlapping paths
    std::vector<grid_vertex_path_t> paths = reconstruct_paths(
        const_cast<reachability_map_t&>(reachability_map));

    // Extract endpoints from paths
    std::vector<vertex_info_t> endpoints;
    endpoints.reserve(paths.size());  // Pre-allocate for efficiency

    for (const auto& path : paths) {
        // Skip empty paths
        if (path.vertices.empty()) {
            continue;
        }

        // Get the endpoint (furthest vertex in the path)
        // After path reversal in reconstruct_paths, the endpoint is the last vertex
        size_t endpoint_vertex = path.vertices.back();
        double endpoint_dist = path.dist_to_ref_vertex.back();

        // Add to endpoints list
        endpoints.push_back({endpoint_vertex, endpoint_dist});
    }

    return endpoints;
}


#if 0
std::vector<xyw_path_t> uniform_grid_graph_t::reconstruct_xyw_paths(
    const reachability_map_t& reachability_map,
    const std::vector<double>& y,
    double radius,
    size_t kernel_type,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    const edge_weights_t& edge_weights) const {

    std::vector<ref_vertex_path_t> paths;
    std::unordered_set<size_t> used_vertices;

    // Process vertices within radius
    for (const auto& vertex_dist_info : reachability_map.sorted_vertices) {
        if (vertex_dist_info.distance > radius) continue;
        if (used_vertices.count(vertex_dist_info.vertex) > 0) continue;

        auto new_path = reconstruct_single_path(
            vertex_dist_info.vertex,
            vertex_dist_info.distance,
            reachability_map
            );

        if (new_path.vertices.size() >= 2) {
            paths.push_back(std::move(new_path));
            for (size_t vertex : new_path.vertices) {
                used_vertices.insert(vertex);
            }
        }
    }

    return process_xyw_paths(
        paths,
        y,
        min_path_size,
        diff_threshold,
        kernel_type,
        dist_normalization_factor,
        edge_weights
        );
}
#endif

// Helper method to reconstruct a single path
ref_vertex_path_t uniform_grid_graph_t::reconstruct_single_path(
    size_t vertex,
    double distance,
    const reachability_map_t& reachability_map) const {

    ref_vertex_path_t path;
    path.target_vertex = reachability_map.ref_vertex;
    path.total_weight = distance;

    size_t curr = vertex;
    while (curr != INVALID_VERTEX) {
        if (is_original_vertex(curr)) {
            path.vertices.push_back(curr);

            auto dist_it = reachability_map.distances.find(curr);
            if (dist_it != reachability_map.distances.end()) {
                path.dist_to_ref_vertex.push_back(dist_it->second);
            }
        }

        auto pred_it = reachability_map.predecessors.find(curr);
        if (pred_it == reachability_map.predecessors.end()) break;
        curr = pred_it->second;
    }

    std::reverse(path.vertices.begin(), path.vertices.end());
    std::reverse(path.dist_to_ref_vertex.begin(), path.dist_to_ref_vertex.end());

    return path;
}

// Helper method to process paths into XYW format
#if 0
std::vector<xyw_path_t> uniform_grid_graph_t::process_xyw_paths(
    const std::vector<ref_vertex_path_t>& paths,
    const std::vector<double>& y,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights) const {

    std::vector<xyw_path_t> result;

    // Handle single path case
    if (paths.size() == 1 && paths[0].vertices.size() >= min_path_size) {
        result.push_back(convert_to_xyw_path(paths[0], y, kernel_type, dist_normalization_factor, edge_weights));
    }
    // Handle multiple paths case
    else if (paths.size() > 1) {
        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = i + 1; j < paths.size(); ++j) {
                if (can_combine_paths(paths[i], paths[j], min_path_size, diff_threshold)) {
                    result.push_back(
                        create_composite_xyw_path(
                            paths[i],
                            paths[j],
                            y,
                            kernel_type,
                            dist_normalization_factor,
                            edge_weights
                            )
                        );
                }
            }
        }
    }

    return result;
}
#endif

// Helper method to check if paths can be combined
bool uniform_grid_graph_t::can_combine_paths(
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    size_t min_path_size,
    size_t diff_threshold) const {

    size_t combined_size = path1.vertices.size() + path2.vertices.size() - 1;
    if (combined_size < min_path_size) return false;

    return paths_explore_different_directions(path1, path2, diff_threshold);
}


// Convert a single reference vertex path to XYW format
#if 0
xyw_path_t uniform_grid_graph_t::convert_to_xyw_path(
    const ref_vertex_path_t& path,
    const std::vector<double>& y,
    size_t kernel_type,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights) const {

    xyw_path_t xyw_path;
    size_t n_vertices = path.vertices.size();

    // Copy vertices and initialize vectors
    xyw_path.vertices = path.vertices;
    xyw_path.x_path.resize(n_vertices);
    xyw_path.y_path.resize(n_vertices);
    xyw_path.w_path.resize(n_vertices);
    // xyw_path.dist_from_ref_vertex = path.dist_to_ref_vertex;

    // Calculate x coordinates (cumulative distances along path)
    xyw_path.x_path[0] = 0.0;
    for (size_t i = 1; i < n_vertices; ++i) {
        double weight = get_parent_graph_edge_weight(
            edge_weights,
            path.vertices[i-1],
            path.vertices[i]
            );
        xyw_path.x_path[i] = xyw_path.x_path[i-1] + weight;
    }

    // Store y values
    for (size_t i = 0; i < n_vertices; ++i) {
        xyw_path.y_path[i] = y[path.vertices[i]];
    }

    // Calculate kernel weights
    compute_kernel_weights(xyw_path, kernel_type, dist_normalization_factor);

    // Store reference vertex position
    xyw_path.x_ref_vertex = xyw_path.x_path[0];  // Since ref_vertex is first in path

    return xyw_path;
}

// Create a composite path by combining two paths
xyw_path_t uniform_grid_graph_t::create_composite_xyw_path(
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    const std::vector<double>& y,
    size_t kernel_type,
    double dist_normalization_factor,
    const edge_weights_t& edge_weights) const {

    xyw_path_t composite_path;

    // Combine vertices: path2 (reversed) + path1 (excluding shared vertex)
    composite_path.vertices.reserve(path1.vertices.size() + path2.vertices.size() - 1);

    // Add reversed path2
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path2.vertices.rbegin(),
        path2.vertices.rend()
        );

    // Add path1 (excluding first vertex if shared)
    size_t offset = (path1.vertices[0] == path2.vertices[0]) ? 1 : 0;
    composite_path.vertices.insert(
        composite_path.vertices.end(),
        path1.vertices.begin() + offset,
        path1.vertices.end()
        );

    // Initialize path vectors
    size_t n_vertices = composite_path.vertices.size();
    composite_path.x_path.resize(n_vertices);
    composite_path.y_path.resize(n_vertices);
    composite_path.w_path.resize(n_vertices);

    // Calculate cumulative distances (x coordinates)
    composite_path.x_path[0] = 0.0;
    for (size_t i = 1; i < n_vertices; ++i) {
        double weight = get_parent_graph_edge_weight(
            edge_weights,
            composite_path.vertices[i-1],
            composite_path.vertices[i]
            );
        composite_path.x_path[i] = composite_path.x_path[i-1] + weight;
    }

    // Store y values
    for (size_t i = 0; i < n_vertices; ++i) {
        composite_path.y_path[i] = y[composite_path.vertices[i]];
    }

    // Calculate distances from reference vertex
    size_t ref_vertex_index = path2.vertices.size() - 1;  // Index where paths join
    composite_path.x_ref_vertex = composite_path.x_path[ref_vertex_index];

    composite_path.dist_from_ref_vertex.resize(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        composite_path.dist_from_ref_vertex[i] =
            std::abs(composite_path.x_path[i] - composite_path.x_ref_vertex);
    }

    // Calculate kernel weights
    compute_kernel_weights(composite_path, kernel_type, dist_normalization_factor);

    return composite_path;
}


// Helper method to compute kernel weights for a path
void uniform_grid_graph_t::compute_kernel_weights(
    xyw_path_t& path,
    size_t kernel_type,
    double dist_normalization_factor) const {

    initialize_kernel(kernel_type, 1.0);

    size_t n_vertices = path.vertices.size();
    std::vector<double> normalized_distances(n_vertices);

    // Find maximum distance for normalization
    double max_dist = 0.0;
    for (size_t i = 0; i < n_vertices; ++i) {
        double dist = std::abs(path.x_path[i] - path.x_ref_vertex);
        max_dist = std::max(max_dist, dist);
    }

    // Normalize distances
    max_dist = std::max(max_dist * dist_normalization_factor, 1.0);
    for (size_t i = 0; i < n_vertices; ++i) {
        normalized_distances[i] =
            std::abs(path.x_path[i] - path.x_ref_vertex) / max_dist;
    }

    // Apply kernel function
    kernel_fn(normalized_distances.data(), n_vertices, path.w_path.data());

    // Normalize weights to sum to 1
    double total_weight = std::accumulate(
        path.w_path.begin(),
        path.w_path.end(),
        0.0
        );

    if (total_weight > 0) {
        for (size_t i = 0; i < n_vertices; ++i) {
            path.w_path[i] /= total_weight;
        }
    }
}
#endif

// Helper method to check if paths explore different directions
bool uniform_grid_graph_t::paths_explore_different_directions(
    const ref_vertex_path_t& path1,
    const ref_vertex_path_t& path2,
    size_t diff_threshold) const {

    std::set<size_t> vertices1, vertices2;

    // Get first diff_threshold vertices from each path (excluding ref vertex)
    for (size_t i = 1; i < path1.vertices.size() && vertices1.size() < diff_threshold; ++i) {
        vertices1.insert(path1.vertices[i]);
    }

    for (size_t i = 1; i < path2.vertices.size() && vertices2.size() < diff_threshold; ++i) {
        vertices2.insert(path2.vertices[i]);
    }

    // Check for intersection
    std::vector<size_t> intersection;
    std::set_intersection(
        vertices1.begin(), vertices1.end(),
        vertices2.begin(), vertices2.end(),
        std::back_inserter(intersection)
        );

    return intersection.empty();
}

#if 0
struct weight_pred_error_t {
    double weight;
    double prediction;
    double error;

    weight_pred_error_t(double w, double p, double e)
        : weight(w), prediction(p), error(e) {}
};


std::vector<path_data_t> generate_paths_for_bandwidth(
    const uniform_grid_graph_t& uniform_grid_graph,
    const reachability_map_t& reachability_map,
    const std::vector<double>& y,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    const edge_weights_t& edge_weights);


std::vector<grid_vertex_path_t> handle_insufficient_path_lengths(
    const std::vector<grid_vertex_path_t>& composite_paths,
    const uniform_grid_graph_t& graph,
    size_t min_path_size,
    size_t grid_vertex,
    double max_bw) {

    // If we have valid paths, return them
    bool has_valid_path = std::any_of(composite_paths.begin(), composite_paths.end(),
        [min_path_size](const grid_vertex_path_t& path) {
            return path.vertices.size() >= min_path_size;
        });

    if (has_valid_path) return composite_paths;

    // Level 1: Try to extend existing paths if possible
    std::vector<grid_vertex_path_t> extended_paths =
        extend_paths(composite_paths, graph, min_path_size);

    if (!extended_paths.empty()) return extended_paths;

    // Level 2: Increase search radius adaptively
    double new_radius = max_bw * 1.5;  // Start with 1.5x increase
    while (new_radius <= graph.estimate_diameter() * 0.5) {  // Limit to half of graph diameter
        auto new_reach_map = graph.compute_reachability_map(grid_vertex, new_radius);
        auto new_paths = graph.reconstruct_paths(new_reach_map);
        auto new_composite_paths = graph.create_composite_paths(new_paths);

        has_valid_path = std::any_of(new_composite_paths.begin(), new_composite_paths.end(),
            [min_path_size](const grid_vertex_path_t& path) {
                return path.vertices.size() >= min_path_size;
            });

        if (has_valid_path) return new_composite_paths;

        new_radius *= 1.5;  // Increase radius by 50% each iteration
    }

    // Level 3: Fall back to relaxed criteria
    return create_fallback_paths(composite_paths, min_path_size);
}

#endif


/**
 * @brief Reconstructs a set of maximal length non-overlapping paths from a reachability map
 *
 * Starting from vertices furthest from the reference vertex, reconstructs paths back to
 * the reference vertex. Each path is constructed to be maximal in length while ensuring
 * no overlap with previously constructed paths (except possibly at the reference vertex
 * if it belongs to the original graph).
 *
 * The paths are reconstructed in order of decreasing distance from the reference vertex,
 * ensuring that longer paths are prioritized. A vertex can belong to at most one path,
 * except for the reference vertex if it's part of the original graph.
 *
 * @param reachability_map The reachability map containing:
 *                        - sorted_vertices: Vertices sorted by distance from reference vertex
 *                        - distances: Map of vertex to distance from reference vertex
 *                        - predecessors: Map of vertex to its predecessor in shortest path
 *                        - ref_vertex: The reference vertex from which paths are constructed
 *
 * @return std::vector<grid_vertex_path_t> A vector of non-overlapping paths, where each path contains:
 *                                        - vertices: Sequence of vertices forming the path
 *                                        - dist_to_ref_vertex: Distance from each vertex to the reference vertex
 *
 * @note Paths are returned with vertices ordered from reference vertex to end vertex
 * @note The reference vertex is included in a path only if it belongs to the original graph
 * @note The function guarantees that returned paths are non-overlapping except possibly
 *       at the reference vertex
 */
std::vector<grid_vertex_path_t> uniform_grid_graph_t::reconstruct_paths(
    reachability_map_t& reachability_map
    ) const {

    std::vector<grid_vertex_path_t> paths;

    // Create a set of vertices to process, ordered in the decreasing order by distance from reference vertex
    // std::set<vertex_info_t, decltype([](const vertex_info_t& a, const vertex_info_t& b) {
    //     return a.distance > b.distance;
    // })> to_process(reachability_map.sorted_vertices.begin(),
    //               reachability_map.sorted_vertices.end());

    struct VertexInfoComparator {
        bool operator()(const vertex_info_t& a, const vertex_info_t& b) const {
            return a.distance > b.distance;
        }
    };

    std::set<vertex_info_t, VertexInfoComparator> to_process(
        reachability_map.sorted_vertices.begin(),
        reachability_map.sorted_vertices.end());

    // Track discovered vertices to ensure non-overlapping paths
    std::vector<bool> discovered(adjacency_list.size(), false);

    while (!to_process.empty()) {
        // Get the furthest unprocessed vertex
        auto current = *to_process.begin();
        to_process.erase(to_process.begin());

        // Skip if this vertex has already been discovered in a previous path
        if (discovered[current.vertex]) {
            continue;
        }

        // Initialize new path
        grid_vertex_path_t path;

        // Always add ref_vertex to the set of grid vertices of the path
        path.grid_vertices.push_back(reachability_map.ref_vertex);
        path.grid_dist_to_ref_vertex.push_back(0);

        size_t vertex = current.vertex;

        // Reconstruct path from current vertex back to reference vertex
        while (vertex != reachability_map.ref_vertex &&
               vertex != INVALID_VERTEX) {
            // Skip if we hit a vertex that's already in another path
            if (vertex != current.vertex && discovered[vertex]) {
                path.vertices.clear();
                path.dist_to_ref_vertex.clear();
                break;
            }

            if (is_original_vertex(vertex)) {
                path.vertices.push_back(vertex);
                path.dist_to_ref_vertex.push_back(reachability_map.distances.at(vertex));
            } else {
                path.grid_vertices.push_back(vertex);
                path.grid_dist_to_ref_vertex.push_back(reachability_map.distances.at(vertex));
            }

            // Mark vertex as discovered
            discovered[vertex] = true;

            // Move to predecessor
            vertex = reachability_map.predecessors.at(vertex);
        }

        // If path reconstruction was successful
        if (!path.vertices.empty()) {
            // Add reference vertex if it's part of the original graph
            if (is_original_vertex(reachability_map.ref_vertex)) {
                path.vertices.push_back(reachability_map.ref_vertex);
                path.dist_to_ref_vertex.push_back(0.0);
            }

            // Reverse the path to start from reference vertex
            std::reverse(path.vertices.begin(), path.vertices.end());
            std::reverse(path.dist_to_ref_vertex.begin(), path.dist_to_ref_vertex.end());

            // Remove vertices from to_process that are now part of this path
            auto it = to_process.begin();
            while (it != to_process.end()) {
                if (discovered[it->vertex]) {
                    it = to_process.erase(it);
                } else {
                    ++it;
                }
            }

            paths.push_back(std::move(path));
        }
    }

    return paths;
}

/**
 * @brief Creates composite paths by combining pairs of grid vertex paths
 *
 * This function generates composite paths by combining input grid vertex paths
 * with the following key characteristics:
 * - Handles both original graph vertices and grid vertices in the paths
 * - Combines paths while avoiding duplicate vertices when paths share a starting point
 * - Preserves distance information for both original and grid vertices
 *
 * Path Combination Logic:
 * - For a single input path, returns it as-is
 * - For multiple paths, creates composite paths by:
 *   1. Reversing the first path of each pair
 *   2. Appending the second path
 *   3. Maintaining relative distances to the reference vertex
 *
 * Complexity:
 * - For n input paths, generates n(n-1)/2 composite paths
 * - Includes both original graph vertices and grid vertices in the composite paths
 *
 * @param grid_vertex_paths Vector of grid vertex paths, each containing
 *        original graph vertices, grid vertices, and their respective distances
 *
 * @return std::vector<compose_path_t> Vector of composite paths created
 *         from combining input paths
 *
 * @note Assumes input paths have consistent vertex and distance vector sizes
 * @note Distances are calculated relative to the reference vertex and path start
 * @note Handles paths with shared starting vertices to avoid duplicates
 *
 * @throws std::runtime_error If input path data is inconsistent
 *         (mismatched vertex and distance vector sizes)
 */
std::vector<compose_path_t> uniform_grid_graph_t::create_composite_paths(
    std::vector<grid_vertex_path_t>& grid_vertex_paths
    ) const {

    std::vector<compose_path_t> composite_paths;

    // Simple Case: Only one path exists
    if (grid_vertex_paths.size() == 1) {
        // Directly convert single grid_vertex_path_t to compose_path_t
        // Ensures no modification of original path while creating composite path
        compose_path_t path;

        path.vertices           = grid_vertex_paths[0].vertices;
        path.dist_to_ref_vertex = grid_vertex_paths[0].dist_to_ref_vertex;
        // Initialize dist_to_from_start as same as dist_to_ref_vertex for single path
        path.dist_to_from_start = grid_vertex_paths[0].dist_to_ref_vertex;

        path.grid_vertices           = grid_vertex_paths[0].grid_vertices;
        path.grid_dist_to_ref_vertex = grid_vertex_paths[0].grid_dist_to_ref_vertex;
        path.grid_dist_to_from_start = grid_vertex_paths[0].grid_dist_to_ref_vertex;

        composite_paths.push_back(path);
        return composite_paths;
    }

    // Comprehensive Path Combination
    for (size_t i = 0; i < grid_vertex_paths.size(); ++i) {
        for (size_t j = i + 1; j < grid_vertex_paths.size(); ++j) {
            // Boundary Condition Check: Validate input paths
            const auto& path1 = grid_vertex_paths[i];
            const auto& path2 = grid_vertex_paths[j];

            // Determine if paths share a starting vertex
            // This helps avoid duplicate vertices in composite path
            size_t offset = (path1.vertices[0] == path2.vertices[0]) ? 1 : 0;

            // Composite Path Construction
            compose_path_t composite_path;

            // Preallocate memory for efficiency
            composite_path.vertices.reserve(path1.vertices.size() + path2.vertices.size() - offset);
            composite_path.dist_to_ref_vertex.reserve(path1.vertices.size() + path2.vertices.size() - offset);

            // Detailed path combination logic
            // 1. Reverse first path (excluding potential shared starting vertex)
            double first_path_length = path1.dist_to_ref_vertex.back();

            #if 0
            // Add original vertices from first path in reverse
            // Restructuring to avoid decrementing beyond zero
            size_t vertices_count = path1.vertices.size();
            for (size_t i = 0; i < vertices_count - offset; ++i) {
                size_t k = vertices_count - 1 - i;  // Calculate the index in reverse order

                composite_path.vertices.push_back(path1.vertices[k]);
                composite_path.dist_to_ref_vertex.push_back(path1.dist_to_ref_vertex[k]);

                // Calculate distance from start of path
                double dist = first_path_length - path1.dist_to_ref_vertex[k];
                composite_path.dist_to_from_start.push_back(dist);
            }

            // Add grid vertices from first path in reverse
            size_t grid_vertices_count = path1.grid_vertices.size();
            if (grid_vertices_count <= offset) {
                // Skip this path or handle appropriately
                continue;
            }
            for (size_t i = 0; i < grid_vertices_count - offset; ++i) {
                size_t k = grid_vertices_count - 1 - i;  // Calculate the index in reverse order

                composite_path.grid_vertices.push_back(path1.grid_vertices[k]);
                composite_path.grid_dist_to_ref_vertex.push_back(path1.grid_dist_to_ref_vertex[k]); // FIXED
                // Calculate distance from start of path
                double dist = first_path_length - path1.grid_dist_to_ref_vertex[k]; // Now consistent
                composite_path.grid_dist_to_from_start.push_back(dist);
            }

            #else

            // Converting to a signed loop counter
            for (int k = static_cast<int>(path1.vertices.size()) - 1; k >= static_cast<int>(offset); --k) {
                composite_path.vertices.push_back(path1.vertices[static_cast<size_t>(k)]);
                composite_path.dist_to_ref_vertex.push_back(path1.dist_to_ref_vertex[static_cast<size_t>(k)]);

                // Calculate distance from start of path
                double dist = first_path_length - path1.dist_to_ref_vertex[static_cast<size_t>(k)];
                composite_path.dist_to_from_start.push_back(dist);
            }

            // Add grid vertices from first path in reverse
            for (int k = static_cast<int>(path1.grid_vertices.size()) - 1; k >= static_cast<int>(offset); --k) {
                composite_path.grid_vertices.push_back(path1.grid_vertices[static_cast<size_t>(k)]);
                composite_path.grid_dist_to_ref_vertex.push_back(path1.grid_dist_to_ref_vertex[static_cast<size_t>(k)]);

                // Calculate distance from start of path
                double dist = first_path_length - path1.grid_dist_to_ref_vertex[static_cast<size_t>(k)];
                composite_path.grid_dist_to_from_start.push_back(dist);
            }
            #endif

            // 2. Add second path forward
            for (size_t k = 0; k < path2.vertices.size(); ++k) {
                composite_path.vertices.push_back(path2.vertices[k]);
                composite_path.dist_to_ref_vertex.push_back(path2.dist_to_ref_vertex[k]);

                // Distance calculation relative to first path's length
                composite_path.dist_to_from_start.push_back(first_path_length + path2.dist_to_ref_vertex[k]);
            }

            // Add grid vertices from second path
            for (size_t k = 0; k < path2.grid_vertices.size(); ++k) {
                composite_path.grid_vertices.push_back(path2.grid_vertices[k]);
                composite_path.grid_dist_to_ref_vertex.push_back(path2.grid_dist_to_ref_vertex[k]);

                // Distance calculation relative to first path's length
                composite_path.grid_dist_to_from_start.push_back(first_path_length + path2.grid_dist_to_ref_vertex[k]);
            }

            // Add validated composite path to results
            composite_paths.push_back(std::move(composite_path));
        }
    }

    return composite_paths;
}


/**
* @brief Checks if there exists a path with at least min_path_size origianal graph vertices and min_num_grid_vertices grid vertices in the given set of paths
*
* @param composite_paths Vector of paths to check
* @param min_path_size Minimum required number of origianl graph vertices in the given path
* @param min_num_grid_vertices Minimum required number of grid graph vertices in the given path
*
* @return true if at least one path has at least min_path_size original graph vertices and at least min_num_grid_vertices grid vertices, false otherwise
*/
bool uniform_grid_graph_t::has_min_size_path(
   const std::vector<compose_path_t>& composite_paths,
   size_t min_path_size,
   size_t min_num_grid_vertices) const {
   return std::any_of(composite_paths.begin(), composite_paths.end(),
       [min_path_size,min_num_grid_vertices](const compose_path_t& path) {
           return (path.vertices.size() >= min_path_size && path.grid_vertices.size() >= min_num_grid_vertices);
       });
}

/**
 * @brief Attempts to find paths that meet minimum size requirement by extending search radius
 *
 * If the current composite paths don't contain any path of minimum required size,
 * this function attempts to find such paths by:
 * 1. Computing the current maximum bandwidth from existing paths
 * 2. Incrementally increasing the search radius and recomputing paths
 * 3. Stopping when either:
 *    - Paths of sufficient size are found
 *    - Search radius exceeds half of the graph diameter
 *
 * @param grid_vertex Original grid vertex from which paths are computed
 * @param min_path_size Minimum required path size
 * @param current_max_bw Current maximum bandwidth used in path computation
 *
 * @return Vector of composite paths, either new paths meeting size requirement or
 *         original paths if new paths couldn't be found
 */
std::vector<compose_path_t> uniform_grid_graph_t::find_min_size_composite_paths(
    size_t grid_vertex,
    size_t min_path_size,
    double current_max_bw,
    double graph_diameter
    ) const {

    // Don't exceed half of the graph's diameter
    double max_allowed_bw = graph_diameter * 0.9;
    // if (current_max_bw >= max_allowed_bw) {
    //     max_allowed_bw = graph_diameter * 0.75;
    // }

    // Start with 1.5 times the current bandwidth
    double new_max_bw = current_max_bw * 1.5;

    while (new_max_bw <= max_allowed_bw) {
        // Compute new reachability map with increased radius
        auto reachability_map = compute_reachability_map(grid_vertex, new_max_bw);

        // Reconstruct base paths
        auto grid_vertex_paths = reconstruct_paths(reachability_map);

        // Create composite paths from base paths
        auto composite_paths = create_composite_paths(grid_vertex_paths);

        // Check if we found paths of sufficient size
        if (has_min_size_path(composite_paths, min_path_size, 2)) {
            return composite_paths;
        }

        // Increase radius for next iteration
        new_max_bw *= 1.5;
    }

    // If we couldn't find paths of sufficient size, log a warning
    REPORT_WARNING("Warning: Could not find paths of size %zu even with extended search radius.\n", min_path_size);

    // Return empty vector to indicate failure
    return std::vector<compose_path_t>();
}


/**
 * @brief Filters composite paths based on bandwidth and minimum size requirements
 *
 * Restricts input composite paths to vertices and grid vertices within a specified
 * distance from the reference vertex, and ensures paths meet minimum size requirements.
 *
 * Filtering Criteria:
 * - Retains only vertices with distance to reference vertex ≤ bandwidth
 * - Requires at least min_path_size original graph vertices per path
 * - Requires at least min_num_grid_vertices grid vertices per path
 * - Maintains all associated distance information
 *
 * @param bw Bandwidth (maximum distance from reference vertex)
 * @param composite_paths Vector of composite paths to filter
 * @param min_path_size Minimum number of original graph vertices required in each path
 * @param min_num_grid_vertices Minimum number of grid vertices required in each path
 *
 * @return std::vector<compose_path_t> Filtered paths meeting all requirements
 *
 * @note Paths failing to meet either bandwidth or size requirements are excluded
 */
std::vector<compose_path_t> uniform_grid_graph_t::get_bw_composite_paths(
    double bw,
    const std::vector<compose_path_t>& composite_paths,
    size_t min_path_size,
    size_t min_num_grid_vertices
) const {
    std::vector<compose_path_t> bw_paths;
    bw_paths.reserve(composite_paths.size());  // Preallocate for efficiency

    for (const auto& path : composite_paths) {
        compose_path_t restricted_path;

        // Process original vertices within bandwidth
        for (size_t i = 0; i < path.vertices.size(); ++i) {
            if (path.dist_to_ref_vertex[i] <= bw) {
                restricted_path.vertices.push_back(path.vertices[i]);
                restricted_path.dist_to_ref_vertex.push_back(path.dist_to_ref_vertex[i]);
                restricted_path.dist_to_from_start.push_back(path.dist_to_from_start[i]);
            }
        }

        // Process grid vertices within bandwidth
        for (size_t i = 0; i < path.grid_vertices.size(); ++i) {
            if (path.grid_dist_to_ref_vertex[i] <= bw) {
                restricted_path.grid_vertices.push_back(path.grid_vertices[i]);
                restricted_path.grid_dist_to_ref_vertex.push_back(path.grid_dist_to_ref_vertex[i]);
                restricted_path.grid_dist_to_from_start.push_back(path.grid_dist_to_from_start[i]);
            }
        }

        // Only keep paths that satisfy both minimum size requirements
        if (restricted_path.vertices.size() >= min_path_size &&
            restricted_path.grid_vertices.size() >= min_num_grid_vertices) {
            bw_paths.push_back(std::move(restricted_path));
        }
    }

    return bw_paths;
}

/**
 * @brief Prepares path data for weighted linear regression
 *
 * Transforms a composite path and response values into a format suitable
 * for local linear regression. This involves:
 *
 * 1. Extracting vertices from the original graph path
 * 2. Computing kernel-weighted distances as regression weights
 * 3. Calculating cumulative distances along the path
 * 4. Selecting corresponding response values
 *
 * Kernel Weighting:
 * - Assigns weights to vertices based on their distance from a reference vertex
 * - Provides more influence to closer vertices in the regression
 *
 * @pre Kernel function must be initialized before calling this method
 * @pre Response vector 'y' must contain values for all vertices in the path
 *
 * @param y Vector of response values for all graph vertices
 * @param path Composite path containing vertex information
 * @param kernel_type Type of kernel function to apply (optional)
 *
 * @return xyw_path_t Prepared data structure for local linear regression
 *
 * @note Assumes consistent sizing between path vertices and response values
 * @warning No input validation is performed
 *
 * @see uniform_grid_graph_t::initialize_kernel() For kernel initialization
 */
xyw_path_t uniform_grid_graph_t::get_xyw_data(
    const std::vector<double>& y,
    compose_path_t& path,
    double dist_normalization_factor
    ) const {

   xyw_path_t result;

   // Copy vertices and distances from reference vertex
   result.vertices      = path.vertices;
   result.grid_vertices = path.grid_vertices;

   int n_vertices = path.vertices.size();
   result.w_path.resize(n_vertices);

   // Normalizing path.dist_to_ref_vertex distances
   double max_dist = 0.0;
   for (int i = 0; i < n_vertices; ++i) {
       max_dist = std::max(max_dist, path.dist_to_ref_vertex[i]);
   }

   if (max_dist == 0) max_dist = 1;
   max_dist *= dist_normalization_factor;

   for (int i = 0; i < n_vertices; ++i) {
       path.dist_to_ref_vertex[i] /= max_dist;
   }

   // it is assumed that initialize_kernel(kernel_type, 1.0); was called by the parent function
   kernel_fn(path.dist_to_ref_vertex.data(), n_vertices, result.w_path.data());
   // <<<--- I am not sure if we really need to normalize result.w_path ???
   double total_w_path = std::accumulate(result.w_path.begin(), result.w_path.end(), 0.0);
   for (int i = 0; i < n_vertices; ++i)
       result.w_path[i] /= total_w_path;


   int n_grid_vertices = path.grid_vertices.size();
   result.grid_w_path.resize(n_grid_vertices);

   // Normalizing path.grid_dist_to_ref_vertex distances
   double max_grid_dist = 0.0;
   for (int i = 0; i < n_grid_vertices; ++i) {
       max_grid_dist = std::max(max_grid_dist, path.grid_dist_to_ref_vertex[i]);
   }

   if (max_grid_dist == 0) max_grid_dist = 1;
   max_grid_dist *= dist_normalization_factor;

   for (int i = 0; i < n_grid_vertices; ++i) {
       path.grid_dist_to_ref_vertex[i] /= max_grid_dist;
   }

   kernel_fn(path.grid_dist_to_ref_vertex.data(), n_grid_vertices, result.grid_w_path.data());
   // <<<--- I am not sure if we really need to normalize result.grid_w_path ???
   double total_grid_w_path = std::accumulate(result.grid_w_path.begin(), result.grid_w_path.end(), 0.0);
   for (int i = 0; i < n_grid_vertices; ++i)
       result.grid_w_path[i] /= total_grid_w_path;

   result.x_path      = path.dist_to_from_start;
   result.grid_x_path = path.grid_dist_to_from_start;

   result.y_path.reserve(n_vertices);
   result.y_path.push_back(y[path.vertices[0]]);
   for (size_t i = 1; i < n_vertices; ++i) {
       result.y_path.push_back(y[path.vertices[i]]);
   }

   return result;
}
