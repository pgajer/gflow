#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"
#include "MS_complex.h"

#include <vector>     // for std::vector
#include <numeric>    // for std::iota
#include <algorithm>  // for std::copy, std::find
#include <set>        // for std::set
#include <map>        // for std::map
#include <utility>    // for std::pair

#include <R.h>
#include <Rinternals.h>

std::vector<int> union_find(const std::vector<std::vector<int>>& adj_vect);
std::unique_ptr<std::unordered_map<int, int>> count_subgraph_components(const std::vector<std::vector<int>>& graph,
                                                                        const std::vector<int>& V);
std::vector<int> shortest_path_by_hops(const std::vector<std::vector<int>>& graph,
                                       int start, int end);


extern "C" {
    SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    SEXP S_graph_MS_cx_with_path_search(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    SEXP S_graph_MS_cx_using_short_h_hops(SEXP s_graph,
                                      SEXP s_hop_list,
                                      SEXP s_core_graph,
                                      SEXP s_Ey);
}

/**
 * @brief Computes connected components of a subgraph induced by a given vertex set
 *
 * @details Given a graph G and a subset of vertices V, this function:
 * 1. Creates a subgraph containing only vertices in V and edges between them
 * 2. Maps original vertex indices to contiguous subgraph indices
 * 3. Computes connected components using union_find
 * 4. Maps component indices back to original vertex indices
 *
 * The returned map groups original vertex indices by their component.
 *
 * Time complexity: O(V'×log(V') + E') where V' is size of input vertex set,
 * E' is number of edges in the induced subgraph
 * Space complexity: O(V' + E')
 *
 * @param graph Adjacency list representation of the full graph
 * @param vertices Set of vertices defining the subgraph
 * @return Map from component ID to set of original vertex indices in that component
 */
std::map<int, std::set<int>> count_subgraph_set_components(
    const std::vector<std::vector<int>>& graph,
    const std::set<int>& vertices) {

    // Create mapping from original vertices to subgraph indices
    std::map<int, int> vertex_to_index;
    std::vector<int> index_to_vertex;
    int idx = 0;
    for (int v : vertices) {
        vertex_to_index[v] = idx++;
        index_to_vertex.push_back(v);
    }

    // Create adjacency list for subgraph
    std::vector<std::vector<int>> subgraph_adj_list(vertices.size());
    for (int v : vertices) {
        int i = vertex_to_index[v];
        for (int neighbor : graph[v]) {
            if (vertices.count(neighbor)) {
                subgraph_adj_list[i].push_back(vertex_to_index[neighbor]);
            }
        }
    }

    // Find components
    auto components = union_find(subgraph_adj_list);

    // Group vertices by component
    std::map<int, std::set<int>> component_groups;
    for (size_t i = 0; i < components.size(); ++i) {
        component_groups[components[i]].insert(index_to_vertex[i]);
    }

    return component_groups;
}

/**
 * @brief Computes the Morse-Smale complex of a function defined on graph vertices
 *
 * @details This function computes gradient flow trajectories and Morse-Smale
 * complex components for a function Ey defined on vertices of a graph G =
 * pIG_k^h(X). For each edge in pIG_k^h(X) the function calls
 * shortest_path_by_hops to find the shortest path between its end points in its
 * core graph pIG_k(X).
 *
 * Key concepts and definitions:
 * 1. Gradient Flow Trajectories:
 *    - For each vertex v, computes ascending and descending trajectories
 *    - Ascending trajectory follows steepest ascent path to local maximum
 *    - Descending trajectory follows steepest descent path to local minimum
 *    - Trajectories are constrained to valid paths in core graph pIG_k(X)
 *
 * 2. Morse-Smale Pro-cells:
 *    - A pro-cell is defined by a pair (M,m) where M is a local maximum and m is a local minimum
 *    - Contains all vertices that lie on gradient trajectories starting at m and ending at M
 *    - Formally: pro-cell(M,m) = {v \in V | exists trajectory through v from m to M}
 *
 * 3. Morse-Smale Cells:
 *    - Connected components of pro-cells after removing their defining extrema
 *    - For pro-cell(M,m), compute components of subgraph induced by:
 *      vertices \in pro-cell(M,m) \ {M,m}
 *    - Unlike classical Morse theory, cells may not be disjoint
 *
 * Implementation Overview:
 * 1. Trajectory Computation:
 *    a. For each vertex v:
 *       - Compute ascending trajectory following steepest ascent
 *       - Compute descending trajectory following steepest descent
 *       - Store endpoints as local extrema
 *       - Update connectivity maps between extrema
 *
 * 2. Pro-cell Construction:
 *    - During trajectory computation, for each vertex v:
 *      - If trajectory connects local minimum m to maximum M
 *      - Add all trajectory vertices to pro-cell(M,m)
 *
 * 3. Cell Decomposition:
 *    - For each pro-cell(M,m):
 *      a. Remove M and m from vertex set
 *      b. Compute connected components of remaining subgraph
 *      c. Each component is a Morse-Smale cell
 *
 * @param s_graph SEXP containing h-th power graph pIG_k^h(X) as adjacency lists
 * @param s_core_graph SEXP containing core graph pIG_k(X) as adjacency lists
 * @param s_Ey SEXP containing function values at vertices
 *
 * @return SEXP List containing:
 *         - trajectories: List of integer vectors with gradient trajectories
 *         - lmax_to_lmin: List mapping local maxima to connected local minima
 *         - lmin_to_lmax: List mapping local minima to connected local maxima
 *         - local_maxima: Integer vector of local maxima indices
 *         - local_minima: Integer vector of local minima indices
 *         - procell_keys: List of (max_idx, min_idx) pairs defining pro-cells
 *         - procells: List of integer vectors containing vertices in each pro-cell
 *         - cells: List of lists, where cells[[i]] contains connected components
 *                 of pro-cell i after removing extrema
 *
 * @throws R Rf_error if input lengths don't match
 */

/**
 * @brief Helper structure to store information about neighboring vertices
 *
 * Used internally by graph_MS_cx to track information about valid neighbors
 * when computing ascending and descending trajectories.
 */
struct neighbor_info_t {
    int vertex;              ///< The neighboring vertex
    std::vector<int> path;   ///< Path to the neighbor
    int path_length;         ///< Length of the path
    double diff_value;       ///< Difference in scalar values between vertices
};


/**
 * @brief Checks if a path is monotonically increasing in scalar values
 *
 * @param path Vector of vertex indices representing the path
 * @param Ey Pointer to array of scalar values
 * @return true if scalar values strictly increase along the path
 * @return false otherwise
 */
bool is_path_monotonic_increasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] <= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_path_monotonic_increasing(const std::vector<int>& path, const std::vector<double>& Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] <= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}


/**
 * @brief Checks if a path is monotonically decreasing in scalar values
 *
 * @param path Vector of vertex indices representing the path
 * @param Ey Pointer to array of scalar values
 * @return true if scalar values strictly decrease along the path
 * @return false otherwise
 */
bool is_path_monotonic_decreasing(const std::vector<int>& path, const double* Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] >= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

bool is_path_monotonic_decreasing(const std::vector<int>& path, const std::vector<double>& Ey) {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        if (Ey[path[i + 1]] >= Ey[path[i]]) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Determines if a vertex is a local maximum
 *
 * @param vertex Index of the vertex to check
 * @param graph Adjacency list representation of the graph
 * @param Ey Pointer to array of scalar values
 * @return true if the vertex is a local maximum
 * @return false otherwise
 */
bool is_local_maximum(int vertex, const std::vector<std::vector<int>>& graph, const double* Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] < Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

bool is_local_maximum(int vertex, const std::vector<std::vector<int>>& graph, const std::vector<double>& Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] < Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Determines if a vertex is a local minimum
 *
 * @param vertex Index of the vertex to check
 * @param graph Adjacency list representation of the graph
 * @param Ey Pointer to array of scalar values
 * @return true if the vertex is a local minimum
 * @return false otherwise
 */
bool is_local_minimum(int vertex, const std::vector<std::vector<int>>& graph, const double* Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] > Ey[neighbor]) {
            return false;
        }
    }
    return true;
}

bool is_local_minimum(int vertex, const std::vector<std::vector<int>>& graph, const std::vector<double>& Ey) {
    for (int neighbor : graph[vertex]) {
        if (Ey[vertex] > Ey[neighbor]) {
            return false;
        }
    }
    return true;
}



SEXP S_graph_MS_cx_with_path_search(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        Rf_error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph      = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<int>> core_graph = convert_adj_list_from_R(s_core_graph);
    double* Ey = REAL(s_Ey);

    // Initialize maps for extrema connectivity
    std::map<int, std::set<int>> lmax_to_lmin_map;
    std::map<int, std::set<int>> lmin_to_lmax_map;
    std::set<int> local_maxima_set;
    std::set<int> local_minima_set;

    // Allocate R objects
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 8));

    // Update names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 8));
        SET_STRING_ELT(names, 0, Rf_mkChar("trajectories"));
        SET_STRING_ELT(names, 1, Rf_mkChar("lmax_to_lmin"));
        SET_STRING_ELT(names, 2, Rf_mkChar("lmin_to_lmax"));
        SET_STRING_ELT(names, 3, Rf_mkChar("local_maxima"));
        SET_STRING_ELT(names, 4, Rf_mkChar("local_minima"));
        SET_STRING_ELT(names, 5, Rf_mkChar("procell_keys"));
        SET_STRING_ELT(names, 6, Rf_mkChar("procells"));
        SET_STRING_ELT(names, 7, Rf_mkChar("cells"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    // Helper function for valid neighbors (same as before)
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : graph[current_vertex]) {

            std::vector<int> path = shortest_path_by_hops(core_graph, current_vertex, neighbor);

            if (path.empty()) {
                Rprintf("\nCRITICAL ERROR: Invalid (graph, core_graph) pair detected!\n");
                Rprintf("Details:\n");
                Rprintf("  - Current vertex: %d\n", current_vertex);
                Rprintf("  - Neighbor vertex: %d\n", neighbor);
                Rprintf("  - Direction: %s\n", ascending ? "ascending" : "descending");
                Rprintf("  - Function values: current=%.6f, neighbor=%.6f\n",
                        Ey[current_vertex], Ey[neighbor]);

                // Print information about vertex connectivity
                Rprintf("\nConnectivity Information:\n");
                Rprintf("  Current vertex %d neighbors in graph: ", current_vertex);
                for (int v : graph[current_vertex]) {
                    Rprintf("%d ", v);
                }
                Rprintf("\n  Current vertex %d neighbors in core_graph: ", current_vertex);
                for (int v : core_graph[current_vertex]) {
                    Rprintf("%d ", v);
                }
                Rprintf("\n  Neighbor vertex %d neighbors in core_graph: ", neighbor);
                for (int v : core_graph[neighbor]) {
                    Rprintf("%d ", v);
                }

                Rf_error("\nFATAL: No path exists in core_graph between vertices that are neighbors in graph.\n"
                      "This violates the fundamental property of (graph, core_graph) pair where\n"
                      "every edge in graph must have a corresponding path in core_graph.\n"
                      "Please verify the construction of your graph and core_graph.");
            }

            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                             is_path_monotonic_increasing(path, Ey) :
                             is_path_monotonic_decreasing(path, Ey))) {
                neighbors.emplace_back(neighbor_info_t{neighbor, path, static_cast<int>(path.size() - 1), diff});
            }
        }

        // Sort by difference value in descending order
        std::sort(neighbors.begin(), neighbors.end(),
                  [](const neighbor_info_t& a, const neighbor_info_t& b) {
                      return a.diff_value > b.diff_value;
                  });

        return neighbors;
    };

    // Add map for pro-cells
    std::map<std::pair<int,int>, std::set<int>> ms_procells;
    std::map<std::pair<int,int>, std::vector<std::set<int>>> ms_cells;

    // 0: trajectories
    {
        SEXP trajectories = PROTECT(Rf_allocVector(VECSXP, n_vertices));

        for (int vertex = 0; vertex < n_vertices; vertex++) {
            bool is_max = is_local_maximum(vertex, graph, Ey);
            bool is_min = is_local_minimum(vertex, graph, Ey);

            if (is_max) {
                local_maxima_set.insert(vertex);
            }
            if (is_min) {
                local_minima_set.insert(vertex);
            }

            // Compute ascending and descending trajectories
            std::vector<int> ascending_trajectory = {vertex};
            std::vector<int> descending_trajectory = {vertex};
            int asc_vertex = vertex;
            int desc_vertex = vertex;

            // Compute ascending trajectory if not at local maximum
            if (!is_max) {
                while (true) {
                    auto valid_neighbors = get_valid_neighbors(asc_vertex, true);

                    if (valid_neighbors.empty()) {
                        // This should never happen for non-extremal vertices
                        Rprintf("\nCRITICAL ERROR: No valid ascending neighbors found for non-maximum vertex!\n");
                        Rprintf("Details:\n");
                        Rprintf("  - Current vertex: %d\n", asc_vertex);
                        Rprintf("  - Function value: %.6f\n", Ey[asc_vertex]);
                        Rprintf("  - Trajectory so far: ");
                        for (int v : ascending_trajectory) {
                            Rprintf("%d ", v);
                        }

                        Rprintf("\n\nNeighborhood Analysis:\n");
                        Rprintf("  Neighbors in graph:");
                        for (int n : graph[asc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rprintf("\n  Neighbors in core_graph:");
                        for (int n : core_graph[asc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rf_error("\nFATAL: Found a non-maximum vertex with no valid ascending neighbors.\n"
                                 "This violates the properties of a valid Morse function on the graph.\n"
                                 "Please verify that:\n"
                                 "1. The function values form a valid Morse function\n"
                                 "2. The graph connectivity is correct\n"
                                 "3. All paths in core_graph preserve monotonicity");
                    }

                    auto path = valid_neighbors[0].path;
                    ascending_trajectory.insert(ascending_trajectory.end(),
                                                path.begin() + 1, path.end());
                    asc_vertex = valid_neighbors[0].vertex;
                }

                // Check if endpoint is local maximum
                if (is_local_maximum(asc_vertex, graph, Ey)) {
                    local_maxima_set.insert(asc_vertex);
                }
            }

            // Compute descending trajectory if not at local minimum
            if (!is_min) {
                while (true) {
                    auto valid_neighbors = get_valid_neighbors(desc_vertex, false);

                    if (valid_neighbors.empty()) {
                        // This should never happen for non-extremal vertices
                        Rprintf("\nCRITICAL ERROR: No valid descending neighbors found for non-minimum vertex!\n");
                        Rprintf("Details:\n");
                        Rprintf("  - Current vertex: %d\n", desc_vertex);
                        Rprintf("  - Function value: %.6f\n", Ey[desc_vertex]);
                        Rprintf("  - Trajectory so far: ");
                        for (int v : descending_trajectory) {
                            Rprintf("%d ", v);
                        }

                        Rprintf("\n\nNeighborhood Analysis:\n");
                        Rprintf("  Neighbors in graph:");
                        for (int n : graph[desc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rprintf("\n  Neighbors in core_graph:");
                        for (int n : core_graph[desc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rf_error("\nFATAL: Found a non-minimum vertex with no valid descending neighbors.\n"
                                 "This violates the properties of a valid Morse function on the graph.\n"
                                 "Please verify that:\n"
                                 "1. The function values form a valid Morse function\n"
                                 "2. The graph connectivity is correct\n"
                                 "3. All paths in core_graph preserve monotonicity");
                    }

                    auto path = valid_neighbors[0].path;
                    descending_trajectory.insert(descending_trajectory.end(),
                                                 path.begin() + 1, path.end());
                    desc_vertex = valid_neighbors[0].vertex;
                }

                // Check if endpoint is local minimum
                if (is_local_minimum(desc_vertex, graph, Ey)) {
                    local_minima_set.insert(desc_vertex);
                }
            }

            if (is_max) {
                // Single vertex trajectory for local maximum
                SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, 1));
                INTEGER(trajectory)[0] = vertex;
                SET_VECTOR_ELT(trajectories, vertex, trajectory);
                UNPROTECT(1);
            } else {
                // Construct full trajectory
                descending_trajectory.erase(descending_trajectory.begin());
                int tr_size = ascending_trajectory.size() + descending_trajectory.size();
                SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, tr_size));
                int* trajectory_ptr = INTEGER(trajectory);

                std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
                std::copy(ascending_trajectory.begin(), ascending_trajectory.end(),
                          trajectory_ptr + descending_trajectory.size());

                SET_VECTOR_ELT(trajectories, vertex, trajectory);
                UNPROTECT(1);
            }

            // Update connectivity maps and pro-cells if trajectory connects local extrema
            if (is_local_minimum(desc_vertex, graph, Ey) &&
                is_local_maximum(asc_vertex, graph, Ey)) {

                // Update connectivity maps
                lmax_to_lmin_map[asc_vertex].insert(desc_vertex);
                lmin_to_lmax_map[desc_vertex].insert(asc_vertex);

                // Update pro-cell
                std::pair<int,int> cell_key(asc_vertex, desc_vertex);

                // Add all vertices from the trajectory to the pro-cell
                ms_procells[cell_key].insert(descending_trajectory.begin(),
                                             descending_trajectory.end());
                ms_procells[cell_key].insert(ascending_trajectory.begin(),
                                             ascending_trajectory.end());
            }
        }

        SET_VECTOR_ELT(result, 0, trajectories);
        UNPROTECT(1); // trajectories
    }

    // 1-4: lmax_to_lmin, lmin_to_lmax, local_maxima, local_minima
    {
        // Convert maps to R lists
        SEXP lmax_to_lmin       = PROTECT(Rf_allocVector(VECSXP, local_maxima_set.size()));
        SEXP lmax_to_lmin_names = PROTECT(Rf_allocVector(STRSXP, local_maxima_set.size()));

        SEXP lmin_to_lmax       = PROTECT(Rf_allocVector(VECSXP, local_minima_set.size()));
        SEXP lmin_to_lmax_names = PROTECT(Rf_allocVector(STRSXP, local_minima_set.size()));

        // Convert local maxima/minima sets to R vectors
        SEXP local_maxima = PROTECT(Rf_allocVector(INTSXP, local_maxima_set.size()));
        SEXP local_minima = PROTECT(Rf_allocVector(INTSXP, local_minima_set.size()));

        // Fill local maxima vector and its connectivity list
        int max_idx = 0;
        for (int lmax : local_maxima_set) {
            INTEGER(local_maxima)[max_idx] = lmax;

            const auto& connected_mins = lmax_to_lmin_map[lmax];
            SEXP mins = PROTECT(Rf_allocVector(INTSXP, connected_mins.size()));
            std::copy(connected_mins.begin(), connected_mins.end(), INTEGER(mins));
            SET_VECTOR_ELT(lmax_to_lmin, max_idx, mins);
            std::string lmax_str = std::to_string(lmax + 1);  // Turning lmax into const char* after turning it into 1-base interger
            const char* lmax_char = lmax_str.c_str();
            SET_STRING_ELT(lmax_to_lmin_names, max_idx, Rf_mkChar(lmax_char));
            UNPROTECT(1);

            max_idx++;
        }
        Rf_setAttrib(lmax_to_lmin, R_NamesSymbol, lmax_to_lmin_names);
        UNPROTECT(1); // lmax_to_lmin_names

        // Fill local minima vector and its connectivity list
        int min_idx = 0;
        for (int lmin : local_minima_set) {
            INTEGER(local_minima)[min_idx] = lmin;

            const auto& connected_maxs = lmin_to_lmax_map[lmin];
            SEXP maxs = PROTECT(Rf_allocVector(INTSXP, connected_maxs.size()));
            std::copy(connected_maxs.begin(), connected_maxs.end(), INTEGER(maxs));
            SET_VECTOR_ELT(lmin_to_lmax, min_idx, maxs);
            std::string lmin_str = std::to_string(lmin + 1);  // Turning lmin into const char*
            const char* lmin_char = lmin_str.c_str();
            SET_STRING_ELT(lmin_to_lmax_names, min_idx, Rf_mkChar(lmin_char));
            UNPROTECT(1);

            min_idx++;
        }
        Rf_setAttrib(lmin_to_lmax, R_NamesSymbol, lmin_to_lmax_names);
        UNPROTECT(1); // lmin_to_lmax_names

        SET_VECTOR_ELT(result, 1, lmax_to_lmin);
        SET_VECTOR_ELT(result, 2, lmin_to_lmax);
        SET_VECTOR_ELT(result, 3, local_maxima);
        SET_VECTOR_ELT(result, 4, local_minima);
        UNPROTECT(4); // lmax_to_lmin, lmin_to_lmax, local_maxima, local_minima
    }

    // Compute MS cells from pro-cells
    for (const auto& [key, procell] : ms_procells) {

        if (procell.size() > 2) {
            // Remove extrema from the set for component calculation
            std::set<int> cell_vertices = procell;

            // print_set(cell_vertices, "cell_vertices");

            cell_vertices.erase(key.first);   // Remove local maximum
            cell_vertices.erase(key.second);  // Remove local minimum

            // print_set(cell_vertices, "cell_vertices");

            // Find connected components
            auto components = count_subgraph_set_components(core_graph, cell_vertices);

            // print_map_to_set(components, "components");

            // Store components as MS cells
            std::vector<std::set<int>> cell_components;
            for (const auto& [comp_id, comp_vertices] : components) {
                std::set<int> new_comp = comp_vertices;  // Make a copy
                new_comp.insert(key.first);
                new_comp.insert(key.second);
                cell_components.push_back(new_comp);
            }
            ms_cells[key] = cell_components;
        } else {
            std::vector<std::set<int>> cell_components;
            cell_components.push_back(procell);
            ms_cells[key] = cell_components;
        }
    }

    // 5-7: procell_keys, procells, cells
    {
        // Convert pro-cells and cells to R structures
        SEXP procells = PROTECT(Rf_allocVector(VECSXP, ms_procells.size()));
        SEXP cells = PROTECT(Rf_allocVector(VECSXP, ms_cells.size()));
        SEXP procell_keys = PROTECT(Rf_allocVector(VECSXP, ms_procells.size()));

        int pcell_idx = 0;
        for (const auto& [key, vertices] : ms_procells) {
            // Create key pair
            SEXP key_pair = PROTECT(Rf_allocVector(INTSXP, 2));
            INTEGER(key_pair)[0] = key.first;
            INTEGER(key_pair)[1] = key.second;
            SET_VECTOR_ELT(procell_keys, pcell_idx, key_pair);

            // Create vertex set
            SEXP vertex_set = PROTECT(Rf_allocVector(INTSXP, vertices.size()));
            std::copy(vertices.begin(), vertices.end(), INTEGER(vertex_set));
            SET_VECTOR_ELT(procells, pcell_idx, vertex_set);

            // Create cells list for this pro-cell
            const auto& cell_components = ms_cells[key];
            SEXP components_list = PROTECT(Rf_allocVector(VECSXP, cell_components.size()));

            for (size_t i = 0; i < cell_components.size(); ++i) {
                SEXP component = PROTECT(Rf_allocVector(INTSXP, cell_components[i].size()));
                std::copy(cell_components[i].begin(), cell_components[i].end(),
                          INTEGER(component));
                SET_VECTOR_ELT(components_list, i, component);
                UNPROTECT(1);
            }

            SET_VECTOR_ELT(cells, pcell_idx, components_list);

            UNPROTECT(3);  // key_pair, vertex_set, components_list
            pcell_idx++;
        }

        SET_VECTOR_ELT(result, 5, procell_keys);
        SET_VECTOR_ELT(result, 6, procells);
        SET_VECTOR_ELT(result, 7, cells);
        UNPROTECT(3); // procell_keys, procells, cells
    }

    UNPROTECT(1);
    return result;
}

/**
* @brief Computes the Morse-Smale complex using shortest hop distances in the graph
*
* @details This function computes the Morse-Smale complex for a function defined on the vertices
* of a graph by finding ascending and descending trajectories that follow paths of steepest
* ascent/descent. The key features of this implementation include:
*
* - Uses pre-computed hop distances from s_hop_list to find nearest valid neighbors
* - Prioritizes neighbors at minimum hop distance with maximum function value difference
* - Validates paths between neighbors using core_graph connectivity
* - Constructs gradient flow trajectories, pro-cells and cells of the MS complex
*
* The algorithm works as follows:
* 1. For each vertex v:
*    - If not a local maximum, compute ascending trajectory by iteratively:
*      * Finding valid neighbors at minimum hop distance with higher function values
*      * Selecting neighbor with maximum value difference
*      * Following shortest path in core_graph to that neighbor
*      * Continuing until a local maximum is reached
*    - Similarly compute descending trajectory to local minimum if needed
*    - Store complete trajectory and update connectivity maps
*
* 2. Construct MS complex components:
*    - Pro-cells: Collection of vertices in trajectories between (max, min) pairs
*    - Cells: Connected components of pro-cells after removing extrema
*    - Update connectivity maps between maxima and minima
*
* @param s_graph SEXP containing h-hop neighbor graph adjacency lists. Each list[i] contains
*               indices of vertices reachable from vertex i within h hops
* @param s_hop_list SEXP containing hop distances corresponding to s_graph adjacency lists.
*                   hop_list[i][j] gives number of hops from vertex i to its j-th neighbor
* @param s_core_graph SEXP containing core graph adjacency lists used for finding valid paths
*                     between neighbor vertices
* @param s_Ey SEXP containing function values at vertices (double array)
*
* @return SEXP List containing:
*   - trajectories: List of integer vectors with gradient trajectories for each vertex
*   - lmax_to_lmin: Named list mapping local maxima to connected local minima
*   - lmin_to_lmax: Named list mapping local minima to connected local maxima
*   - local_maxima: Integer vector of local maxima indices
*   - local_minima: Integer vector of local minima indices
*   - procell_keys: List of (max_idx, min_idx) pairs defining pro-cells
*   - procells: List of integer vectors containing vertices in each pro-cell
*   - cells: List of lists containing connected components for each pro-cell
*
* @throws R Rf_error if:
*   - Input lengths don't match
*   - No valid path exists between neighbor vertices in core_graph
*   - Non-extremal vertex has no valid ascending/descending neighbors
*
* @note The function assumes:
*   - s_graph represents an h-hop neighbor graph
*   - s_hop_list contains valid hop distances matching s_graph
*   - s_core_graph provides valid paths between neighbor vertices
*   - Function values in s_Ey form a valid Morse function on the graph
*
* @see shortest_path_by_hops() for path finding in core graph
* @see is_local_maximum(), is_local_minimum() for extrema detection
* @see count_subgraph_set_components() for cell decomposition
*/
SEXP S_graph_MS_cx_using_short_h_hops(SEXP s_graph,
                                      SEXP s_hop_list,
                                      SEXP s_core_graph,
                                      SEXP s_Ey) {
    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_hop_list) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        Rf_error("All input lengths (s_graph, s_hop_list, s_core_graph, s_Ey) must be %d", n_vertices);
    }

    // Convert input data
    std::vector<std::vector<int>> graph      = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<int>> core_graph = convert_adj_list_from_R(s_core_graph);
    std::vector<std::vector<int>> hop_list   = convert_adj_list_from_R(s_hop_list);
    double* Ey = REAL(s_Ey);

    // Initialize maps for extrema connectivity
    std::map<int, std::set<int>> lmax_to_lmin_map;
    std::map<int, std::set<int>> lmin_to_lmax_map;
    std::set<int> local_maxima_set;
    std::set<int> local_minima_set;

    // Update result list to include pro-cells and cells
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 8));

    // Update names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 8));
        SET_STRING_ELT(names, 0, Rf_mkChar("trajectories"));
        SET_STRING_ELT(names, 1, Rf_mkChar("lmax_to_lmin"));
        SET_STRING_ELT(names, 2, Rf_mkChar("lmin_to_lmax"));
        SET_STRING_ELT(names, 3, Rf_mkChar("local_maxima"));
        SET_STRING_ELT(names, 4, Rf_mkChar("local_minima"));
        SET_STRING_ELT(names, 5, Rf_mkChar("procell_keys"));
        SET_STRING_ELT(names, 6, Rf_mkChar("procells"));
        SET_STRING_ELT(names, 7, Rf_mkChar("cells"));
        Rf_setAttrib(result, R_NamesSymbol, names);
        UNPROTECT(1); // names
    }

    // Helper function to get valid neighbors using hop distances
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<hop_neighbor_info_t> neighbors;
        int min_valid_hop = -1;
        std::vector<int> path;

        // First pass: find minimum hop distance with valid neighbors
        for (size_t i = 0; i < graph[current_vertex].size(); i++) {
            int neighbor = graph[current_vertex][i];
            int hop_dist = hop_list[current_vertex][i];

            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0) {
                if (min_valid_hop == -1 || hop_dist < min_valid_hop) {
                    min_valid_hop = hop_dist;
                }
            }
        }

        // Second pass: collect all neighbors at minimum hop distance
        if (min_valid_hop != -1) {
            for (size_t i = 0; i < graph[current_vertex].size(); i++) {
                int neighbor = graph[current_vertex][i];
                int hop_dist = hop_list[current_vertex][i];

                if (hop_dist == min_valid_hop) {
                    double diff = ascending ?
                        Ey[neighbor] - Ey[current_vertex] :
                        Ey[current_vertex] - Ey[neighbor];

                    if (diff > 0) {
                        path = shortest_path_by_hops(core_graph, current_vertex, neighbor);
                        if (path.empty()) {
                            Rprintf("\nCRITICAL ERROR: No path exists in core_graph between vertices %d and %d\n",
                                    current_vertex, neighbor);
                            Rf_error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                        }
                        neighbors.emplace_back(hop_neighbor_info_t{neighbor, hop_dist, diff, path});
                    }
                }
            }

            // Sort by difference value in descending order
            std::sort(neighbors.begin(), neighbors.end(),
                     [](const hop_neighbor_info_t& a, const hop_neighbor_info_t& b) {
                         return a.diff_value > b.diff_value;
                     });
        }

        return neighbors;
    };

    // Add map for pro-cells
    std::map<std::pair<int,int>, std::set<int>> ms_procells;
    std::map<std::pair<int,int>, std::vector<std::set<int>>> ms_cells;

    // 0:
    {
        SEXP trajectories = PROTECT(Rf_allocVector(VECSXP, n_vertices));
        for (int i = 0; i < n_vertices; i++) {
            SET_VECTOR_ELT(trajectories, i, R_NilValue);
        }

        for (int vertex = 0; vertex < n_vertices; vertex++) {
            bool is_max = is_local_maximum(vertex, graph, Ey);
            bool is_min = is_local_minimum(vertex, graph, Ey);

            if (is_max) {
                local_maxima_set.insert(vertex);
            }
            if (is_min) {
                local_minima_set.insert(vertex);
            }

            // Compute ascending and descending trajectories
            std::vector<int> ascending_trajectory = {vertex};
            std::vector<int> descending_trajectory = {vertex};
            int asc_vertex = vertex;
            int desc_vertex = vertex;

            // Compute ascending trajectory if not at local maximum
            if (!is_max) {
                while (true) {
                    auto valid_neighbors = get_valid_neighbors(asc_vertex, true);

                    if (valid_neighbors.empty()) {
                        // This should never happen for non-extremal vertices
                        Rprintf("\nCRITICAL ERROR: No valid ascending neighbors found for non-maximum vertex!\n");
                        Rprintf("Details:\n");
                        Rprintf("  - Current vertex: %d\n", asc_vertex);
                        Rprintf("  - Function value: %.6f\n", Ey[asc_vertex]);
                        Rprintf("  - Trajectory so far: ");
                        for (int v : ascending_trajectory) {
                            Rprintf("%d ", v);
                        }

                        //Rprintf("\n\nNeighborhood Analysis:\n");
                        Rprintf("  Neighbors in graph:");
                        for (int n : graph[asc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rf_error("\nFATAL: Found a non-maximum vertex with no valid ascending neighbors.\n"
                                 "This violates the properties of a valid Morse function on the graph.\n"
                                 "Please verify that:\n"
                                 "1. The function values form a valid Morse function\n"
                                 "2. The graph connectivity is correct");
                    }

                    auto path = valid_neighbors[0].path;
                    ascending_trajectory.insert(ascending_trajectory.end(),
                                                path.begin() + 1, path.end());
                    asc_vertex = valid_neighbors[0].vertex;
                    if (is_local_maximum(asc_vertex, graph, Ey)) {
                        break;
                    }
                }
            }

            // Compute descending trajectory if not at local minimum
            if (!is_min) {
                while (true) {
                    auto valid_neighbors = get_valid_neighbors(desc_vertex, false);

                    if (valid_neighbors.empty()) {
                        // This should never happen for non-extremal vertices
                        Rprintf("\nCRITICAL ERROR: No valid descending neighbors found for non-minimum vertex!\n");
                        Rprintf("Details:\n");
                        Rprintf("  - Current vertex: %d\n", desc_vertex);
                        Rprintf("  - Function value: %.6f\n", Ey[desc_vertex]);
                        Rprintf("  - Trajectory so far: ");
                        for (int v : descending_trajectory) {
                            Rprintf("%d ", v);
                        }

                        //Rprintf("\n\nNeighborhood Analysis:\n");
                        Rprintf("  Neighbors in graph:");
                        for (int n : graph[desc_vertex]) {
                            Rprintf(" %d(%.6f)", n, Ey[n]);
                        }

                        Rf_error("\nFATAL: Found a non-minimum vertex with no valid descending neighbors.\n"
                                 "This violates the properties of a valid Morse function on the graph.\n"
                                 "Please verify that:\n"
                                 "1. The function values form a valid Morse function\n"
                                 "2. The graph connectivity is correct");
                    }

                    auto path = valid_neighbors[0].path;
                    descending_trajectory.insert(descending_trajectory.end(),
                                                 path.begin() + 1, path.end());
                    desc_vertex = valid_neighbors[0].vertex;
                    if (is_local_minimum(desc_vertex, graph, Ey)) {
                        break;
                    }
                }
            }

            if (is_max) {
                // Single vertex trajectory for local maximum
                SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, 1));
                INTEGER(trajectory)[0] = vertex;
                SET_VECTOR_ELT(trajectories, vertex, trajectory);
                UNPROTECT(1);
            } else {
                // Construct full trajectory
                descending_trajectory.erase(descending_trajectory.begin());
                int tr_size = ascending_trajectory.size() + descending_trajectory.size();
                SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, tr_size));
                int* trajectory_ptr = INTEGER(trajectory);

                std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
                std::copy(ascending_trajectory.begin(), ascending_trajectory.end(),
                          trajectory_ptr + descending_trajectory.size());

                SET_VECTOR_ELT(trajectories, vertex, trajectory);
                UNPROTECT(1);
            }

            // Update connectivity maps and pro-cells if trajectory connects local extrema
            if (is_local_minimum(desc_vertex, graph, Ey) &&
                is_local_maximum(asc_vertex, graph, Ey)) {

                // Update connectivity maps
                lmax_to_lmin_map[asc_vertex].insert(desc_vertex);
                lmin_to_lmax_map[desc_vertex].insert(asc_vertex);

                // Update pro-cell
                std::pair<int,int> cell_key(asc_vertex, desc_vertex);

                // Add all vertices from the trajectory to the pro-cell
                ms_procells[cell_key].insert(descending_trajectory.begin(),
                                             descending_trajectory.end());
                ms_procells[cell_key].insert(ascending_trajectory.begin(),
                                             ascending_trajectory.end());
            }
        }

        SET_VECTOR_ELT(result, 0, trajectories);
        UNPROTECT(1); // trajectories
    }

    // 1-4: lmax_to_lmin, lmin_to_lmax, local_maxima, local_minima
    {
        // Convert maps to R lists
        SEXP lmax_to_lmin       = PROTECT(Rf_allocVector(VECSXP, local_maxima_set.size()));
        SEXP lmax_to_lmin_names = PROTECT(Rf_allocVector(STRSXP, local_maxima_set.size()));

        SEXP lmin_to_lmax       = PROTECT(Rf_allocVector(VECSXP, local_minima_set.size()));
        SEXP lmin_to_lmax_names = PROTECT(Rf_allocVector(STRSXP, local_minima_set.size()));

        // Convert local maxima/minima sets to R vectors
        SEXP local_maxima = PROTECT(Rf_allocVector(INTSXP, local_maxima_set.size()));
        SEXP local_minima = PROTECT(Rf_allocVector(INTSXP, local_minima_set.size()));

        // Fill local maxima vector and its connectivity list
        int max_idx = 0;
        for (int lmax : local_maxima_set) {
            INTEGER(local_maxima)[max_idx] = lmax;

            const auto& connected_mins = lmax_to_lmin_map[lmax];
            SEXP mins = PROTECT(Rf_allocVector(INTSXP, connected_mins.size()));
            std::copy(connected_mins.begin(), connected_mins.end(), INTEGER(mins));
            SET_VECTOR_ELT(lmax_to_lmin, max_idx, mins);
            std::string lmax_str = std::to_string(lmax + 1);  // Turning lmax into const char* after turning it into 1-base interger
            const char* lmax_char = lmax_str.c_str();
            SET_STRING_ELT(lmax_to_lmin_names, max_idx, Rf_mkChar(lmax_char));
            UNPROTECT(1);

            max_idx++;
        }
        Rf_setAttrib(lmax_to_lmin, R_NamesSymbol, lmax_to_lmin_names);
        UNPROTECT(1); // lmax_to_lmin_names

        // Fill local minima vector and its connectivity list
        int min_idx = 0;
        for (int lmin : local_minima_set) {
            INTEGER(local_minima)[min_idx] = lmin;

            const auto& connected_maxs = lmin_to_lmax_map[lmin];
            SEXP maxs = PROTECT(Rf_allocVector(INTSXP, connected_maxs.size()));
            std::copy(connected_maxs.begin(), connected_maxs.end(), INTEGER(maxs));
            SET_VECTOR_ELT(lmin_to_lmax, min_idx, maxs);
            std::string lmin_str = std::to_string(lmin + 1);  // Turning lmin into const char*
            const char* lmin_char = lmin_str.c_str();
            SET_STRING_ELT(lmin_to_lmax_names, min_idx, Rf_mkChar(lmin_char));
            UNPROTECT(1);

            min_idx++;
        }
        Rf_setAttrib(lmin_to_lmax, R_NamesSymbol, lmin_to_lmax_names);
        UNPROTECT(1); // lmin_to_lmax_names

        SET_VECTOR_ELT(result, 1, lmax_to_lmin);
        SET_VECTOR_ELT(result, 2, lmin_to_lmax);
        SET_VECTOR_ELT(result, 3, local_maxima);
        SET_VECTOR_ELT(result, 4, local_minima);
        UNPROTECT(4);
    }

    // Compute MS cells from pro-cells
    for (const auto& [key, procell] : ms_procells) {

        if (procell.size() > 2) {
            // Remove extrema from the set for component calculation
            std::set<int> cell_vertices = procell;

            // print_set(cell_vertices, "cell_vertices");

            cell_vertices.erase(key.first);   // Remove local maximum
            cell_vertices.erase(key.second);  // Remove local minimum

            // print_set(cell_vertices, "cell_vertices");

            // Find connected components
            auto components = count_subgraph_set_components(core_graph, cell_vertices);
            if (components.empty() && !cell_vertices.empty()) {
                Rf_warning("Empty components found for non-empty cell vertices");
            }
            // print_map_to_set(components, "components");

            // Store components as MS cells
            std::vector<std::set<int>> cell_components;
            for (const auto& [comp_id, comp_vertices] : components) {
                std::set<int> new_comp = comp_vertices;  // Make a copy
                new_comp.insert(key.first);
                new_comp.insert(key.second);
                cell_components.push_back(new_comp);
            }
            ms_cells[key] = cell_components;
        } else {
            std::vector<std::set<int>> cell_components;
            cell_components.push_back(procell);
            ms_cells[key] = cell_components;
        }
    }

    // 5-7: procell_keys, procells, cells
    {
        // Convert pro-cells and cells to R structures
        SEXP procells = PROTECT(Rf_allocVector(VECSXP, ms_procells.size()));
        SEXP cells = PROTECT(Rf_allocVector(VECSXP, ms_cells.size()));
        SEXP procell_keys = PROTECT(Rf_allocVector(VECSXP, ms_procells.size()));

        int pcell_idx = 0;
        for (const auto& [key, vertices] : ms_procells) {
            // Create key pair
            SEXP key_pair = PROTECT(Rf_allocVector(INTSXP, 2));
            INTEGER(key_pair)[0] = key.first;
            INTEGER(key_pair)[1] = key.second;
            SET_VECTOR_ELT(procell_keys, pcell_idx, key_pair);

            // Create vertex set
            SEXP vertex_set = PROTECT(Rf_allocVector(INTSXP, vertices.size()));
            std::copy(vertices.begin(), vertices.end(), INTEGER(vertex_set));
            SET_VECTOR_ELT(procells, pcell_idx, vertex_set);

            // Create cells list for this pro-cell
            const auto& cell_components = ms_cells[key];
            SEXP components_list = PROTECT(Rf_allocVector(VECSXP, cell_components.size()));

            for (size_t i = 0; i < cell_components.size(); ++i) {
                SEXP component = PROTECT(Rf_allocVector(INTSXP, cell_components[i].size()));
                std::copy(cell_components[i].begin(), cell_components[i].end(),
                          INTEGER(component));
                SET_VECTOR_ELT(components_list, i, component);
                UNPROTECT(1);
            }

            SET_VECTOR_ELT(cells, pcell_idx, components_list);

            UNPROTECT(3);  // key_pair, vertex_set, components_list
            pcell_idx++;
        }

        SET_VECTOR_ELT(result, 5, procell_keys);
        SET_VECTOR_ELT(result, 6, procells);
        SET_VECTOR_ELT(result, 7, cells);
        UNPROTECT(3); // procell_keys, procells, cells
    }

    UNPROTECT(1);
    return result;
}



/**
 * @brief Computes constrained gradient flow trajectories on a graph while preventing basin-jumping
 *
 * @param s_graph SEXP containing h-th power graph pIG_k^h(X) as adjacency lists
 * @param s_core_graph SEXP containing core graph pIG_k(X) as adjacency lists
 * @param s_Ey SEXP containing function values at vertices
 *
 * @return SEXP List of integer vectors, each containing complete trajectory
 *         (descending path reversed + ascending path) for corresponding vertex
 *
 * @throws R Rf_error if input lengths don't match
 */
SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {
    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        Rf_error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph      = convert_adj_list_from_R(s_graph);
    std::vector<std::vector<int>> core_graph = convert_adj_list_from_R(s_core_graph);
    double* Ey = REAL(s_Ey);

    SEXP trajectories = PROTECT(Rf_allocVector(VECSXP, n_vertices));

    // Helper function to get valid neighbors with their paths
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : graph[current_vertex]) {
            std::vector<int> path = shortest_path_by_hops(core_graph, current_vertex, neighbor);
            if (path.empty()) continue;

            double diff = ascending ?
                         Ey[neighbor] - Ey[current_vertex] :
                         Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                is_path_monotonic_increasing(path, Ey) :
                is_path_monotonic_decreasing(path, Ey))) {
                neighbors.emplace_back(neighbor_info_t{neighbor, path, static_cast<int>(path.size() - 1), diff});
            }
        }

        // Sort by difference value in descending order
        std::sort(neighbors.begin(), neighbors.end(),
                 [](const neighbor_info_t& a, const neighbor_info_t& b) {
                     return a.diff_value > b.diff_value;
                 });

        return neighbors;
    };

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        bool is_max = is_local_maximum(vertex, graph, Ey);

        // Ascending trajectory
        std::vector<int> ascending_trajectory = {vertex};
        int current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, true);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                ascending_trajectory.insert(ascending_trajectory.end(),
                                         path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }
        }

        // Descending trajectory
        std::vector<int> descending_trajectory = {vertex};
        current_vertex = vertex;

        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(current_vertex, false);
                if (valid_neighbors.empty()) break;
                auto path = valid_neighbors[0].path;
                descending_trajectory.insert(descending_trajectory.end(),
                                          path.begin() + 1, path.end());
                current_vertex = valid_neighbors[0].vertex;
            }
        }

        if (is_max) {
            // Single vertex trajectory for local maximum
            SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
            // Construct full trajectory
            descending_trajectory.erase(descending_trajectory.begin());
            int tr_size = ascending_trajectory.size() + descending_trajectory.size();
            SEXP trajectory = PROTECT(Rf_allocVector(INTSXP, tr_size));
            int* trajectory_ptr = INTEGER(trajectory);

            std::copy(descending_trajectory.rbegin(), descending_trajectory.rend(), trajectory_ptr);
            std::copy(ascending_trajectory.begin(), ascending_trajectory.end(),
                     trajectory_ptr + descending_trajectory.size());

            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        }
    }

    UNPROTECT(1);
    return trajectories;
}


// First approach: Trajectory clustering using efficient similarity computation
struct TrajectoryCluster {
    std::vector<size_t> trajectory_indices;
    std::vector<int> representative_path;
    double average_similarity;
};

class MSCellComputer {
private:
    // Efficient data structure for trajectory overlap computation
    struct TrajectoryIndex {
        std::unordered_map<int, std::set<size_t>> vertex_to_trajectories;

        void build(const std::vector<std::vector<int>>& trajectories) {
            for (size_t i = 0; i < trajectories.size(); i++) {
                for (int vertex : trajectories[i]) {
                    vertex_to_trajectories[vertex].insert(i);
                }
            }
        }

        // Get trajectories that share vertices with given trajectory
        std::set<size_t> get_potential_matches(const std::vector<int>& trajectory) {
            std::set<size_t> candidates;
            for (int vertex : trajectory) {
                const auto& matching_trajectories = vertex_to_trajectories[vertex];
                candidates.insert(matching_trajectories.begin(), matching_trajectories.end());
            }
            return candidates;
        }
    };

    // Compute Jaccard similarity efficiently
    double compute_similarity(const std::vector<int>& traj1,
                            const std::vector<int>& traj2) {
        std::set<int> union_set(traj1.begin(), traj1.end());
        union_set.insert(traj2.begin(), traj2.end());

        std::set<int> intersection;
        std::set_intersection(traj1.begin(), traj1.end(),
                            traj2.begin(), traj2.end(),
                            std::inserter(intersection, intersection.begin()));

        return static_cast<double>(intersection.size()) / union_set.size();
    }

public:
    // Main function to compute cells using trajectory clustering
    std::vector<std::set<int>> compute_cells(
        const MS_complex_t& ms_cx,
        const std::pair<int, int>& extrema_pair,
        double similarity_threshold = 0.9) {

        // Get relevant trajectories for this cell
        std::vector<std::vector<int>> cell_trajectories;
        std::vector<size_t> trajectory_indices;

        for (const auto& [cell, traj_indices] : ms_cx.cell_trajectories) {
            if (cell.lmin == extrema_pair.first &&
                cell.lmax == extrema_pair.second) {
                for (size_t idx : traj_indices) {
                    cell_trajectories.push_back(ms_cx.unique_trajectories[idx]);
                    trajectory_indices.push_back(idx);
                }
            }
        }

        // Build index for efficient similarity computation
        TrajectoryIndex traj_index;
        traj_index.build(cell_trajectories);

        // Cluster trajectories
        std::vector<TrajectoryCluster> clusters;
        std::vector<bool> processed(cell_trajectories.size(), false);

        for (size_t i = 0; i < cell_trajectories.size(); i++) {
            if (processed[i]) continue;

            TrajectoryCluster new_cluster;
            new_cluster.trajectory_indices.push_back(trajectory_indices[i]);
            new_cluster.representative_path = cell_trajectories[i];

            // Find similar trajectories efficiently
            auto candidates = traj_index.get_potential_matches(cell_trajectories[i]);

            for (size_t j : candidates) {
                if (j <= i || processed[j]) continue;

                double similarity = compute_similarity(
                    cell_trajectories[i], cell_trajectories[j]);

                if (similarity >= similarity_threshold) {
                    new_cluster.trajectory_indices.push_back(trajectory_indices[j]);
                    processed[j] = true;

                    // Update representative path (can be optimized further)
                    std::set<int> merged_path(new_cluster.representative_path.begin(),
                                           new_cluster.representative_path.end());
                    merged_path.insert(cell_trajectories[j].begin(),
                                     cell_trajectories[j].end());
                    new_cluster.representative_path = std::vector<int>(
                        merged_path.begin(), merged_path.end());
                }
            }

            clusters.push_back(new_cluster);
            processed[i] = true;
        }

        // Convert clusters to cells
        std::vector<std::set<int>> cells;
        for (const auto& cluster : clusters) {
            std::set<int> component(cluster.representative_path.begin(),
                                  cluster.representative_path.end());
            cells.push_back(component);
        }

        return cells;
    }
};

/**
 * @brief Computes the Morse-Smale complex for a graph with associated scalar values
 *
 * This function analyzes a graph with scalar values assigned to its vertices to identify
 * its Morse-Smale complex structure. It identifies local maxima and minima, computes
 * gradient trajectories from each vertex to its associated extrema, and determines the
 * connectivity between critical points. It also computes proto-cells, their decomposition
 * into Morse-Smale cells, and associates trajectories with each cell.
 *
 * @param adj_list The adjacency list representation of the graph. adj_list[i] contains
 *                 the vertices adjacent to vertex i
 * @param hop_list A list containing hop distances between vertices. hop_list[i][j] contains
 *                the hop distance between vertex i and its j-th neighbor in adj_list[i]
 * @param core_adj_list A possibly reduced adjacency list used for component calculations
 *                     in cell decomposition
 * @param shortest_paths A map storing pre-computed shortest paths between vertex pairs.
 *                      The key is a pair of vertices (source, target), and the value
 *                      is the vector of vertices representing the shortest path
 * @param Ey A vector of scalar values associated with each vertex. Ey[i] is the scalar
 *          value at vertex i
 *
 * @return A MS_complex_t structure containing:
 *         - Local maxima and minima sets
 *         - Connectivity maps between extrema
 *         - Proto-cells mapping extrema pairs to vertex sets
 *         - Decomposed Morse-Smale cells
 *         - Unique gradient trajectories
 *         - Mapping of cells to their associated trajectories
 *
 * @note The function uses a two-pass neighbor selection strategy:
 *       1. First pass identifies the minimum hop distance with valid neighbors
 *       2. Second pass collects all neighbors at that minimum hop distance
 *       This ensures consistent neighbor selection based on both graph distance
 *       and function value differences.
 *
 * @details The algorithm proceeds in several steps:
 * 1. For each vertex:
 *    - Determines if it is a local maximum or minimum
 *    - Computes a complete gradient trajectory by:
 *      * Following steepest ascending path to a local maximum
 *      * Following steepest descending path to a local minimum
 *      * For each step:
 *        - Finds valid neighbors at minimum hop distance
 *        - Selects neighbor with maximum function value difference
 *        - Uses pre-computed shortest path to that neighbor
 *    - Updates the connectivity maps between extrema
 *    - Constructs proto-cells from complete trajectories
 *
 * 2. For each proto-cell:
 *    - Removes the extrema vertices
 *    - Identifies connected components in the remaining vertices
 *    - Creates separate cells for each component
 *    - Associates relevant trajectories with each cell
 *
 * The neighbor selection process:
 * - Does not impose monotonicity conditions on paths
 * - Considers both hop distance and function value differences
 * - Selects neighbors at minimum hop distance with maximum function value difference
 * - Uses pre-computed shortest paths between selected neighbors
 *
 * @Rf_warning The function assumes that:
 * - All vertex indices are valid (0 to adj_list.size()-1)
 * - The graph is undirected
 * - Ey contains valid scalar values for all vertices
 * - shortest_paths contains valid paths when needed
 * - hop_list contains valid hop distances for all adjacent vertices
 *
 * Example usage:
 * @code
 * std::vector<std::vector<int>> adj_list = {...};  // Graph adjacency list
 * std::vector<std::vector<int>> hop_list = {...};  // Hop distances
 * std::vector<std::vector<int>> core_adj_list = {...};  // Core graph adjacency list
 * std::map<std::pair<int,int>, std::vector<int>> shortest_paths = {...};  // Pre-computed paths
 * std::vector<double> Ey = {...};  // Scalar values
 *
 * MS_complex_t ms_complex = graph_MS_cx(adj_list, hop_list, core_adj_list, shortest_paths, Ey);
 *
 * // Access basic MS complex information
 * for (int max_vertex : ms_complex.local_maxima) {
 *     const auto& connected_minima = ms_complex.lmax_to_lmin[max_vertex];
 *     // Process connected minima...
 * }
 *
 * // Access cells and their trajectories
 * for (const auto& [cell_key, cell_paths] : ms_complex.cell_trajectories) {
 *     int lmax = cell_key.lmax;
 *     int lmin = cell_key.lmin;
 *     int cell_idx = cell_key.cell_index;
 *
 *     // Access all trajectories in this cell
 *     for (size_t traj_idx : cell_paths) {
 *         const auto& trajectory = ms_complex.unique_trajectories[traj_idx];
 *         // Process trajectory...
 *     }
 * }
 * @endcode
 *
 * @see MS_complex_t
 * @see cell_t
 * @see is_local_maximum()
 * @see is_local_minimum()
 * @see count_subgraph_set_components()
 */
MS_complex_t graph_MS_cx(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<int>>& hop_list,
    const std::vector<std::vector<int>>& core_adj_list,
    const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
    const std::vector<double>& Ey) {

    if (adj_list.size() != Ey.size() || hop_list.size() != Ey.size() ||
        core_adj_list.size() != Ey.size()) {
        Rf_error("Input vector sizes don't match");
    }

    #define DEBUG__graph_MS_cx 0

    MS_complex_t ms_cx;
    std::set<std::vector<int>> unique_trajectories_set;

    // First, precompute all strict local maxima and minima
    for (size_t vertex = 0; vertex < adj_list.size(); vertex++) {
        if (is_local_maximum(vertex, adj_list, Ey.data())) {
            ms_cx.local_maxima.insert(vertex);
        }
        if (is_local_minimum(vertex, adj_list, Ey.data())) {
            ms_cx.local_minima.insert(vertex);
        }
    }

    #if DEBUG__graph_MS_cx
    print_set(ms_cx.local_maxima, "local_maxima");
    print_set(ms_cx.local_minima, "local_minima");
    #endif

    // Helper function to add a trajectory and update MS complex
    auto add_trajectory = [&](const std::vector<int>& ascending,
                              const std::vector<int>& descending) {
        // Construct complete trajectory: descending -> ascending (excluding middle Rf_duplicate vertex)
        std::vector<int> complete_trajectory;
        complete_trajectory.reserve(descending.size() + ascending.size() - 1);

        std::copy(descending.rbegin(), descending.rend() - 1,
                 std::back_inserter(complete_trajectory));
        std::copy(ascending.begin(), ascending.end(),
                 std::back_inserter(complete_trajectory));

        // Check if this is a new unique trajectory
        auto [it, inserted] = unique_trajectories_set.insert(complete_trajectory);

        if (inserted) {
            ms_cx.unique_trajectories.push_back(complete_trajectory);
            size_t traj_idx = ms_cx.unique_trajectories.size() - 1;

            int local_min = complete_trajectory.front();
            if (!ms_cx.local_minima.count(local_min)) {
                Rprintf("\nCRITICAL ERROR: The last vertex of a trajectory is not a local minimum - end vertex: %d!\n", local_min);
                Rf_error("\nFATAL: Found a non-miniumu vertex that is the end of a trajectory.\n");
            }

            int local_max = complete_trajectory.back();
            if (!ms_cx.local_maxima.count(local_max)) {
                Rprintf("\nCRITICAL ERROR: The first vertex of a trajectory is not a local maximum - start vertex: %d!\n", local_max);
                Rf_error("\nFATAL: Found a non-maximum vertex that is the start of a trajectory.\n");
            }

            ms_cx.lmax_to_lmin[local_max].insert(local_min);
            ms_cx.lmin_to_lmax[local_min].insert(local_max);

            std::pair<int,int> cell_key(local_min, local_max);
            ms_cx.procells[cell_key].insert(complete_trajectory.begin(),
                                            complete_trajectory.end());

            // Store trajectory index for later cell assignment
            return std::make_pair(true, traj_idx);
        }
        return std::make_pair(false, size_t(0));
    };

    // Helper function to get valid neighbors with their paths
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<hop_neighbor_info_t> neighbors;
        int min_valid_hop = -1;

        // First pass: find minimum hop distance with valid neighbors
        for (size_t i = 0; i < adj_list[current_vertex].size(); i++) {
            int neighbor = adj_list[current_vertex][i];
            int hop_dist = hop_list[current_vertex][i];

            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0) {
                if (min_valid_hop == -1 || hop_dist < min_valid_hop) {
                    min_valid_hop = hop_dist;
                }
            }
        }

        // Second pass: collect all neighbors at minimum hop distance
        if (min_valid_hop != -1) {
            for (size_t i = 0; i < adj_list[current_vertex].size(); i++) {
                int neighbor = adj_list[current_vertex][i];
                int hop_dist = hop_list[current_vertex][i];

                if (hop_dist == min_valid_hop) {
                    double diff = ascending ?
                        Ey[neighbor] - Ey[current_vertex] :
                        Ey[current_vertex] - Ey[neighbor];

                    if (diff > 0) {
                        auto path_key = std::make_pair(
                            std::min(current_vertex, neighbor),
                            std::max(current_vertex, neighbor)
                            );

                        auto path_it = shortest_paths.find(path_key);
                        if (path_it == shortest_paths.end()) {
                            Rprintf("\nCRITICAL ERROR: No path found between vertices %d and %d in shortest_paths map\n",
                                    current_vertex, neighbor);
                            Rf_error("\nFATAL: Missing shortest path record in shortest_paths map.\n");
                        }

                        const auto& path = path_it->second;
                        if (path.empty()) {
                            Rprintf("\nCRITICAL ERROR: No path exists in core_adj_list between vertices %d and %d\n",
                                    current_vertex, neighbor);
                            Rf_error("\nFATAL: We found a pair of vertices for which the shortest path connecting them is empty.\n");
                        }

                        neighbors.emplace_back(hop_neighbor_info_t{
                                neighbor,
                                hop_dist,
                                diff,
                                path_it->second
                            });
                    }
                }
            }

            if (neighbors.size() > 1) {
                std::sort(neighbors.begin(), neighbors.end(),
                          [](const hop_neighbor_info_t& a, const hop_neighbor_info_t& b) {
                              return a.diff_value > b.diff_value;
                          });
            }
        }

        return neighbors;
    };

    #if DEBUG__graph_MS_cx
    Rprintf("\nEntering the main loop\n");
    #endif
    // Process vertices and compute trajectories
    std::vector<size_t> valid_trajectory_indices;
    for (size_t vertex = 0; vertex < adj_list.size(); vertex++) {
        bool is_max = is_local_maximum(vertex, adj_list, Ey.data());
        bool is_min = is_local_minimum(vertex, adj_list, Ey.data());

        #if DEBUG__graph_MS_cx
        Rprintf("vertex: %d\n", (int)vertex);
        Rprintf("is_max: %s\n", is_max ? "TRUE" : "FALSE");
        Rprintf("is_min: %s\n", is_min ? "TRUE" : "FALSE");
        #endif

        if (is_max) {
            if (!ms_cx.local_maxima.count(vertex)) {
                Rprintf("\nCRITICAL ERROR: We found a new local maximum - vertex: %d!\n", (int)vertex);
                Rf_error("\nFATAL: Found a non-maximum vertex that appears to be a local maximum.\n");
            }
            continue;  // Skip to next vertex
        }

        if (is_min) {
            if (!ms_cx.local_minima.count(vertex)) {
                Rprintf("\nCRITICAL ERROR: We found a new local minimum - vertex: %d!\n", (int)vertex);
                Rf_error("\nFATAL: Found a non-minimum vertex that appears to be a local minimum.\n");
            }
            continue;  // Skip to next vertex
        }

        std::vector<int> ascending = {static_cast<int>(vertex)};
        std::vector<int> descending = {static_cast<int>(vertex)};
        int asc_vertex = vertex;
        int desc_vertex = vertex;

        // Compute trajectories if not at extrema
        if (!is_max) {
            while (true) {

                #if DEBUG__graph_MS_cx
                Rprintf("asc_vertex: %d\n", asc_vertex);
                #endif

                if (ms_cx.local_maxima.count(asc_vertex)) { // Check if we've reached a local maximum
                    ascending.push_back(asc_vertex);
                    break;
                }

                auto valid_neighbors = get_valid_neighbors(asc_vertex, true);
                if (valid_neighbors.empty()) {
                    // This should never happen for non-extremal vertices
                    Rprintf("\nCRITICAL ERROR: No valid ascending neighbors found for non-maximum vertex!\n");
                    Rprintf("Details:\n");
                    Rprintf("  - Current vertex: %d\n", asc_vertex);
                    Rprintf("  - Function value: %.6f\n", Ey[asc_vertex]);
                    Rprintf("  - Trajectory so far: ");
                    for (int v : ascending) {
                        Rprintf("%d ", v);
                    }

                    //Rprintf("\n\nNeighborhood Analysis:\n");
                    Rprintf("  Neighbors in graph:");
                    for (int n : adj_list[asc_vertex]) {
                        Rprintf(" %d(%.6f)", n, Ey[n]);
                    }

                    Rf_error("\nFATAL: Found a non-maximum vertex with no valid ascending neighbors.\n"
                          "This violates the properties of a valid Morse function on the graph.\n"
                          "Please verify that:\n"
                          "1. The function values form a valid Morse function\n"
                          "2. The graph connectivity is correct");
                }

                auto& path = valid_neighbors[0].path;
                ascending.insert(ascending.end(), path.begin() + 1, path.end());
                asc_vertex = valid_neighbors[0].vertex;
            }
        }

        if (!is_min) {
            while (true) {

                #if DEBUG__graph_MS_cx
                Rprintf("desc_vertex: %d\n", desc_vertex);
                #endif

                if (ms_cx.local_minima.count(desc_vertex)) {
                    descending.push_back(desc_vertex);
                    break;
                }

                auto valid_neighbors = get_valid_neighbors(desc_vertex, false);
                if (valid_neighbors.empty()) {
                    // This should never happen for non-extremal vertices
                    Rprintf("\nCRITICAL ERROR: No valid descending neighbors found for non-minimum vertex!\n");
                    Rprintf("Details:\n");
                    Rprintf("  - Current vertex: %d\n", desc_vertex);
                    Rprintf("  - Function value: %.6f\n", Ey[desc_vertex]);
                    Rprintf("  - Trajectory so far: ");
                    for (int v : descending) {
                        Rprintf("%d ", v);
                    }

                    //Rprintf("\n\nNeighborhood Analysis:\n");
                    Rprintf("  Neighbors in graph:");
                    for (int n : adj_list[desc_vertex]) {
                        Rprintf(" %d(%.6f)", n, Ey[n]);
                    }

                    Rf_error("\nFATAL: Found a non-minimum vertex with no valid descending neighbors.\n"
                          "This violates the properties of a valid Morse function on the graph.\n"
                          "Please verify that:\n"
                          "1. The function values form a valid Morse function\n"
                          "2. The graph connectivity is correct");
                }

                auto& path = valid_neighbors[0].path;
                descending.insert(descending.end(), path.begin() + 1, path.end());
                desc_vertex = valid_neighbors[0].vertex;
            }
        }

        // Add trajectory and collect valid indices
        auto [is_valid, traj_idx] = add_trajectory(ascending, descending);
        if (is_valid) {
            valid_trajectory_indices.push_back(traj_idx);
        }
    }

    #if 0
    // Compute MS cells and assign trajectories to cells
    MSCellComputer cell_computer;
    for (const auto& [key, procell] : ms_cx.procells) {
        ms_cx.cells[key] = cell_computer.compute_cells(ms_cx, key);
    }


    for (const auto& [key, procell] : ms_cx.procells) {
        if (procell.size() > 2) {
            // Remove extrema for component calculation
            std::set<int> cell_vertices = procell;
            cell_vertices.erase(key.first);   // Remove local minimum
            cell_vertices.erase(key.second);  // Remove local maximum

            // Find connected components
            auto components = count_subgraph_set_components(core_adj_list, cell_vertices);

            // Create cells and assign trajectories
            std::vector<std::set<int>> cell_components;
            int cell_idx = 0;
            for (const auto& [comp_id, comp_vertices] : components) {
                std::set<int> new_comp = comp_vertices;
                new_comp.insert(key.first);
                new_comp.insert(key.second);
                cell_components.push_back(new_comp);

                // Create cell identifier
                cell_t cell{key.first, key.second, cell_idx};

                // Find trajectories that belong to this cell
                for (size_t traj_idx : valid_trajectory_indices) {
                    const auto& trajectory = ms_cx.unique_trajectories[traj_idx];
                    bool belongs_to_cell = true;

                    // Check if all trajectory vertices are in the cell
                    for (int vertex : trajectory) {
                        if (vertex != key.first && vertex != key.second &&
                            comp_vertices.find(vertex) == comp_vertices.end()) {
                            belongs_to_cell = false;
                            break;
                        }
                    }

                    if (belongs_to_cell) {
                        ms_cx.cell_trajectories[cell].insert(traj_idx);
                    }
                }
                cell_idx++;
            }
            ms_cx.cells[key] = cell_components;
        } else {
            // Handle simple cells (just two vertices)
            std::vector<std::set<int>> cell_components;
            cell_components.push_back(procell);
            ms_cx.cells[key] = cell_components;

            cell_t cell{key.first, key.second, 0};

            // Find direct trajectories between the extrema
            for (size_t traj_idx : valid_trajectory_indices) {
                const auto& trajectory = ms_cx.unique_trajectories[traj_idx];
                if (trajectory.front() == key.first && trajectory.back() == key.second) {
                    ms_cx.cell_trajectories[cell].insert(traj_idx);
                }
            }
        }
    }
    #endif

    return ms_cx;
}
