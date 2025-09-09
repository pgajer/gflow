//
//  Morse-Smale complexe utility functions
//

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>     // for std::vector
#include <numeric>    // for std::iota
#include <algorithm>  // for std::copy, std::find
#include <set>        // for std::set
#include <map>        // for std::map
#include <utility>    // for std::pair

#include "SEXP_cpp_conversion_utils.h"
#include "cpp_utils.h"
#include "MS_complex.h"

std::unique_ptr<std::unordered_map<int, int>> count_subgraph_components(const std::vector<std::vector<int>>& graph,
                                                                        const std::vector<int>& V);
std::vector<int> shortest_path_by_hops(const std::vector<std::vector<int>>& graph,
                                       int start, int end);


extern "C" {
    SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
    SEXP S_graph_MS_cx_with_path_search(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
}


/**
 * @brief Computes connected components of a graph using Union-Find (disjoint set) data structure
 *
 * @details Implements weighted quick union with path compression for optimal performance.
 * Each vertex starts in its own component, and components are merged when edges are processed.
 * The function returns a vector where each entry i contains the representative (root) vertex
 * of the component containing vertex i.
 *
 * Time complexity: O(α(V)×E) where α is the inverse Ackermann function, V is number of vertices,
 * and E is number of edges.
 * Space complexity: O(V)
 *
 * @param adj_vect Adjacency list representation of the graph
 * @return Vector where entry i is the component identifier for vertex i
 */
std::vector<int> union_find(const std::vector<std::vector<int>>& adj_vect) {
    int n_vertices = adj_vect.size();
    std::vector<int> parent(n_vertices);
    std::iota(parent.begin(), parent.end(), 0);

    auto find_root = [&](int i) {
        while (i != parent[i]) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };

    std::vector<int> rank(n_vertices, 0);
    auto union_sets = [&](int x, int y) {
        int root_x = find_root(x);
        int root_y = find_root(y);
        if (root_x != root_y) {
            if (rank[root_x] < rank[root_y]) {
                parent[root_x] = root_y;
            } else if (rank[root_x] > rank[root_y]) {
                parent[root_y] = root_x;
            } else {
                parent[root_y] = root_x;
                rank[root_x]++;
            }
        }
    };

    for (int i = 0; i < n_vertices; ++i) {
        for (int adj : adj_vect[i]) {
            union_sets(i, adj);
        }
    }

    std::vector<int> component_vect(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        component_vect[i] = find_root(i);
    }
    return component_vect;
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
 *    - Formally: pro-cell(M,m) = {v ∈ V | ∃ trajectory through v from m to M}
 *
 * 3. Morse-Smale Cells:
 *    - Connected components of pro-cells after removing their defining extrema
 *    - For pro-cell(M,m), compute components of subgraph induced by:
 *      vertices ∈ pro-cell(M,m) \ {M,m}
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
 * @throws R error if input lengths don't match
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
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);

    // Initialize maps for extrema connectivity
    std::map<int, std::set<int>> lmax_to_lmin_map;
    std::map<int, std::set<int>> lmin_to_lmax_map;
    std::set<int> local_maxima_set;
    std::set<int> local_minima_set;

    // Allocate R objects
    int nprot = 0;
    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices)); nprot++;

    // Helper function for valid neighbors (same as before)
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

    // Add map for pro-cells
    std::map<std::pair<int,int>, std::set<int>> ms_procells;
    std::map<std::pair<int,int>, std::vector<std::set<int>>> ms_cells;

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
                if (valid_neighbors.empty()) break;
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
                if (valid_neighbors.empty()) break;
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
            SEXP trajectory = PROTECT(allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
            // Construct full trajectory
            descending_trajectory.erase(descending_trajectory.begin());
            int tr_size = ascending_trajectory.size() + descending_trajectory.size();
            SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
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

    // Convert maps to R lists
    SEXP lmax_to_lmin       = PROTECT(allocVector(VECSXP, local_maxima_set.size())); nprot++;
    SEXP lmax_to_lmin_names = PROTECT(allocVector(STRSXP, local_maxima_set.size())); nprot++;

    SEXP lmin_to_lmax       = PROTECT(allocVector(VECSXP, local_minima_set.size())); nprot++;
    SEXP lmin_to_lmax_names = PROTECT(allocVector(STRSXP, local_minima_set.size())); nprot++;

    // Convert local maxima/minima sets to R vectors
    SEXP local_maxima = PROTECT(allocVector(INTSXP, local_maxima_set.size())); nprot++;
    SEXP local_minima = PROTECT(allocVector(INTSXP, local_minima_set.size())); nprot++;

    // Fill local maxima vector and its connectivity list
    int max_idx = 0;
    for (int lmax : local_maxima_set) {
        INTEGER(local_maxima)[max_idx] = lmax;

        const auto& connected_mins = lmax_to_lmin_map[lmax];
        SEXP mins = PROTECT(allocVector(INTSXP, connected_mins.size()));
        std::copy(connected_mins.begin(), connected_mins.end(), INTEGER(mins));
        SET_VECTOR_ELT(lmax_to_lmin, max_idx, mins);
        std::string lmax_str = std::to_string(lmax + 1);  // Turning lmax into const char* after turning it into 1-base interger
        const char* lmax_char = lmax_str.c_str();
        SET_STRING_ELT(lmax_to_lmin_names, max_idx, mkChar(lmax_char));
        UNPROTECT(1);

        max_idx++;
    }
    setAttrib(lmax_to_lmin, R_NamesSymbol, lmax_to_lmin_names);

    // Fill local minima vector and its connectivity list
    int min_idx = 0;
    for (int lmin : local_minima_set) {
        INTEGER(local_minima)[min_idx] = lmin;

        const auto& connected_maxs = lmin_to_lmax_map[lmin];
        SEXP maxs = PROTECT(allocVector(INTSXP, connected_maxs.size()));
        std::copy(connected_maxs.begin(), connected_maxs.end(), INTEGER(maxs));
        SET_VECTOR_ELT(lmin_to_lmax, min_idx, maxs);
        std::string lmin_str = std::to_string(lmin + 1);  // Turning lmin into const char*
        const char* lmin_char = lmin_str.c_str();
        SET_STRING_ELT(lmin_to_lmax_names, min_idx, mkChar(lmin_char));
        UNPROTECT(1);

        min_idx++;
    }
    setAttrib(lmin_to_lmax, R_NamesSymbol, lmin_to_lmax_names);

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

    // Convert pro-cells and cells to R structures
    SEXP procells = PROTECT(allocVector(VECSXP, ms_procells.size())); nprot++;
    SEXP cells = PROTECT(allocVector(VECSXP, ms_cells.size())); nprot++;
    SEXP procell_keys = PROTECT(allocVector(VECSXP, ms_procells.size())); nprot++;

    int pcell_idx = 0;
    for (const auto& [key, vertices] : ms_procells) {
        // Create key pair
        SEXP key_pair = PROTECT(allocVector(INTSXP, 2));
        INTEGER(key_pair)[0] = key.first;
        INTEGER(key_pair)[1] = key.second;
        SET_VECTOR_ELT(procell_keys, pcell_idx, key_pair);

        // Create vertex set
        SEXP vertex_set = PROTECT(allocVector(INTSXP, vertices.size()));
        std::copy(vertices.begin(), vertices.end(), INTEGER(vertex_set));
        SET_VECTOR_ELT(procells, pcell_idx, vertex_set);

        // Create cells list for this pro-cell
        const auto& cell_components = ms_cells[key];
        SEXP components_list = PROTECT(allocVector(VECSXP, cell_components.size()));

        for (size_t i = 0; i < cell_components.size(); ++i) {
            SEXP component = PROTECT(allocVector(INTSXP, cell_components[i].size()));
            std::copy(cell_components[i].begin(), cell_components[i].end(),
                     INTEGER(component));
            SET_VECTOR_ELT(components_list, i, component);
            UNPROTECT(1);
        }

        SET_VECTOR_ELT(cells, pcell_idx, components_list);

        UNPROTECT(3);  // key_pair, vertex_set, components_list
        pcell_idx++;
    }

    // Update result list to include pro-cells and cells
    SEXP result = PROTECT(allocVector(VECSXP, 8)); nprot++;
    SET_VECTOR_ELT(result, 0, trajectories);
    SET_VECTOR_ELT(result, 1, lmax_to_lmin);
    SET_VECTOR_ELT(result, 2, lmin_to_lmax);
    SET_VECTOR_ELT(result, 3, local_maxima);
    SET_VECTOR_ELT(result, 4, local_minima);
    SET_VECTOR_ELT(result, 5, procell_keys);
    SET_VECTOR_ELT(result, 6, procells);
    SET_VECTOR_ELT(result, 7, cells);

    // Update names
    SEXP names = PROTECT(allocVector(STRSXP, 8)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("trajectories"));
    SET_STRING_ELT(names, 1, mkChar("lmax_to_lmin"));
    SET_STRING_ELT(names, 2, mkChar("lmin_to_lmax"));
    SET_STRING_ELT(names, 3, mkChar("local_maxima"));
    SET_STRING_ELT(names, 4, mkChar("local_minima"));
    SET_STRING_ELT(names, 5, mkChar("procell_keys"));
    SET_STRING_ELT(names, 6, mkChar("procells"));
    SET_STRING_ELT(names, 7, mkChar("cells"));
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(nprot);
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
 * @throws R error if input lengths don't match
 */
SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey) {
    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);

    SEXP trajectories = PROTECT(allocVector(VECSXP, n_vertices));

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
            SEXP trajectory = PROTECT(allocVector(INTSXP, 1));
            INTEGER(trajectory)[0] = vertex;
            SET_VECTOR_ELT(trajectories, vertex, trajectory);
            UNPROTECT(1);
        } else {
            // Construct full trajectory
            descending_trajectory.erase(descending_trajectory.begin());
            int tr_size = ascending_trajectory.size() + descending_trajectory.size();
            SEXP trajectory = PROTECT(allocVector(INTSXP, tr_size));
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
 * @note The function assumes that shortest_paths contains valid paths between connected
 *       vertices in the graph. Paths should be monotonic with respect to the scalar values
 *       in Ey.
 *
 * @details The algorithm proceeds in several steps:
 * 1. For each vertex:
 *    - Determines if it is a local maximum or minimum
 *    - Computes a complete gradient trajectory by:
 *      * Following steepest ascending path to a local maximum
 *      * Following steepest descending path to a local minimum
 *    - Updates the connectivity maps between extrema
 *    - Constructs proto-cells from complete trajectories
 *
 * 2. For each proto-cell:
 *    - Removes the extrema vertices
 *    - Identifies connected components in the remaining vertices
 *    - Creates separate cells for each component
 *    - Associates relevant trajectories with each cell
 *
 * The monotonicity of paths ensures that:
 * - Ascending paths strictly increase in scalar value
 * - Descending paths strictly decrease in scalar value
 *
 * @warning The function assumes that:
 * - All vertex indices are valid (0 to adj_list.size()-1)
 * - The graph is undirected
 * - Ey contains valid scalar values for all vertices
 * - shortest_paths contains valid paths when needed
 *
 * Example usage:
 * @code
 * std::vector<std::vector<int>> adj_list = {...};  // Graph adjacency list
 * std::vector<std::vector<int>> core_adj_list = {...};  // Core graph adjacency list
 * std::map<std::pair<int,int>, std::vector<int>> shortest_paths = {...};  // Pre-computed paths
 * std::vector<double> Ey = {...};  // Scalar values
 *
 * MS_complex_t ms_complex = graph_MS_cx(adj_list, core_adj_list, shortest_paths, Ey);
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
 * @see is_path_monotonic_increasing()
 * @see is_path_monotonic_decreasing()
 * @see count_subgraph_set_components()
 */
MS_complex_t graph_MS_cx(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<int>>& core_adj_list,
    const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
    const std::vector<double>& Ey) {

    MS_complex_t ms_cx;
    std::set<std::vector<int>> unique_trajectories_set;

    // Helper function to add a trajectory and update MS complex
    auto add_trajectory = [&](const std::vector<int>& ascending,
                            const std::vector<int>& descending) {
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

            int local_max = complete_trajectory.front();
            int local_min = complete_trajectory.back();

            ms_cx.lmax_to_lmin[local_max].insert(local_min);
            ms_cx.lmin_to_lmax[local_min].insert(local_max);

            std::pair<int,int> cell_key(local_max, local_min);
            ms_cx.procells[cell_key].insert(complete_trajectory.begin(),
                                            complete_trajectory.end());

            // Store trajectory index for later cell assignment
            return std::make_pair(true, traj_idx);
        }
        return std::make_pair(false, size_t(0));
    };

    // Helper function to get valid neighbors with their paths
    auto get_valid_neighbors = [&](int current_vertex, bool ascending) {
        std::vector<neighbor_info_t> neighbors;

        for (int neighbor : adj_list[current_vertex]) {
            // Use pre-computed shortest paths
            auto path_key = std::make_pair(current_vertex, neighbor);
            auto path_it = shortest_paths.find(path_key);
            if (path_it == shortest_paths.end()) continue;

            const auto& path = path_it->second;
            if (path.empty()) continue;

            double diff = ascending ?
                Ey[neighbor] - Ey[current_vertex] :
                Ey[current_vertex] - Ey[neighbor];

            if (diff > 0 && (ascending ?
                             is_path_monotonic_increasing(path, Ey) :
                             is_path_monotonic_decreasing(path, Ey))) {
                neighbors.emplace_back(neighbor_info_t{
                        neighbor,
                        path,
                        static_cast<int>(path.size() - 1),
                        diff
                    });
            }
        }

        // Sort by difference value in descending order
        std::sort(neighbors.begin(), neighbors.end(),
                  [](const neighbor_info_t& a, const neighbor_info_t& b) {
                      return a.diff_value > b.diff_value;
                  });

        return neighbors;
    };

    // Process vertices and compute trajectories
    std::vector<size_t> valid_trajectory_indices;
    for (size_t vertex = 0; vertex < adj_list.size(); vertex++) {
        bool is_max = is_local_maximum(vertex, adj_list, Ey.data());
        bool is_min = is_local_minimum(vertex, adj_list, Ey.data());

        if (is_max) ms_cx.local_maxima.insert(vertex);
        if (is_min) ms_cx.local_minima.insert(vertex);

        std::vector<int> ascending = {static_cast<int>(vertex)};
        std::vector<int> descending = {static_cast<int>(vertex)};
        int asc_vertex = vertex;
        int desc_vertex = vertex;

        // Compute trajectories if not at extrema
        if (!is_max) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(asc_vertex, true);
                if (valid_neighbors.empty()) break;
                auto& path = valid_neighbors[0].path;
                ascending.insert(ascending.end(), path.begin() + 1, path.end());
                asc_vertex = valid_neighbors[0].vertex;
            }

            if (is_local_maximum(asc_vertex, adj_list, Ey.data())) {
                ms_cx.local_maxima.insert(asc_vertex);
            }
        }

        if (!is_min) {
            while (true) {
                auto valid_neighbors = get_valid_neighbors(desc_vertex, false);
                if (valid_neighbors.empty()) break;
                auto& path = valid_neighbors[0].path;
                descending.insert(descending.end(), path.begin() + 1, path.end());
                desc_vertex = valid_neighbors[0].vertex;
            }

            if (is_local_minimum(desc_vertex, adj_list, Ey.data())) {
                ms_cx.local_minima.insert(desc_vertex);
            }
        }

        // Add trajectory and collect valid indices
        auto [is_valid, traj_idx] = add_trajectory(ascending, descending);
        if (is_valid) {
            valid_trajectory_indices.push_back(traj_idx);
        }
    }

    // Compute MS cells and assign trajectories to cells
    for (const auto& [key, procell] : ms_cx.procells) {
        if (procell.size() > 2) {
            // Remove extrema for component calculation
            std::set<int> cell_vertices = procell;
            cell_vertices.erase(key.first);   // Remove local maximum
            cell_vertices.erase(key.second);  // Remove local minimum

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

    return ms_cx;
}

// I don't think I would ever your this function in isolation from path graph construct and kernel mean estimates
#if 0
SEXP S_graph_MS_cx(SEXP s_graph,
                   SEXP s_core_graph,
                   SEXP s_shortest_paths,
                   SEXP s_Ey) {

    int n_vertices = LENGTH(s_graph);

    if (LENGTH(s_Ey) != n_vertices || LENGTH(s_core_graph) != n_vertices) {
        error("The lengths of s_graph, s_core_graph, and s_Ey must be the same");
    }

    // Convert input data
    std::vector<std::vector<int>> graph = Rgraph_to_vector(s_graph);
    std::vector<std::vector<int>> core_graph = Rgraph_to_vector(s_core_graph);
    double* Ey = REAL(s_Ey);


}
#endif

/**
 * @brief Computes the Morse-Smale complex with extensive validation checks
 *
 * This function is a validated version of graph_MS_cx that performs extensive checks
 * to ensure the correctness of the computed Morse-Smale complex. It includes validation
 * of trajectories, cell structure, and topological relationships while maintaining the
 * same functionality as graph_MS_cx.
 *
 * @param adj_list The adjacency list representation of the graph. adj_list[i] contains
 *                 the vertices adjacent to vertex i
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
 * @throws std::runtime_error in the following cases:
 *         - Empty trajectory detected
 *         - Non-monotonic trajectory path found
 *         - Trajectory endpoints are not proper extrema
 *         - Disconnected trajectory vertices found
 *         - Cell missing local maximum or minimum
 *         - Disconnected cell detected
 *         - Cell without valid trajectories found
 *         - Invalid simple cell size (not equal to 2)
 *
 * @note This function performs the same computation as graph_MS_cx but with additional
 *       validation checks. Use this version when:
 *       - Debugging complex issues in Morse-Smale complex computation
 *       - Working with potentially problematic or noisy input data
 *       - Needing guaranteed validity of the output structure
 *       - During development and testing phases
 *
 * @warning The validation checks may significantly impact performance. For production
 *          use with well-tested data, consider using graph_MS_cx instead.
 *
 * @pre The following preconditions are checked:
 *      - Non-empty input graph (adj_list and core_adj_list)
 *      - Valid vertex indices in shortest_paths
 *      - Monotonic paths in shortest_paths
 *      - Connected graph structure in core_adj_list
 *
 * @post The following postconditions are guaranteed:
 *      - All trajectories are monotonic and connect proper extrema
 *      - All cells are connected and contain their extrema
 *      - Each cell has at least one valid trajectory
 *      - Cell decomposition preserves topological structure
 *
 * @details The validation process includes:
 *
 * 1. Trajectory Validation:
 *    - Non-empty check for all paths
 *    - Strict monotonicity verification along paths
 *    - Proper extrema at trajectory endpoints
 *    - Path connectivity in core_adj_list
 *
 * 2. Cell Structure Validation:
 *    - Presence of extrema in each cell
 *    - Cell connectivity using depth-first search
 *    - Existence of valid trajectories
 *    - Proper decomposition of procells
 *
 * 3. Topological Validation:
 *    - Consistent extrema connections
 *    - Valid cell-trajectory associations
 *    - Proper simple cell structure
 *
 * Example usage:
 * @code
 * try {
 *     std::vector<std::vector<int>> adj_list = {...};  // Graph adjacency list
 *     std::vector<std::vector<int>> core_adj_list = {...};  // Core graph adjacency list
 *     std::map<std::pair<int,int>, std::vector<int>> shortest_paths = {...};  // Pre-computed paths
 *     std::vector<double> Ey = {...};  // Scalar values
 *
 *     MS_complex_t ms_complex = graph_MS_cx_with_validation(
 *         adj_list, core_adj_list, shortest_paths, Ey
 *     );
 *
 *     // Process validated MS complex...
 *
 * } catch (const std::runtime_error& e) {
 *     std::cerr << "MS complex validation failed: " << e.what() << std::endl;
 *     // Handle validation failure...
 * }
 * @endcode
 *
 * Performance considerations:
 * - The validation checks add significant computational overhead
 * - Memory usage is identical to graph_MS_cx
 * - Consider using graph_MS_cx for performance-critical applications
 *
 * @see graph_MS_cx() The non-validated version of this function
 * @see MS_complex_t The structure containing the computed complex
 * @see cell_t The structure identifying individual cells
 * @see is_local_maximum()
 * @see is_local_minimum()
 * @see is_path_monotonic_increasing()
 * @see is_path_monotonic_decreasing()
 * @see count_subgraph_set_components()
 *
 * @since Version 2.0
 * @author [Your Name]
 * @date [Current Date]
 */

#if 0
MS_complex_t graph_MS_cx_with_validation(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<int>>& core_adj_list,
    const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
    const std::vector<double>& Ey) {

    MS_complex_t ms_cx;
    std::set<std::vector<int>> unique_trajectories_set;

    // Helper function to validate trajectory
    auto validate_trajectory = [&](const std::vector<int>& trajectory, bool check_extrema = true) {
        if (trajectory.empty()) {
            throw std::runtime_error("Empty trajectory detected");
        }

        // Check monotonicity
        for (size_t i = 0; i < trajectory.size() - 1; ++i) {
            if (Ey[trajectory[i + 1]] <= Ey[trajectory[i]]) {
                throw std::runtime_error(
                    "Non-monotonic trajectory detected: " +
                    std::to_string(trajectory[i]) + "(" + std::to_string(Ey[trajectory[i]]) + ") -> " +
                    std::to_string(trajectory[i + 1]) + "(" + std::to_string(Ey[trajectory[i + 1]]) + ")"
                );
            }
        }

        if (check_extrema) {
            // Check if endpoints are proper extrema
            if (!is_local_maximum(trajectory.front(), adj_list, Ey.data())) {
                throw std::runtime_error(
                    "Trajectory doesn't start at local maximum: " +
                    std::to_string(trajectory.front())
                );
            }
            if (!is_local_minimum(trajectory.back(), adj_list, Ey.data())) {
                throw std::runtime_error(
                    "Trajectory doesn't end at local minimum: " +
                    std::to_string(trajectory.back())
                );
            }
        }

        // Check path connectivity in core_adj_list
        for (size_t i = 0; i < trajectory.size() - 1; ++i) {
            bool connected = false;
            for (int neighbor : core_adj_list[trajectory[i]]) {
                if (neighbor == trajectory[i + 1]) {
                    connected = true;
                    break;
                }
            }
            if (!connected) {
                throw std::runtime_error(
                    "Disconnected trajectory vertices: " +
                    std::to_string(trajectory[i]) + " -> " +
                    std::to_string(trajectory[i + 1])
                );
            }
        }
    };

    // Helper function to add a trajectory and update MS complex
    auto add_trajectory = [&](const std::vector<int>& ascending,
                            const std::vector<int>& descending) {
        // Validate individual paths before merging
        validate_trajectory(ascending, false);  // Don't check extrema for half-paths
        validate_trajectory(descending, false);

        std::vector<int> complete_trajectory;
        complete_trajectory.reserve(descending.size() + ascending.size() - 1);

        // Create complete trajectory
        std::copy(descending.rbegin(), descending.rend() - 1,
                 std::back_inserter(complete_trajectory));
        std::copy(ascending.begin(), ascending.end(),
                 std::back_inserter(complete_trajectory));

        // Validate complete trajectory
        validate_trajectory(complete_trajectory);

        // Check if this is a new unique trajectory
        auto [it, inserted] = unique_trajectories_set.insert(complete_trajectory);

        if (inserted) {
            size_t traj_idx = ms_cx.unique_trajectories.size();
            ms_cx.unique_trajectories.push_back(complete_trajectory);

            int local_max = complete_trajectory.front();
            int local_min = complete_trajectory.back();

            if (is_local_minimum(local_min, adj_list, Ey.data()) &&
                is_local_maximum(local_max, adj_list, Ey.data())) {

                ms_cx.lmax_to_lmin[local_max].insert(local_min);
                ms_cx.lmin_to_lmax[local_min].insert(local_max);

                std::pair<int,int> cell_key(local_max, local_min);
                ms_cx.procells[cell_key].insert(complete_trajectory.begin(),
                                              complete_trajectory.end());

                return std::make_pair(true, traj_idx);
            }
        }
        return std::make_pair(false, size_t(0));
    };

    // ... rest of the function remains the same until cell computation ...

    // Helper function to validate cell structure
    auto validate_cell = [&](const std::set<int>& cell_vertices,
                           int local_max,
                           int local_min) {
        // Check extrema presence
        if (cell_vertices.find(local_max) == cell_vertices.end()) {
            throw std::runtime_error(
                "Cell missing local maximum: " + std::to_string(local_max)
            );
        }
        if (cell_vertices.find(local_min) == cell_vertices.end()) {
            throw std::runtime_error(
                "Cell missing local minimum: " + std::to_string(local_min)
            );
        }

        // Check cell connectivity
        std::set<int> visited;
        std::function<void(int)> dfs = [&](int vertex) {
            visited.insert(vertex);
            for (int neighbor : core_adj_list[vertex]) {
                if (cell_vertices.find(neighbor) != cell_vertices.end() &&
                    visited.find(neighbor) == visited.end()) {
                    dfs(neighbor);
                }
            }
        };
        dfs(local_max);

        if (visited != cell_vertices) {
            throw std::runtime_error("Disconnected cell detected");
        }
    };

    // Compute MS cells and assign trajectories to cells
    for (const auto& [key, procell] : ms_cx.procells) {
        if (procell.size() > 2) {
            std::set<int> cell_vertices = procell;
            cell_vertices.erase(key.first);   // Remove local maximum
            cell_vertices.erase(key.second);  // Remove local minimum

            auto components = count_subgraph_set_components(core_adj_list, cell_vertices);

            std::vector<std::set<int>> cell_components;
            int cell_idx = 0;
            for (const auto& [comp_id, comp_vertices] : components) {
                std::set<int> new_comp = comp_vertices;
                new_comp.insert(key.first);
                new_comp.insert(key.second);

                // Validate cell structure
                validate_cell(new_comp, key.first, key.second);

                cell_components.push_back(new_comp);
                cell_t cell{key.first, key.second, cell_idx};

                // Track cell validity
                bool has_valid_trajectory = false;

                // Find trajectories that belong to this cell
                for (size_t traj_idx : valid_trajectory_indices) {
                    const auto& trajectory = ms_cx.unique_trajectories[traj_idx];
                    bool belongs_to_cell = true;

                    for (int vertex : trajectory) {
                        if (vertex != key.first && vertex != key.second &&
                            comp_vertices.find(vertex) == comp_vertices.end()) {
                            belongs_to_cell = false;
                            break;
                        }
                    }

                    if (belongs_to_cell) {
                        ms_cx.cell_trajectories[cell].insert(traj_idx);
                        has_valid_trajectory = true;
                    }
                }

                // Each cell should have at least one trajectory
                if (!has_valid_trajectory) {
                    throw std::runtime_error(
                        "Cell without valid trajectories detected: max=" +
                        std::to_string(key.first) + ", min=" +
                        std::to_string(key.second) + ", idx=" +
                        std::to_string(cell_idx)
                    );
                }

                cell_idx++;
            }
            ms_cx.cells[key] = cell_components;
        } else {
            // Validate simple cells
            if (procell.size() != 2) {
                throw std::runtime_error("Invalid simple cell size");
            }

            std::vector<std::set<int>> cell_components;
            cell_components.push_back(procell);
            ms_cx.cells[key] = cell_components;

            cell_t cell{key.first, key.second, 0};
            bool has_valid_trajectory = false;

            for (size_t traj_idx : valid_trajectory_indices) {
                const auto& trajectory = ms_cx.unique_trajectories[traj_idx];
                if (trajectory.front() == key.first && trajectory.back() == key.second) {
                    ms_cx.cell_trajectories[cell].insert(traj_idx);
                    has_valid_trajectory = true;
                }
            }

            if (!has_valid_trajectory) {
                throw std::runtime_error(
                    "Simple cell without valid trajectory detected: max=" +
                    std::to_string(key.first) + ", min=" +
                    std::to_string(key.second)
                );
            }
        }
    }

    return ms_cx;
}
#endif
