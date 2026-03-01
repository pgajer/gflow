#include "fns_over_graphs_r.h"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_diffusion_smoother.hpp"
#include "stats_utils.h"

#include <vector>
#include <queue>
#include <memory>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <stack>
#include <numeric>
#include <unordered_map>
#include <random>
#include <chrono>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

std::unique_ptr<std::unordered_map<int, int>> count_subgraph_components(const std::vector<std::vector<int>>& graph,
                                                                        const std::vector<int>& V);

extern "C" {
    SEXP S_loc_const_vertices(SEXP Rgraph, SEXP Ry, SEXP Rprec);
}

// ------------------------------------------------------------------------------------------
//
// proportion of neighbors with smaller y values for each vertex in a graph
//
// ------------------------------------------------------------------------------------------

/**
 * @brief Calculate the proportion of neighbors with smaller y values for each vertex in a graph.
 *
 * This function computes, for each vertex in the graph, the proportion of its neighbors
 * that have a smaller y value than the vertex itself.
 *
 * @param adj_list A vector of vectors representing the adjacency list of the graph.
 *              adj_list[i] contains the indices of the neighbors of vertex i.
 * @param y A vector of double values associated with each vertex in the adj_list.
 *
 * @return A unique pointer to a vector of doubles. Each element i in the returned vector
 *         represents the proportion of neighbors of vertex i that have a smaller y value
 *         than y[i]. If a vertex has no neighbors, its proportion is set to -1.0.
 *
 * @note The function assumes that the size of y matches the number of vertices in the graph,
 *       and that all neighbor indices in graph are valid (i.e., less than y.size()).
 *
 * @throws std::out_of_range If any neighbor index in graph is out of range for y.
 *
 * Time Complexity: O(E), where E is the total number of edges in the graph.
 * Space Complexity: O(V), where V is the number of vertices in the graph.
 */
std::unique_ptr<std::vector<double>> prop_nbhrs_with_smaller_y(const std::vector<std::vector<int>>& adj_list,
                                                               const std::vector<double>& y ) {
    int n_vertices = y.size();
    auto prop_nbhrs = std::make_unique<std::vector<double>>(n_vertices, -1.0);
    for (int vertex = 0; vertex < n_vertices; vertex++) {
        if (!adj_list[vertex].empty()) {  // Check if the vertex has any neighbors
            (*prop_nbhrs)[vertex] = 0;
            for (const auto& neighbor : adj_list[vertex]) {
                if (y[neighbor] < y[vertex]) {
                    (*prop_nbhrs)[vertex]++;
                }
            }
            (*prop_nbhrs)[vertex] /= adj_list[vertex].size();
        }
        // If adj_list[vertex] is empty, prop_nbhrs[vertex] remains -1.0
    }
    return prop_nbhrs;
}

/**
 * @brief R interface for calculating the proportion of neighbors with smaller y values for each vertex in a graph.
 *
 * This function serves as an interface between R and the C++ prop_nbhrs_with_smaller_y function.
 * It computes, for each vertex in the graph, the proportion of its neighbors that have a smaller y value
 * than the vertex itself.
 *
 * @param r_adj_list An R list representing the adjacency list of the graph. Each element of the list
 *               should be an integer vector containing the indices of neighboring vertices.
 *               The indices in R are assumed to be 1-based and are converted to 0-based for C++.
 * @param Ry An R numeric vector of values associated with each vertex in the graph.
 *
 * @return An R numeric vector (REALSXP) where each element i represents the proportion of neighbors
 *         of vertex i that have a smaller y value than y[i]. If a vertex has no neighbors,
 *         its proportion is set to 0.0.
 *
 * @note This function uses PROTECT/UNPROTECT for proper memory management in R.
 *       It assumes that the length of Ry matches the number of vertices in the graph,
 *       and that all neighbor indices in r_adj_list are valid (i.e., not greater than Rf_length(Ry)).
 *
 * @see prop_nbhrs_with_smaller_y (the underlying C++ function)
 *
 * @example
 * In R:
 * result <- .Call("S_prop_nbhrs_with_smaller_y", graph, y)
 * where graph is a list of integer vectors and y is a numeric vector.
 *
 * @throws R_ERROR If memory allocation fails or if invalid data is provided.
 *
 * Time Complexity: O(E), where E is the total number of edges in the graph.
 * Space Complexity: O(V), where V is the number of vertices in the graph.
 */
SEXP S_prop_nbhrs_with_smaller_y(SEXP r_adj_list, SEXP Ry) {
    // Convert R graph to C++ vector of vectors
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(r_adj_list);

    // Convert R y to C++ vector
    int y_length = LENGTH(Ry);
    std::vector<double> y_vec(REAL(Ry), REAL(Ry) + y_length);

    // Call the C++ function
    std::unique_ptr<std::vector<double>> result = prop_nbhrs_with_smaller_y(adj_list, y_vec);

    // Convert result to SEXP and return
    SEXP Rresult = PROTECT(Rf_allocVector(REALSXP, result->size()));
    for (size_t i = 0; i < result->size(); ++i) {
        REAL(Rresult)[i] = (*result)[i];
    }
    UNPROTECT(1);
    return Rresult;
}


// ------------------------------------------------------------------------------------------
//
// calculate_smallest_difference
//
// ------------------------------------------------------------------------------------------

/**
 * Calculates the smallest difference between any two adjacent unique values in a vector.
 *
 * This function takes a vector of doubles and finds the smallest difference between any two
 * adjacent unique values in the sorted order of the unique values. If the vector contains
 * less than two unique values, the function returns positive infinity.
 *
 * @param y The input vector of doubles.
 * @return The smallest difference between any two adjacent unique values in the vector.
 *         If the vector contains less than two unique values, positive infinity is returned.
 *
 * @note The time complexity of this function is O(n log n), where n is the size of the input vector,
 *       due to the sorting operation performed by std::set. The space complexity is O(n) to store
 *       the unique values in the set.
 *
 * @example
 * std::vector<double> y = {1.5, 2.3, 0.8, 1.5, 3.1, 2.3, 4.2};
 * double smallest_diff = calculate_smallest_difference(y);
 * // smallest_diff will be 0.8, which is the smallest difference between 0.8 and 1.5 in the sorted order.
 */
double calculate_smallest_difference(const std::vector<double>& y) {
    std::set<double> unique_values(y.begin(), y.end());

    if (unique_values.size() < 2) {
        return std::numeric_limits<double>::infinity();
    }

    double min_difference = std::numeric_limits<double>::infinity();
    auto it = unique_values.begin();
    double prev = *it;
    ++it;

    while (it != unique_values.end()) {
        min_difference = std::min(min_difference, *it - prev);
        prev = *it;
        ++it;
    }

    return min_difference;
}


// ------------------------------------------------------------------------------------------
//
// Function for detecting vertices at which a function over a graph is locally constant
//
// ------------------------------------------------------------------------------------------

/**
 * Detects vertices at which a function over a graph is locally constant.
 *
 * A vertex is considered locally constant if the absolute difference between
 * the value of the graph function at the vertex and the value of the graph
 * function at its every neighbor is less than or equal to the specified
 * precision.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference
 *             between the value at a vertex and any of its neighbors is greater than this threshold,
 *             the vertex is not considered locally constant.
 *
 * @return A unique pointer to a vector containing the indices of the vertices
 * at which the function is locally constant. If none of the vertices satisfy
 * the condition of being locally constant (i.e., the absolute difference
 * between the value at each vertex and its neighbors is always greater than the
 * precision), then the function returns a unique pointer to an empty vector,
 * not a null pointer.
 */
std::unique_ptr<std::vector<int>> loc_const_vertices(const std::vector<std::vector<int>>& adj_list,
                                                     const std::vector<double>& y,
                                                     double prec = 1e-8) {
    auto non_const_vertices = std::vector<int>(); // vector of indices of vertices at which y is not constant
    int n_vertices = adj_list.size();

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        for (int neighbor : adj_list[vertex]) {
            if (fabs(y[vertex] - y[neighbor]) > prec) {
                non_const_vertices.push_back(vertex);
                break;
            }
        }
    }

    // Computing the set difference between all indices and non_const_vertices
    std::vector<int> all_indices(n_vertices);
    std::iota(all_indices.begin(), all_indices.end(), 0);

    auto const_vertices = std::make_unique<std::vector<int>>(); // pointer to a vector of indices of vertices at which y is constant

    std::set_difference(all_indices.begin(), all_indices.end(),
                        non_const_vertices.begin(), non_const_vertices.end(),
                        std::back_inserter(*const_vertices));

    return const_vertices;
}

/**
 * An R interface to the loc_const_vertices function.
 *
 * This function takes an R list representing a graph adjacency list, a numeric vector
 * of function values at each vertex, and a precision threshold. It calls the C++
 * loc_const_vertices function to detect vertices at which the function is locally constant
 * and returns the indices of those vertices as an R integer vector.
 *
 * @param r_adj_list An R list representing the graph adjacency list. Each element of the list
 *               should be an integer vector specifying the neighboring vertex indices for
 *               a given vertex.
 * @param Ry A numeric vector representing the function values at each vertex of the graph.
 * @param Rprec A numeric scalar specifying the precision threshold used to determine local
 *              constancy. If the absolute difference between the value at a vertex and any
 *              of its neighbors is greater than this threshold, the vertex is not considered
 *              locally constant.
 * @return An R integer vector containing the indices of the vertices at which the function
 *         is locally constant.
 */
SEXP S_loc_const_vertices(SEXP r_adj_list, SEXP Ry, SEXP Rprec) {

    // ---- Validate inputs to avoid UB when accessing REAL()/INTEGER() ----
    if (!Rf_isVectorList(r_adj_list)) {
        Rf_error("`adj.list` must be a list (of integer vectors).");
    }
    if (!Rf_isNumeric(Ry)) {
        Rf_error("`y` must be a numeric vector.");
    }
    if (!Rf_isNumeric(Rprec) || Rf_length(Rprec) < 1) {
        Rf_error("`prec` must be a numeric scalar.");
    }

    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(r_adj_list);
    const size_t ny = (size_t) LENGTH(Ry);
    std::vector<double> y(REAL(Ry), REAL(Ry) + ny);
    // Coerce scalar and check positivity as the R wrapper enforces
    double prec = Rf_asReal(Rprec);
    if (!(prec > 0)) {
        Rf_error("`prec` must be positive.");
    }

    if ((size_t)ny != adj_list.size()) {
        Rf_error("Length of `y` (%ld) must equal number of vertices in the graph (%zu).",
                 (long)ny, adj_list.size());
    }

    // Compute
    std::unique_ptr<std::vector<int>> locs;
    try {
        locs = loc_const_vertices(adj_list, y, prec);
    } catch (const std::exception& e) {
        Rf_error("loc_const_vertices() failed: %s", e.what());
    } catch (...) {
        Rf_error("loc_const_vertices() failed with an unknown error.");
    }
    if (!locs) {
        Rf_error("Internal error: computation returned null result.");
    }

    const int nres = (int)locs->size(); // safe under LENGTH-first
    SEXP result = PROTECT(Rf_allocVector(INTSXP, nres));
    std::copy(locs->begin(), locs->end(), INTEGER(result));

    UNPROTECT(1);
    return result;
}

/**
 * Finds edges in a graph where the function values at the edge vertices are the same.
 *
 * This function takes a graph represented as an adjacency list and a vector of function values
 * at each vertex. It identifies the edges where the function values at the connected vertices
 * are equal and returns a set of these edges.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph.
 * @return A unique pointer to a set of edges (represented as pairs or sorted sets of integers)
 *         where the function values at the connected vertices are the same.
 *
 * @note The function assumes that the graph is undirected, meaning that if vertex u is connected
 *       to vertex v, then vertex v is also connected to vertex u. The function only considers
 *       edges in one direction to avoid duplicates in the result set.
 */
std::unique_ptr<std::set<std::pair<int, int>>> loc_const_edges(const std::vector<std::vector<int>>& graph,
                                                               const std::vector<double>& y) {
    auto const_edges = std::make_unique<std::set<std::pair<int, int>>>();
    int n_vertices = graph.size();

    for (int vertex = 0; vertex < n_vertices; ++vertex) {
        for (int neighbor : graph[vertex]) {
            if (vertex < neighbor && y[vertex] == y[neighbor]) {
                const_edges->insert(std::make_pair(vertex, neighbor));
            }
        }
    }

    return const_edges;
}

// ------------------------------------------------------------------------------------------
//
// Making the response variable locally non-constant
//
// ------------------------------------------------------------------------------------------

/**
 * Checks if there exists a vertex in the graph at which the function y is locally constant.
 *
 * A vertex is considered locally constant if the absolute difference between the value of y
 * at that vertex and the values of y at its neighboring vertices is within the specified
 * precision threshold.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference
 *             between the value of y at a vertex and any of its neighbors is greater than this threshold,
 *             the vertex is not considered locally constant.
 * @return true if there exists a vertex in the graph at which y is locally constant, false otherwise.
 *
 * @note The time complexity of this function is O(V + E), where V is the number of vertices and E is the
 *       number of edges in the graph. It iterates over each vertex and its neighbors until a locally constant
 *       vertex is found or all vertices have been visited.
 */
bool has_constant_vertex(const std::vector<std::vector<int>>& graph,
                         const std::vector<double>& y,
                         double prec) {

    int n_vertices = y.size();
    bool found_const_vertex = false;
    int vertex = 0;
    while (!found_const_vertex && vertex < n_vertices) {

        bool is_non_const_vertex = false;
        for (int neighbor : graph[vertex]) {
            if (fabs(y[vertex] - y[neighbor]) > prec) {
                is_non_const_vertex = true;
                break;
            }
        }

        if (!is_non_const_vertex)
            found_const_vertex = true;

        ++vertex;
    }

    return found_const_vertex;
}

// overloaded version
bool has_constant_vertex(const std::vector<std::vector<int>>& graph,
                         const std::vector<double>& y,
                         double prec,
                         const std::vector<int>& indices_to_check) {
    bool found_const_vertex = false;

    // Loop over vertices in indices_to_check instead of all vertices
    for (int vertex : indices_to_check) {

        bool is_non_const_vertex = false;
        for (int neighbor : graph[vertex]) {
            if (fabs(y[vertex] - y[neighbor]) > prec) {
                is_non_const_vertex = true;
                break;
            }
        }

        if (!is_non_const_vertex) {
            found_const_vertex = true;
            break;
        }
    }

    return found_const_vertex;
}

/**
 * Makes a function over graph vertices locally non-constant using diffusion smoothing.
 *
 * This function applies diffusion smoothing to a graph to make the response (function values) locally non-constant
 * at each vertex. It iteratively smooths the function values at each vertex based on the weighted average of its
 * neighboring vertices until no locally constant vertices remain or a maximum number of iterations is reached.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector represents a vertex,
 *              and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the initial function values at each vertex of the graph.
 * @param weights A vector of vertex weights used in the weighted smoothing process. The length of `weights` must be
 *                equal to the number of vertices in the graph.
 * @param step_factor A positive number used to calculate the step size. The step size is the product of `step_factor`
 *                    and the smallest difference between the unique values of `y`.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference between the value
 *             at a vertex and any of its neighbors is greater than this threshold, the vertex is not considered
 *             locally constant. Default value is 1e-8.
 * @param n_itrs The maximum number of iterations for the smoothing process. Default value is 1000.
 * @param mean_adjust A boolean indicating whether to adjust the mean of the smoothed function values at each iteration
 *                    to match the mean of the initial function values. If true, the mean adjustment is performed;
 *                    otherwise, it is skipped. Default value is false.
 * @return A unique pointer to a vector of vectors representing the trajectory of the smoothed function values over time.
 *         Each inner vector corresponds to the function values at a specific iteration of the smoothing process.
 */
std::unique_ptr<std::vector<std::vector<double>>> make_response_locally_non_const(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y,
                                                                                  const std::vector<double>& weights,
                                                                                  double step_factor,
                                                                                  double prec = 1e-8,
                                                                                  int n_itrs = 1000,
                                                                                  bool mean_adjust = false) {
    int n_vertices = y.size();

    if (n_vertices != (int)graph.size())
        Rf_error("The lengths of graph and y must be the same.");

    if (n_vertices != (int)weights.size())
        Rf_error("The lengths of weights and y must be the same.");

    if (prec < 0)
        Rf_error("prec has to be a positive number.");

    if (n_itrs < 1)
        Rf_error("n_itrs has to be greater than or equal to 1.");

    double step_size = calculate_smallest_difference(y);
    step_size *= step_factor;

    // computing the mean of y
    double Ey = 0;
    if (mean_adjust) {
        for (int vertex = 0; vertex < n_vertices; ++vertex)
            Ey += y[vertex];
        Ey /= n_vertices;
    }

    auto y_traj = std::make_unique<std::vector<std::vector<double>>>();
    auto y_t = y; // Initialize y_t to y

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    double average_neighbor_value = 0;
    double Ey_t = 0;
    double Delta;

    // Identifying constant vertices
    std::unique_ptr<std::vector<int>> loc_const_vertices_vect = loc_const_vertices(graph, y, prec);

    std::vector<int> sorted_vertices(n_vertices);

    if ((*loc_const_vertices_vect).empty()) {
        y_traj->push_back(y_t);
        return y_traj;
    } else {
        // Sorting indices so that constant vertices are at the beginning of the sorted vector
        auto const_vertices = *loc_const_vertices_vect;

        std::vector<int> all_indices(n_vertices);
        std::iota(all_indices.begin(), all_indices.end(), 0);
        std::vector<int> non_const_indices;

        std::set_difference(all_indices.begin(), all_indices.end(),
                            const_vertices.begin(), const_vertices.end(),
                            std::back_inserter(non_const_indices));

        // Concatenate the vectors, maintaining the order where constant vertices come first
        std::copy(const_vertices.begin(), const_vertices.end(), sorted_vertices.begin());
        std::copy(non_const_indices.begin(), non_const_indices.end(), sorted_vertices.begin() + const_vertices.size());
    }

    int itr = 0;
    while (itr < n_itrs && has_constant_vertex(graph, y_t, prec, sorted_vertices)) {
        y_traj->push_back(y_t);
        for (int vertex = 0; vertex < n_vertices; ++vertex) {
            average_neighbor_value = 0;
            for (int neighbor : graph[vertex])
                average_neighbor_value += y_t[neighbor];
            average_neighbor_value /= n_neighbors[vertex];

            y_t[vertex] += weights[vertex] * step_size * (average_neighbor_value - y_t[vertex]);
        }

        if (mean_adjust) {
            // Adjusting y_t so that its mean = mean(y)
            Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += y_t[vertex];
            Ey_t /= n_vertices;

            Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                y_t[vertex] += Delta;
        }

        itr++;
    }

    y_traj->push_back(y_t);

    return y_traj;
}


/**
 * An R interface to the make_response_locally_non_const function.
 *
 * This function takes an R list representing a graph adjacency list, a numeric vector
 * of function values at each vertex, a numeric vector of vertex weights, and various
 * parameters for the smoothing process. It calls the C++ weighted_heat_smoother_trajectory
 * function to smooth the function values using the weighted heat diffusion equation and
 * returns a list containing the trajectory of the smoothed function values over time,
 * along with the local minima and local maxima at each iteration.
 *
 * @param r_adj_list An R list representing the graph adjacency list. Each element of the list
 *               should be an integer vector specifying the neighboring vertex indices for
 *               a given vertex.
 * @param Ry A numeric vector representing the initial function values at each vertex of the graph.
 * @param Rweights A numeric vector representing the weights of each vertex in the graph.
 *                 The length of Rweights must be equal to the number of vertices.
 * @param Rstep_factor A numeric scalar specifying the step factor used to calculate the step size
 *                     in the smoothing process. The step size is the product of the step factor
 *                     and the smallest difference between the unique values of the function.
 * @param Rprec A numeric scalar specifying the precision threshold used to determine local
 *              constancy. If the absolute difference between the value at a vertex and any
 *              of its neighbors is greater than this threshold, the vertex is not considered
 *              locally constant.
 * @param Rn_itrs An integer scalar specifying the maximum number of iterations for the smoothing
 *                process. This parameter is used only when the stopping criterion is set to 2.
 * @param Rmean_adjust A boolean scalar specifying whether to adjust the mean of the smoothed function
 *                     values at each iteration to match the mean of the initial function values.
 *                     If Rmean_adjust is true, the mean adjustment is performed; otherwise, it is skipped.
 *
 * @return An R matrix containing the trajectory of the smoothed function values over time, where each
 *         row represents the function values at a specific iteration.
 */
SEXP S_make_response_locally_non_const(SEXP r_adj_list,
                                       SEXP Ry,
                                       SEXP Rweights,
                                       SEXP Rstep_factor,
                                       SEXP Rprec,
                                       SEXP Rn_itrs,
                                       SEXP Rmean_adjust) {

    // ---- Basic validation (cheap checks; no allocations here) ----
    if (!Rf_isVectorList(r_adj_list)) {
        Rf_error("`adj.list` must be a list of integer vectors.");
    }
    if (!Rf_isNumeric(Ry))       Rf_error("`y` must be a numeric vector.");
    if (!Rf_isNumeric(Rweights)) Rf_error("`weights` must be a numeric vector.");

    // Scalars: use Rf_as* (no extra PROTECT needed)
    const double step_factor = Rf_asReal(Rstep_factor);
    const double prec        = Rf_asReal(Rprec);
    const int    n_itrs      = Rf_asInteger(Rn_itrs);
    const int    mean_adjust = Rf_asLogical(Rmean_adjust);

    if (!(prec > 0))        Rf_error("`prec` must be positive.");
    if (!(n_itrs >= 1))     Rf_error("`n_itrs` must be >= 1.");
    if (mean_adjust == NA_LOGICAL) Rf_error("`mean_adjust` must be TRUE or FALSE.");

    // ---- Convert adjacency (pure C++ conversion) ----
    std::vector<std::vector<int>> graph = convert_adj_list_from_R(r_adj_list);

    const size_t n_vertices = (size_t) LENGTH(Ry);
    std::vector<double> y(REAL(Ry), REAL(Ry) + n_vertices);
    const size_t n_weights = (size_t) LENGTH(Rweights);
    std::vector<double> weights(REAL(Rweights), REAL(Rweights) + n_weights);

    // ---- Cross-check sizes ----
    if (n_vertices != graph.size()) {
        Rf_error("Length of `y` (%zu) must equal number of vertices in the graph (%zu).",
                 n_vertices, graph.size());
    }
    if (n_weights != n_vertices) {
        Rf_error("Length of `weights` (%zu) must equal length of `y` (%zu).",
                 n_weights, n_vertices);
    }

    // ---- Compute trajectory ----
    std::unique_ptr<std::vector<std::vector<double>>> y_traj;
    try {
        y_traj = make_response_locally_non_const(
            graph, y, weights, step_factor, prec, n_itrs, (mean_adjust == TRUE)
            );
    } catch (const std::exception& e) {
        Rf_error("make_response_locally_non_const() failed: %s", e.what());
    } catch (...) {
        Rf_error("make_response_locally_non_const() failed with an unknown error.");
    }
    if (!y_traj) {
        Rf_error("Internal error: computation returned null result.");
    }

    // ---- Assemble result matrix: n_vertices x n_iterations (column-major) ----
    const size_t n_iterations = y_traj->size();
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, n_vertices, n_iterations));
    double* result_ptr = REAL(result);
    for (size_t i = 0; i < n_iterations; ++i) {
        const std::vector<double>& col = (*y_traj)[(size_t)i];
        std::copy(col.begin(), col.end(), result_ptr + i * n_vertices);
    }

    UNPROTECT(1); // result
    return result;
}


// ------------------------------------------------------------------------------------------
//
// Local extrema of functions over verties of a graph
//
// ------------------------------------------------------------------------------------------

/**
 * Finds the local extrema (minima and maxima) of a function defined over the vertices of a graph.
 *
 * This function takes a graph represented as an adjacency list and a vector of function values
 * at each vertex. It identifies the vertices at which the function attains local minima and maxima.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph. The size of `y` should
 *          be equal to the number of vertices in the graph.
 * @return A unique pointer to a pair of vectors, where the first vector contains the indices of vertices
 *         at which the function attains local minima, and the second vector contains the indices of vertices
 *         at which the function attains local maxima.
 */
std::unique_ptr<std::pair<std::vector<int>, std::vector<int>>> find_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y) {
    int n_vertices = graph.size();
    std::vector<int> lmin; // indices of local minima
    std::vector<int> lmax; // indices of local maxima

    for (int vertex = 0; vertex < n_vertices; vertex++) {
        // Checking if y attains a local maximum at vertex
        bool is_lmax = true;
        for (int neighbor : graph[vertex]) {
            if (!(y[vertex] > y[neighbor])) {
                is_lmax = false;
                break;
            }
        }
        if (is_lmax) { // Vertex is a local maximum of y
            lmax.push_back(vertex);
        } else {
            // Checking if y attains a local minimum at vertex
            bool is_lmin = true;
            for (int neighbor : graph[vertex]) {
                if (!(y[vertex] < y[neighbor])) {
                    is_lmin = false;
                    break;
                }
            }
            if (is_lmin) { // vertex is a local minimum of y
                lmin.push_back(vertex);
            }
        }
    }

    auto res = std::make_unique<std::pair<std::vector<int>, std::vector<int>>>(std::make_pair(lmin, lmax));

    return res;
}

/**
 * Finds the connected components of local extrema (minima and maxima) of a function defined over the vertices of a graph.
 *
 * This function takes a graph represented as an adjacency list and a vector of function values at each vertex.
 * It identifies the vertices at which the function attains local minima and maxima, and then finds the connected
 * components of the local minimum vertices and the local maximum vertices separately.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector represents a vertex,
 *              and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph. The size of `y` should be equal to
 *          the number of vertices in the graph.
 * @return A unique pointer to a `conn_comps_of_loc_extr_t` struct containing two unordered maps:
 *         - `lmin_cc_map`: An unordered map where the keys are the connected component indices, and the values are
 *                         vectors of indices of vertices at which the function attains local minima within that component.
 *         - `lmax_cc_map`: An unordered map where the keys are the connected component indices, and the values are
 *                         vectors of indices of vertices at which the function attains local maxima within that component.
 *
 * @note The connected components are identified separately for local minima and local maxima.
 */
std::unique_ptr<conn_comps_of_loc_extr_t> find_conn_comps_of_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                           const std::vector<double>& y) {

    auto lextr_res = find_local_extrema(graph, y);
    auto lmin = lextr_res->first;   // indices of local minima
    auto lmax = lextr_res->second;  // indices of local maxima

    std::unique_ptr<std::unordered_map<int, int>> lmin_cc_map_uptr = count_subgraph_components(graph, lmin);
    std::unique_ptr<std::unordered_map<int, int>> lmax_cc_map_uptr = count_subgraph_components(graph, lmax);

    auto res = std::make_unique<conn_comps_of_loc_extr_t>();

    // Creating a map assigning to each connected component a vector of local minima of that connected component
    for (const auto& pair : *lmin_cc_map_uptr)
        res->lmin_cc_map[pair.second].push_back(pair.first);

    // Doing the same for local maxima
    for (const auto& pair : *lmax_cc_map_uptr)
        res->lmax_cc_map[pair.second].push_back(pair.first);

    return res;
}

/**
 * Checks if all local extrema (minima and maxima) of a graph function are isolated.
 *
 * Local extrema of a graph function are considered isolated if each connected component
 * of local extrema consists of a single vertex. This function takes the result of a call
 * to the `find_conn_comps_of_local_extrema` function and checks if all local minima and
 * maxima are isolated.
 *
 * @param res A unique pointer to a `conn_comps_of_loc_extr_t` struct, which is the result
 *            of a call to the `find_conn_comps_of_local_extrema` function. The struct
 *            contains the connected component information for local minima and maxima.
 *
 * @return `true` if all local extrema (minima and maxima) are isolated, meaning that each
 *         connected component of local extrema consists of a single vertex; `false`
 *         otherwise.
 *
 * @note This function assumes that the input `res` struct is properly populated with the
 *       connected component information for local minima and maxima.
 */
bool are_local_extrema_isolated(std::unique_ptr<conn_comps_of_loc_extr_t>& res) {
    // Checking if all local minima are isolated
    for (const auto& [key, vect] : res->lmin_cc_map) {
        if (vect.size() > 1) {
            return false;
        }
    }

    // Checking if all local maxima are isolated
    for (const auto& [key, vect] : res->lmax_cc_map) {
        if (vect.size() > 1) {
            return false;
        }
    }

    return true;
}

/**
 * Adds a small amount of noise to a random vertex from each connected component
 * of local extrema of a graph function that contains more than one vertex.
 *
 * When finding gradient flow of y trajectories, it is optimal if the local
 * extrema of y are isolated. This function makes them isolated by slightly
 * perturbing the function values at selected vertices.
 *
 * @param res A unique pointer to a `conn_comps_of_loc_extr_t` struct, which is the result
 *            of a call to the `find_conn_comps_of_local_extrema` function. The struct
 *            contains the connected component information for local minima and maxima.
 *            - `res->lmin_cc_map`: An unordered map where the keys are the connected component indices,
 *                                  and the values are vectors of indices of vertices at which the function
 *                                  attains local minima within that component.
 *            - `res->lmax_cc_map`: An unordered map where the keys are the connected component indices,
 *                                  and the values are vectors of indices of vertices at which the function
 *                                  attains local maxima within that component.
 * @param y A vector representing the function values at each vertex of the graph. The size
 *          of `y` should be equal to the number of vertices in the graph.
 * @param scaling_factor A scaling factor used to define the amount of noise added to y at
 *                       selected vertices. The noise is `eps = scaling_factor * d`, where
 *                       `d` is the smallest difference between sorted y values. Default
 *                       value is 0.2.
 *
 * @return A unique pointer to a new vector containing the modified function values, where
 *         a small amount of noise is added to a random vertex from each connected component
 *         of local extrema that contains more than one vertex.
 *
 * @note This function modifies the input function values `y` by adding noise to selected
 *       vertices. The amount of noise added is determined by the `scaling_factor` parameter.
 *       The noise is added to a random vertex from each connected component of local extrema
 *       that contains more than one vertex, in order to make the local extrema isolated.
 *
 * @note The function assumes that the input `res` struct is properly populated with the
 *       connected component information for local minima and maxima, as obtained from the
 *       `find_conn_comps_of_local_extrema` function.
 */
std::unique_ptr<std::vector<double>> jitter_non_isolated_loc_extrema(std::unique_ptr<conn_comps_of_loc_extr_t>& res,
                                                                     const std::vector<double>& y,
                                                                     double scaling_factor = 0.2) {
    double d = calculate_smallest_difference(y);
    double eps = scaling_factor * d;

    std::vector<double> jittered_y(y);

    // For each connected component of local minima containing more than one vertex
    GetRNGstate();
    for (const auto& pair : res->lmin_cc_map) {
        if (pair.second.size() > 1) {
            // Pick a random vertex from that component
            int random_index = (int)(unif_rand() * pair.second.size());
            int vertex = pair.second[random_index];

            // Modify the value of jittered_y at that vertex by subtracting eps
            jittered_y[vertex] -= eps;
        }
    }

    // For each connected component of local maxima containing more than one vertex
    for (const auto& pair : res->lmax_cc_map) {
        if (pair.second.size() > 1) {
            // Pick a random vertex from that component
            int random_index = (int)(unif_rand() * pair.second.size());
            int vertex = pair.second[random_index];

            // Modify the value of jittered_y at that vertex by adding eps
            jittered_y[vertex] += eps;
        }
    }
    PutRNGstate();

    return std::make_unique<std::vector<double>>(jittered_y);
}
