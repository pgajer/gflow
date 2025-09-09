#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <set>
#include <unordered_set>
#include <memory>
#include <algorithm>

#include "SEXP_cpp_conversion_utils.h"
#include "cpp_utils.h"

extern "C" {
    SEXP S_create_intersection_graph(SEXP R_graph, SEXP R_p_thld, SEXP R_n_itrs);
}

/**
 * Creates an intersection graph based on the input graph and a threshold.
 *
 * This function constructs a new graph where an edge exists between two nodes
 * if the size of the intersection of their neighborhood sets exceeds a calculated
 * threshold. The function can perform multiple iterations of this graph
 * intersection operation.
 *
 * @param graph The input graph represented as an adjacency list.
 * @param p_thld A relative threshold parameter between 0 and 1.
 * @param n_itrs The number of iterations to perform (default is 1).
 *
 * @return A unique pointer to a vector of sets representing the intersection graph.
 *
 * Operation:
 * 1. For each pair of nodes (i, j), calculate the threshold:
 *    - First iteration: thld = p_thld * min(graph[i].size(), graph[j].size())
 *    - Subsequent iterations: thld = p_thld * min(res[i].size(), res[j].size())
 *      where res is the graph generated in the previous iteration.
 * 2. Create an edge (i, j) if |neighborhoods(i) ∩ neighborhoods(j)| > thld
 * 3. Repeat for n_itrs iterations, using the result of the previous iteration
 *    as input for the next.
 *
 * In each iteration:
 * - First iteration uses the input graph to calculate intersections.
 * - Subsequent iterations use the result from the previous iteration.
 * - Neighborhood sizes for threshold calculation are updated based on the previous result.
 *
 * Time Complexity: O(n_itrs * V^2 * (V + E)), where V is the number of vertices
 * and E is the number of edges in the input graph.
 * Space Complexity: O(V^2) for the result, O(V) for temporary sets and intersections.
 */
std::unique_ptr<std::vector<std::set<int>>> create_intersection_graph(
    const std::vector<std::vector<int>>& graph,
    double p_thld,
    int n_itrs = 1)
{
    auto result = std::make_unique<std::vector<std::set<int>>>(graph.size());
    std::set<int> set_i, set_j;

    std::vector<int> intersection;
    double thld;
    for (int itr = 0; itr < n_itrs; itr++) {

        for (size_t i = 0; i < graph.size(); ++i) {
            set_i.clear();

            if (itr == 0) {
                for (int neighbor : graph[i])
                    set_i.insert(neighbor);
            } else {
                for (int neighbor : (*result)[i])
                    set_i.insert(neighbor);
            }

            for (size_t j = i + 1; j < graph.size(); ++j) {
                set_j.clear();

                if (itr == 0) {
                    for (int neighbor : graph[j])
                        set_j.insert(neighbor);
                } else {
                    for (int neighbor : (*result)[j])
                    set_j.insert(neighbor);
                }

                // Calculate the threshold for this pair
                if (itr == 0) {
                    thld = p_thld * std::min(graph[i].size(), graph[j].size());
                } else {
                    thld = p_thld * std::min((*result)[i].size(), (*result)[j].size());
                }

                intersection.clear();
                std::set_intersection(set_i.begin(), set_i.end(),
                                      set_j.begin(), set_j.end(),
                                      std::back_inserter(intersection));

                if ((int)intersection.size() > thld) {
                    (*result)[i].insert(j);
                    (*result)[j].insert(i);
                }
            }
        }
    }

    return result;
}

/**
 * @brief Creates an intersection graph from an input graph based on neighborhood intersections.
 *
 * This function constructs a new graph where an edge exists between two nodes if their
 * neighborhoods in the input graph have a significant intersection. The significance
 * is determined by a threshold that's proportional to the size of the smaller neighborhood.
 *
 * @param R_graph SEXP (LIST) The input graph represented as a list of integer vectors.
 *                Each vector contains the indices of the neighbors for the corresponding node.
 *                The vectors should be sorted in ascending order.
 * @param R_p_thld SEXP (REAL) A threshold parameter between 0 and 1. The actual threshold
 *                 for each pair of nodes (i, j) is calculated as:
 *                 thld = p_thld * min(length(graph[i]), length(graph[j]))
 *
 * @return SEXP (LIST) A new graph (list of integer vectors) where an edge (i, j) exists
 *         if |graph[i] ∩ graph[j]| >= thld. The returned indices are 1-based (R convention).
 *
 * @note The time complexity of this function is O(n^2 * k), where n is the number of nodes
 *       and k is the average neighborhood size.
 *
 * @note This function assumes that the input vectors are sorted in ascending order.
 *       If they are not, unexpected results may occur.
 *
 * @warning The function does not extensively check if p_thld is between 0 and 1.
 *          Using values outside this range may lead to unexpected results.
 *
 * @example
 * In R:
 * result <- .Call("S_create_intersection_graph",
 *                 list(c(1L,2L), c(0L,2L,3L), c(0L,1L), c(1L)),
 *                 0.5)
 *
 * @see For background on intersection graphs, refer to:
 *      "Visualizing structure and transitions in high-dimensional biological data"
 *      by Moon, et al. (2019)
 */
SEXP S_create_intersection_graph(SEXP R_graph, SEXP R_p_thld, SEXP R_n_itrs) {

    std::vector<std::vector<int>> graph = Rgraph_to_vector(R_graph);
    int n = Rf_length(R_graph);
    double p_thld = REAL(R_p_thld)[0];
    int n_itrs = INTEGER(R_n_itrs)[0];

    auto igraph = std::move(*create_intersection_graph(graph, p_thld, n_itrs));

    // Create result list
    SEXP R_result = PROTECT(allocVector(VECSXP, n));

    // Populate R_result using igraph
    for (int i = 0; i < n; i++) {
        int size = igraph[i].size();
        SEXP R_neighbors = PROTECT(allocVector(INTSXP, size));
        int* neighbors_ptr = INTEGER(R_neighbors);

        int j = 0;
        for (int neighbor : igraph[i]) {
            neighbors_ptr[j++] = neighbor + 1;  // Convert to 1-based indexing for R
        }

        SET_VECTOR_ELT(R_result, i, R_neighbors);
        UNPROTECT(1);  // Unprotect R_neighbors
    }

    UNPROTECT(1);  // Unprotect R_result
    return R_result;
}



