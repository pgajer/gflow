#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm> // Required for std::set_intersection and std::sort
#include <memory>
#include <utility>
#include <limits>

#include "msr2.h"
#include "kNN.h"         // for struct IkNN_vertex_tt and kNN_search
#include "cpp_utils.h"

extern "C" {
    SEXP S_mst_kNN_graph(SEXP RX, SEXP Rk);
}

/*!
 * @brief Constructs a k-Nearest Neighbors (kNN) graph based on a Minimum Spanning Tree (MST)
 *
 * This function implements a kNN graph construction algorithm using the following steps:
 * 1. Compute the k-nearest neighbors for all points using Euclidean distance.
 * 2. Construct a Minimum Spanning Tree (MST) using the k-nearest neighbors information.
 * 3. For each point, perform a breadth-first search on the MST to find its k-nearest neighbors
 *    based on the MST distance (sum of edge weights along the path).
 *
 * The algorithm provides an approximation of the true k-nearest neighbors, which can be
 * computed more efficiently than an exact kNN search, especially for high-dimensional data.
 *
 * @param RX An R matrix of points, where each column represents a point in the dataset.
 * @param Rk An R integer specifying the number of nearest neighbors to find.
 *
 * @return An R list containing two components:
 *         - adj_list: A list of integer vectors, where each vector contains the indices
 *                     of the k-nearest neighbors for the corresponding point.
 *         - dist_list: A list of numeric vectors, where each vector contains the MST-based
 *                      distances to the k-nearest neighbors for the corresponding point.
 *
 * @note This function uses R's C API and assumes that the input matrix is in column-major order.
 * @warning The function modifies the global state of the ANN library. Make sure to call
 *          annClose() when you're done using ANN functions in your R session.
 */
SEXP S_mst_kNN_graph(SEXP RX, SEXP Rk) {

    PROTECT(RX = coerceVector(RX, REALSXP));
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    UNPROTECT(1);
    int nr_X = dimX[0];

    PROTECT(Rk = coerceVector(Rk, INTSXP));
    int k = INTEGER(Rk)[0];
    UNPROTECT(1);

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));
    UNPROTECT(1);

    // Transposing indices to nn_i
    // Transposing distances to nn_d
    int nn_cols = std::min(k, nr_X - 1);

    std::vector<std::vector<int>> nn_i(nr_X, std::vector<int>(nn_cols));
    std::vector<std::vector<double>> nn_d(nr_X, std::vector<double>(nn_cols));

    for (int i = 0; i < nr_X; ++i) {
        for (int j = 0; j < nn_cols; ++j) {
            nn_i[i][j] = indices[j * nr_X + i];
            nn_d[i][j] = distances[j * nr_X + i];
        }
    }

    std::vector<int> nn_i_flat(nr_X * nn_cols);
    std::vector<double> nn_d_flat(nr_X * nn_cols);

    for (int i = 0; i < nr_X; ++i) {
        for (int j = 0; j < nn_cols; ++j) {
            nn_i_flat[i * nn_cols + j] = nn_i[i][j];
            nn_d_flat[i * nn_cols + j] = nn_d[i][j];
        }
    }


    // Computing ldist - a value greater than all distances in nn_d (used as infinity).
    double ldist = *std::max_element(distances, distances + nr_X * k) + 1.0;

    // Allocating memory for edges and edge_lengs
    int nr_X_minus_one = nr_X - 1;
    std::vector<int> edges_vec(nr_X_minus_one * 2);
    std::vector<double> edge_lengs_vec(nr_X_minus_one);

    // Constructing the MST of X
    int init = 0;
    C_mstree(&init,
             nn_i_flat.data(),
             nn_d_flat.data(),
             &ldist,
             &nr_X,
             edges_vec.data(),
             edge_lengs_vec.data());

    // Creating an adjacency 'list' ,mst_adj_list, of the MST
    std::vector<std::vector<int>> mst_adj_list(nr_X);
    for (int i = 0; i < nr_X_minus_one; ++i) {
        int u = edges_vec[2 * i];
        int v = edges_vec[2 * i + 1];
        mst_adj_list[u].push_back(v);
        mst_adj_list[v].push_back(u);
    }

    // Creating an edge lenghts 'list', mst_dist_list, of the MST
    std::vector<std::vector<double>> mst_dist_list(nr_X);
    for (int i = 0; i < nr_X_minus_one; ++i) {
        int u = edges_vec[2 * i];
        int v = edges_vec[2 * i + 1];
        double len = edge_lengs_vec[i];
        mst_dist_list[u].push_back(len);
        mst_dist_list[v].push_back(len);
    }

    /**
     * @brief Finds the k-nearest neighbors of a given vertex in the Minimum Spanning Tree (MST)
     *
     * This lambda function implements a modified Breadth-First Search (BFS) algorithm to find
     * the k-nearest neighbors of a specified vertex in the MST. The distance metric used is
     * the sum of edge weights along the path in the MST, rather than Euclidean distance.
     *
     * The algorithm works as follows:
     * 1. Initialize a vector to store candidate neighbors with their distances.
     * 2. Initialize a distance array to keep track of the shortest known distance to each vertex.
     * 3. Initialize a hop_count array to track the number of edges traversed to reach each vertex.
     * 4. Start with the given vertex, setting its distance to 0 and hop count to 0.
     * 5. Use a queue for BFS traversal of the MST.
     * 6. While the queue is not empty:
     *    a. Dequeue a vertex.
     *    b. If the hop count for this vertex is >= k, skip it.
     *    c. If this is not the starting vertex, add it to the candidates list.
     *    d. For each adjacent vertex in the MST:
     *       - Calculate the new distance and hop count to this adjacent vertex.
     *       - If the new hop count is less than the known hop count:
     *         * Update the distance and hop count.
     *         * Enqueue the adjacent vertex.
     * 7. Sort the candidates by distance.
     * 8. Return the first k candidates (or all if there are fewer than k).
     *
     * @param vertex The index of the vertex for which to find the k-nearest neighbors
     *
     * @return A vector of pairs, where each pair contains:
     *         - The distance to the neighbor along the MST
     *         - The index of the neighbor
     *         The vector is sorted by distance in ascending order and contains at most k elements.
     *
     * @note This function assumes the existence of the following variables in its scope:
     *       - nr_X: The total number of vertices in the MST
     *       - k: The number of nearest neighbors to find
     *       - mst_adj_list: A vector of vectors representing the adjacency list of the MST
     *       - mst_dist_list: A vector of vectors representing the edge weights in the MST
     *
     * @warning This function modifies no external state and is thread-safe, but it does
     *          depend on the correct initialization of the MST structure (mst_adj_list and mst_dist_list).
     *
     * Time Complexity: O(V log V), where V is the number of vertices in the MST.
     * Space Complexity: O(V), for the candidates vector, distance array, and hop_count array.
     */
    auto find_mst_kNNs = [&](int vertex) {
        std::vector<std::pair<double, int>> candidates;
        std::vector<double> distances(nr_X, std::numeric_limits<double>::max());
        std::vector<int> hop_count(nr_X, std::numeric_limits<int>::max());
        std::queue<int> q;

        q.push(vertex);
        distances[vertex] = 0;
        hop_count[vertex] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            if (hop_count[u] >= k) continue;

            if (u != vertex) {
                candidates.push_back({distances[u], u});
            }

            for (size_t i = 0; i < mst_adj_list[u].size(); ++i) {
                int v = mst_adj_list[u][i];
                double new_dist = distances[u] + mst_dist_list[u][i];
                int new_hop = hop_count[u] + 1;

                if (new_hop < hop_count[v]) {
                    distances[v] = new_dist;
                    hop_count[v] = new_hop;
                    q.push(v);
                }
            }
        }

        std::sort(candidates.begin(), candidates.end());
        candidates.resize(std::min(k, static_cast<int>(candidates.size())));

        return candidates;
    };


    // The main loop computing MST kNN's of each point of X
    std::vector<std::vector<int>> mst_knn_indices(nr_X);
    std::vector<std::vector<double>> mst_knn_distances(nr_X);

    for (int i = 0; i < nr_X; ++i) {
        auto knn = find_mst_kNNs(i);
        for (const auto& [idx, dist] : knn) {
            mst_knn_indices[i].push_back(idx);
            mst_knn_distances[i].push_back(dist);
        }
    }

    // Preparing adjacency list of the constructed MST kNN graph
    SEXP adj_list = PROTECT(allocVector(VECSXP, nr_X));
    for (int i = 0; i < nr_X; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, mst_knn_indices[i].size()));
        int* A = INTEGER(RA);
        for (const auto& neighbor : mst_knn_indices[i])
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1);
    }
    UNPROTECT(1);

    // Preparing distance list of the constructed MST kNN graph
    SEXP dist_list = PROTECT(allocVector(VECSXP, nr_X));
    for (int i = 0; i < nr_X; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, mst_knn_distances[i].size()));
        double* D = REAL(RD);
        for (const auto& dist : mst_knn_distances[i])
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
        UNPROTECT(1);
    }
    UNPROTECT(1);


    // Creating a return list with two components
    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, dist_list);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("dist_list"));
    setAttrib(res, R_NamesSymbol, names);
    UNPROTECT(2);

    return res;
}
