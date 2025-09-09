#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Undefine conflicting macros after including R headers
#undef length

#include <execution>     // For parallel algorithms
#include <mutex>         // For thread safety
#include <numeric>       // For std::iota
#include <algorithm>     // For std::sort, std::set_intersection

#include "msr2.h"
#include "kNN.h"         // for struct iknn_vertex_tt and kNN_search
#include "iknn_graphs.hpp"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

/**
 * @brief Creates an Inverse Weight Distance k-Nearest Neighbors graph using parallel processing
 *
 * @param RX SEXP object containing input data matrix
 * @param Rk SEXP object containing k value for k-NN
 * @return std::unique_ptr to vector of vectors containing iknn_vertex_t objects
 */
std::unique_ptr<std::vector<std::vector<iknn_vertex_t>>> create_iknn_graph_parallel(SEXP RX, SEXP Rk) {
    // Protect R objects and get dimensions
    PROTECT(RX = coerceVector(RX, REALSXP));
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int n_points = dimX[0];
    UNPROTECT(1);

    // Get k value
    PROTECT(Rk = coerceVector(Rk, INTSXP));
    int k = INTEGER(Rk)[0];
    UNPROTECT(1);

    // Get kNN results
    SEXP knn_res = PROTECT(S_kNN(RX, Rk));
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    // Create result container
    auto res = std::make_unique<std::vector<std::vector<iknn_vertex_t>>>(n_points);

    // Create vector of indices for parallel processing
    std::vector<int> point_indices(n_points - 1);
    std::iota(point_indices.begin(), point_indices.end(), 0);

    // Mutex vector for thread-safe writing to result vectors
    std::vector<std::mutex> vertex_mutexes(n_points);

    // Parallel processing of points
    std::for_each(std::execution::par_unseq,
                  point_indices.begin(),
                  point_indices.end(),
                  [&](int pt_i) {
        // Thread-local storage
        std::vector<int> nn_i(k);
        std::vector<int> nn_j(k);
        std::vector<int> sorted_nn_i(k);
        std::vector<int> sorted_nn_j(k);
        std::vector<int> intersection;
        intersection.reserve(k); // Pre-allocate for efficiency

        // Get and sort kNN indices for point i
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + n_points * j];
            sorted_nn_i[j] = nn_i[j];
        }
        std::sort(sorted_nn_i.begin(), sorted_nn_i.end());

        // Process connections with later points
        for (int pt_j = pt_i + 1; pt_j < n_points; pt_j++) {
            // Get and sort kNN indices for point j
            for (int j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + n_points * j];
                sorted_nn_j[j] = nn_j[j];
            }
            std::sort(sorted_nn_j.begin(), sorted_nn_j.end());

            // Find intersection
            intersection.clear();
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(),
                                sorted_nn_j.begin(), sorted_nn_j.end(),
                                std::back_inserter(intersection));

            int common_count = intersection.size();
            if (common_count > 0) {
                // Find minimum distance through common neighbors
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    int idx_i = static_cast<int>(std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin());
                    int idx_j = static_cast<int>(std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin());
                    double dist_i_k = distances[pt_i + n_points * idx_i];
                    double dist_j_k = distances[pt_j + n_points * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                // Create new vertex
                iknn_vertex_t new_vertex_i{pt_j, common_count, min_dist};
                iknn_vertex_t new_vertex_j{pt_i, common_count, min_dist};

                // Safely add edges using mutexes
                {
                    std::lock_guard<std::mutex> lock_i(vertex_mutexes[pt_i]);
                    (*res)[pt_i].push_back(new_vertex_i);
                }
                {
                    std::lock_guard<std::mutex> lock_j(vertex_mutexes[pt_j]);
                    (*res)[pt_j].push_back(new_vertex_j);
                }
            }
        }
    });

    UNPROTECT(1); // knn_res
    return res;
}
