#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>

// Undefine conflicting macros after including R headers
#undef length

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

#include <omp.h>

#include <Eigen/Sparse>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

extern "C" {
    SEXP S_diffusion_graph(SEXP R_neighbor_index,
                           SEXP R_neighbor_dist,
                           SEXP R_alpha,
                           SEXP R_t,
                           SEXP R_affinity_threshold,
                           SEXP R_threshold_strategy,
                           SEXP R_percentile);
}

struct diffusion_graph_t {
    Eigen::SparseMatrix<double> D_alpha;
    std::vector<double> affinities;
};

/**
 * @brief Computes the diffusion graph based on the PHATE algorithm.
 *
 * This function constructs a diffusion graph from a k-nearest neighbors graph
 * representation of the data. It applies adaptive bandwidth, thresholding,
 * and a diffusion process to capture both local and global structure of the data.
 *
 * @param neighbor_index Vector of vectors containing indices of neighbors for each point.
 * @param neighbor_dist Vector of vectors containing distances to neighbors for each point.
 * @param alpha Alpha decay factor for the diffusion process (0 < alpha <= 1).
 *        - Lower values emphasize global structure more.
 *        - Higher values preserve more local structure.
 *        - Typical values range from 0.5 to 1.
 * @param n_diff_steps Number of time steps for the diffusion process.
 *        - Lower values capture more local structure.
 *        - Higher values capture more global structure.
 *        - Typical values range from 1 to 100.
 * @param affinity_threshold Threshold for affinity values (used when threshold_strategy = 0).
 *        - Affinities below this value are set to zero.
 *        - Higher values result in a sparser graph.
 * @param threshold_strategy Strategy for thresholding (0: fixed, 1: percentile-based).
 *        - 0: Use the fixed affinity_threshold value.
 *        - 1: Use a percentile-based threshold.
 * @param percentile Percentile to use when threshold_strategy = 1 (0 <= percentile <= 100).
 *        - Determines the cutoff point for affinities.
 *        - Lower values result in a denser graph, higher values in a sparser graph.
 *        - E.g., percentile = 5 means the lowest 5% of affinities are set to zero.
 *
 * @return std::unique_ptr<diffusion_graph_t> Pointer to a struct containing:
 *         - D_alpha: The resulting diffusion graph as a sparse matrix.
 *         - affinities: Vector of all computed affinities before thresholding.
 */
std::unique_ptr<diffusion_graph_t> diffusion_graph(
    const std::vector<std::vector<int>>& neighbor_index,
    const std::vector<std::vector<double>>& neighbor_dist,
    double alpha,
    int n_diff_steps,
    double affinity_threshold = 0.0,
    int threshold_strategy = 0,
    double percentile = 5.0)
{
    int n_vertices = neighbor_index.size();

    // Compute adaptive bandwidth
    std::vector<double> sigma(n_vertices);
    for (int i = 0; i < n_vertices; ++i) {
        sigma[i] = *std::max_element(neighbor_dist[i].begin(), neighbor_dist[i].end());
    }

    // Collect all affinities
    std::vector<double> all_affinities;
    for (int i = 0; i < n_vertices; ++i) {
        for (size_t j = 0; j < neighbor_index[i].size(); ++j) {
            double dist = neighbor_dist[i][j];
            double affinity = std::exp(-std::pow(dist, 2) / (2 * std::pow(sigma[i], 2)));
            all_affinities.push_back(affinity);
        }
    }

    // Apply threshold strategy
    if (threshold_strategy == 1) {
        auto percentile_threshold = [](std::vector<double>& affinities, double percentile) {
            std::sort(affinities.begin(), affinities.end());
            int index = static_cast<int>(percentile / 100.0 * affinities.size());
            return affinities[index];
        };
        affinity_threshold = percentile_threshold(all_affinities, percentile);
    }

    // Construct sparse affinity matrix
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> tripletList;
    for (int i = 0; i < n_vertices; ++i) {
        for (size_t j = 0; j < neighbor_index[i].size(); ++j) {
            int idx = neighbor_index[i][j];
            double dist = neighbor_dist[i][j];
            double affinity = std::exp(-std::pow(dist, 2) / (2 * std::pow(sigma[i], 2)));
            if (affinity > affinity_threshold) {
                tripletList.emplace_back(i, idx, affinity);
            }
        }
    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    // Row-normalize to create Markov transition matrix
    Eigen::VectorXd row_sums = A * Eigen::VectorXd::Ones(n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            it.valueRef() /= row_sums(it.row());
        }
    }

    // Compute diffusion operator
    Eigen::SparseMatrix<double> D = A;
    Eigen::SparseMatrix<double> P_t = A;
    for (int i = 1; i < n_diff_steps; ++i) {
        P_t = P_t * A;
        D += P_t;
    }

    // Apply alpha decay
    for (int k = 0; k < D.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D, k); it; ++it) {
            it.valueRef() = std::pow(it.value(), alpha);
        }
    }

    // Create and return the result
    auto result = std::make_unique<diffusion_graph_t>();
    result->D_alpha = std::move(D);
    result->affinities = std::move(all_affinities);

    return result;
}

/**
 * @brief R interface for the diffusion graph computation based on the PHATE algorithm.
 *
 * This function serves as an interface between R and the C++ implementation of the
 * diffusion graph computation. It takes R objects as input, converts them to C++
 * data structures, calls the diffusion_graph function, and converts the result back
 * to R objects.
 *
 * @param R_neighbor_index SEXP (LIST) List of integer vectors containing neighbor indices.
 * @param R_neighbor_dist SEXP (LIST) List of numeric vectors containing neighbor distances.
 * @param R_alpha SEXP (REAL) Alpha decay factor for the diffusion process (0 < alpha <= 1).
 *        Lower values emphasize global structure more, higher values preserve more local structure.
 * @param R_n_diff_steps SEXP (INTEGER) Number of time steps for the diffusion process.
 *        Lower values capture more local structure, higher values capture more global structure.
 * @param R_affinity_threshold SEXP (REAL) Threshold for affinity values (used when threshold_strategy = 0).
 *        Affinities below this value are set to zero.
 * @param R_threshold_strategy SEXP (INTEGER) Strategy for thresholding (0: fixed, 1: percentile-based).
 * @param R_percentile SEXP (REAL) Percentile to use when threshold_strategy = 1 (0 <= percentile <= 100).
 *        Determines the cutoff point for affinities.
 *
 * @return SEXP A list containing two elements:
 *         1. A sparse matrix (dgCMatrix) representing the diffusion graph (D_alpha).
 *         2. A numeric vector of all computed affinities before thresholding.
 *
 * @note This function assumes that the input lists R_neighbor_index and R_neighbor_dist
 *       have the same length and that each corresponding pair of vectors in these lists
 *       have the same length.
 */
SEXP S_diffusion_graph(SEXP R_neighbor_index,
                       SEXP R_neighbor_dist,
                       SEXP R_alpha,
                       SEXP R_n_diff_steps,
                       SEXP R_affinity_threshold,
                       SEXP R_threshold_strategy,
                       SEXP R_percentile) {

    // Convert R lists to C++ vectors
    int n_vertices = Rf_length(R_neighbor_index);
    std::vector<std::vector<int>> neighbor_index(n_vertices);
    std::vector<std::vector<double>> neighbor_dist(n_vertices);

    for (int i = 0; i < n_vertices; ++i) {
        SEXP R_index_vec = VECTOR_ELT(R_neighbor_index, i);
        SEXP R_dist_vec = VECTOR_ELT(R_neighbor_dist, i);
        int m = Rf_length(R_index_vec);

        neighbor_index[i].resize(m);
        neighbor_dist[i].resize(m);

        int *index_ptr = INTEGER(R_index_vec);
        double *dist_ptr = REAL(R_dist_vec);

        for (int j = 0; j < m; ++j) {
            neighbor_index[i][j] = index_ptr[j] - 1;  // Convert to 0-based indexing
            neighbor_dist[i][j] = dist_ptr[j];
        }
    }

    // Extract other parameters
    double alpha = REAL(R_alpha)[0];
    int n_diff_steps = INTEGER(R_n_diff_steps)[0];
    double affinity_threshold = REAL(R_affinity_threshold)[0];
    int threshold_strategy = INTEGER(R_threshold_strategy)[0];
    double percentile = REAL(R_percentile)[0];

    // Call the C++ function
    std::unique_ptr<diffusion_graph_t> result =
        diffusion_graph(neighbor_index,
                        neighbor_dist,
                        alpha,
                        n_diff_steps,
                        affinity_threshold,
                        threshold_strategy,
                        percentile);

    // Create the return list
    SEXP R_result = PROTECT(allocVector(VECSXP, 2));

    // Convert D_alpha to an R sparse matrix (dgCMatrix)
    SEXP R_D_alpha = PROTECT(allocVector(VECSXP, 5));

    // i slot (row indices)
    std::vector<int> row_indices;
    std::vector<int> col_indices;
    std::vector<double> values;
    for (int k = 0; k < result->D_alpha.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(result->D_alpha, k); it; ++it) {
            row_indices.push_back(it.row() + 1);  // Convert to 1-based indexing
            col_indices.push_back(it.col() + 1);  // Convert to 1-based indexing
            values.push_back(it.value());
        }
    }

    SEXP R_row_indices = PROTECT(allocVector(INTSXP, row_indices.size()));
    memcpy(INTEGER(R_row_indices), row_indices.data(), row_indices.size() * sizeof(int));
    SET_VECTOR_ELT(R_D_alpha, 0, R_row_indices);

    // p slot (column pointers)
    SEXP R_col_indices = PROTECT(allocVector(INTSXP, col_indices.size()));
    memcpy(INTEGER(R_col_indices), col_indices.data(), col_indices.size() * sizeof(int));
    SET_VECTOR_ELT(R_D_alpha, 1, R_col_indices);

    // x slot (values)
    SEXP R_values = PROTECT(allocVector(REALSXP, values.size()));
    memcpy(REAL(R_values), values.data(), values.size() * sizeof(double));
    SET_VECTOR_ELT(R_D_alpha, 2, R_values);

    // Dims slot
    SEXP R_dims = PROTECT(allocVector(INTSXP, 2));
    INTEGER(R_dims)[0] = result->D_alpha.rows();
    INTEGER(R_dims)[1] = result->D_alpha.cols();
    SET_VECTOR_ELT(R_D_alpha, 3, R_dims);

    // Dimnames slot
    SET_VECTOR_ELT(R_D_alpha, 4, R_NilValue);

    // Set class
    SEXP R_class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(R_class, 0, mkChar("dgCMatrix"));
    setAttrib(R_D_alpha, R_ClassSymbol, R_class);

    // Set D_alpha as the first element of the return list
    SET_VECTOR_ELT(R_result, 0, R_D_alpha);

    // Convert affinities to an R numeric vector
    SEXP R_affinities = PROTECT(allocVector(REALSXP, result->affinities.size()));
    memcpy(REAL(R_affinities), result->affinities.data(), result->affinities.size() * sizeof(double));

    // Set affinities as the second element of the return list
    SET_VECTOR_ELT(R_result, 1, R_affinities);

    UNPROTECT(8);

    return R_result;
}
