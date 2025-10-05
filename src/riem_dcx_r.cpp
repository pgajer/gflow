#include <Rcpp.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

#include "riem_dcx.hpp"
#include "error_utils.h"

/**
 * @brief Build nerve complex from k-nearest neighbors with Riemannian structure
 *
 * @details
 * This function constructs a simplicial complex from the nerve of a k-nearest
 * neighbor covering of a point cloud, equipped with a Riemannian metric derived
 * from neighborhood overlaps. The construction addresses a fundamental geometric
 * question: how can we capture both the topology (connectivity patterns) and
 * local density (geometric structure) of data through the intersection patterns
 * of neighborhoods?
 *
 * The answer emerges by treating neighborhood intersections as a measure space.
 * Consider a finite point cloud X = {x_1, ..., x_n} ⊂ ℝ^d. For each point x_i,
 * we form its closed k-nearest neighbor set N̂_k(x_i) = N_k(x_i) ∪ {x_i}, where
 * N_k(x_i) contains the k nearest neighbors of x_i. These neighborhoods cover
 * the point cloud, and their intersection pattern encodes topological and
 * geometric information.
 *
 * NERVE CONSTRUCTION (Combinatorial Structure):
 * The nerve complex K has vertices corresponding to points in X, and a simplex
 * {i_0, ..., i_p} belongs to K if and only if the neighborhoods have non-empty
 * intersection:
 *     N̂_k(x_{i_0}) ∩ ... ∩ N̂_k(x_{i_p}) ≠ ∅
 *
 * This purely combinatorial construction determines which simplices exist. The
 * nerve theorem from algebraic topology ensures that, under mild conditions
 * (contractible intersections), the nerve K is homotopy equivalent to the union
 * of the neighborhoods, meaning they share the same topological features.
 *
 * RIEMANNIAN STRUCTURE (Geometric Weights):
 * To enable discrete calculus (gradients, Laplacians, diffusion), we equip the
 * complex with geometry by introducing a measure μ on X. Simplex weights arise
 * from intersection measures, creating a Riemannian structure that supports
 * discrete differential operators respecting local data geometry.
 *
 * For vertex v_i, the mass is:
 *     m_i = μ(N̂_k(x_i))
 *
 * For edge (i,j), the weight is:
 *     m_ij = μ(N̂_k(x_i) ∩ N̂_k(x_j))
 *
 * For higher-dimensional simplices, weights come from higher-order intersections.
 * These weights populate diagonal metric matrices M[p] at each dimension p, which
 * determine the Hodge Laplacians L[p]. The construction ensures positive semi-
 * definiteness by design, as all inner products arise from L²(μ) pairings.
 *
 * MEASURE OPTIONS:
 * Two approaches define the measure μ:
 *
 * 1) Counting Measure (use_counting_measure = true):
 *    Each point receives unit weight: μ({x}) = 1
 *    Measure of a set: μ(A) = |A| (cardinality)
 *    Appropriate when data are approximately uniformly sampled
 *    Treats all points as equally representative
 *
 * 2) Density-Based Measure (use_counting_measure = false):
 *    Points receive weights inversely proportional to k-NN distance:
 *        w(x) = (ε + d_k(x))^{-α}
 *    where d_k(x) is distance to k-th nearest neighbor,
 *          ε = 10^{-10} (regularization),
 *          α = 1.5 (exponent in [1,2])
 *    Measure of a set: μ(A) = Σ_{x ∈ A} w(x)
 *    Points in dense regions (small d_k) receive large weights
 *    Valuable in moderate-to-high dimensions where kernel density estimation
 *    becomes unstable; provides robust local density surrogate without
 *    bandwidth selection
 *
 * COMPUTATIONAL COMPLEXITY:
 * The algorithm proceeds through several phases:
 *   - k-NN computation: O(n log n) using ANN library with kd-trees
 *   - Edge construction: O(n² k) worst case, typically O(nk²)
 *   - Triangle construction: O(m k) where m = number of edges
 *   - Higher simplices: Combinatorial explosion limits practical max_p
 *
 * Memory scales as O(nk + m) for complex structure, where m is typically O(nk).
 * For max_p > 2, memory requirements increase substantially with higher-
 * dimensional simplex counts.
 *
 * IMPLEMENTATION DETAILS:
 * This is a C interface function callable via .Call() from R. It handles:
 *   1. Extraction and validation of all SEXP parameters
 *   2. Conversion of R matrices (dense or dgCMatrix) to Eigen::SparseMatrix
 *   3. Conversion of R vectors to Eigen::VectorXd
 *   4. Construction of riem_dcx_t via build_knn_riem_dcx() method
 *   5. Wrapping result in external pointer with finalizer
 *   6. Setting class attribute for S4 dispatch in R
 *
 * The function performs extensive validation:
 *   - Matrix dimensions must be positive
 *   - k must satisfy 2 ≤ k < n
 *   - max_p must satisfy 1 ≤ max_p < n
 *   - If max_p > 10, issues warning about memory usage
 *   - Response vector length must match number of points (if provided)
 *   - density_normalization must be non-negative
 *
 * ADAPTIVE DIMENSION REDUCTION:
 * If no simplices of dimension p can be constructed (empty intersections),
 * the function automatically reduces max_p to p-1 and continues. This ensures
 * the complex terminates at the highest dimension where simplices exist.
 *
 * MATRIX FORMAT HANDLING:
 * Accepts both dense R matrices (REALSXP) and sparse matrices (dgCMatrix from
 * Matrix package). For dense matrices, converts to sparse format internally by
 * building triplet list and filtering zero entries. For sparse matrices,
 * extracts dgCMatrix slots (@i, @p, @x, @Dim) and rebuilds Eigen sparse matrix.
 *
 * ERROR HANDLING:
 * Uses R's error mechanism (Rf_error, Rf_warning) for all error conditions.
 * Ensures proper cleanup of allocated riem_dcx_t object if construction fails.
 * All C++ exceptions are caught and converted to R errors with descriptive
 * messages.
 *
 * MEMORY MANAGEMENT:
 * The riem_dcx_t object is allocated on the heap and wrapped in an R external
 * pointer with registered finalizer. When R garbage collects the external
 * pointer, the finalizer deletes the C++ object, preventing memory leaks.
 * The finalizer is registered with R_RegisterCFinalizerEx() with on_exit=TRUE
 * to ensure cleanup even if R session terminates.
 *
 * @param s_X SEXP matrix (REALSXP or S4 dgCMatrix) of dimension n × d where
 *            each row represents a point in d-dimensional space. Features
 *            should be appropriately scaled (consider standardization if
 *            features have different units).
 *
 * @param s_y SEXP numeric vector (REALSXP) of length n containing response
 *            values at each point, or R_NilValue for topology-only analysis.
 *            When provided, response is stored for subsequent signal processing
 *            (smoothing, regression with geometric regularization).
 *
 * @param s_k SEXP integer scalar giving number of nearest neighbors.
 *            Must satisfy 2 ≤ k < n. Larger k produces denser complexes but
 *            may smooth over fine-scale features. Smaller k captures local
 *            structure but may fragment complex. Typical values: k ∈ [5, 20].
 *            Rough guideline: k ≈ log(n) as starting point.
 *
 * @param s_max_p SEXP integer scalar giving maximum simplex dimension.
 *                Complex includes all simplices from dimension 0 (vertices)
 *                through max_p. Must satisfy 1 ≤ max_p < n. For most
 *                applications, max_p = 2 suffices. Higher dimensions increase
 *                computational cost dramatically and are rarely needed unless
 *                studying higher-dimensional topological features (e.g.,
 *                voids, cavities).
 *
 * @param s_use_counting_measure SEXP logical scalar. If TRUE, use counting
 *                                measure (unit weights). If FALSE, use density-
 *                                based weights from k-NN distances. See formula
 *                                in detailed description above.
 *
 * @param s_density_normalization SEXP numeric scalar, non-negative. Specifies
 *                                 target sum for normalized vertex weights.
 *                                 If 0 (default), weights normalized to sum
 *                                 to n (maintains counting-measure interpretation).
 *                                 If positive, weights sum to specified value.
 *                                 No effect when use_counting_measure = TRUE.
 *
 * @return SEXP external pointer (EXTPTRSXP) to riem_dcx_t object with class
 *         attribute "Rcpp_riem_dcx". The wrapped object contains:
 *           - Simplex tables S[0], ..., S[max_p] recording vertices
 *           - Metric matrices g.M[0], ..., g.M[max_p] encoding geometry
 *           - Boundary operators L.B[1], ..., L.B[max_p] for discrete calculus
 *           - Hodge Laplacians L.L[0], ..., L.L[max_p] for diffusion
 *           - Initial density distributions rho.rho[0], ..., rho.rho[max_p]
 *           - Signal state sig.y (response) and sig.y_hat_hist (fitted values)
 *           - Star tables stars[0], ..., stars[max_p-1] for upward incidence
 *
 *         This object can be passed to other package functions for:
 *           - Signal smoothing via heat diffusion (heat_apply)
 *           - Regression with geometric regularization (iteration_step)
 *           - Spectral analysis of Hodge Laplacian
 *           - Discrete Morse theory computations
 *
 * @note Registered in init.c as:
 *       {"S_build_nerve_from_knn", (DL_FUNC) &S_build_nerve_from_knn, 7}
 *       Called from R via:
 *       .Call("S_build_nerve_from_knn", X, y, k, max_p,
 *             use_counting_measure, density_normalization,
 *             PACKAGE = "gflow")
 *
 * @warning For large n or large k, simplex count grows rapidly, especially at
 *          higher dimensions. Function issues warning if max_p > 10. Monitor
 *          memory usage for large-scale problems.
 *
 * @warning Function automatically reduces max_p if no simplices of given
 *          dimension can be constructed. For example, if no triangles exist,
 *          complex terminates at dimension 1 regardless of requested max_p.
 *
 * @see riem_dcx_t::build_knn_riem_dcx() - Internal C++ construction method
 * @see S_riem_dcx_summary() - Extract summary information from complex
 * @see build.nerve.from.knn() - R wrapper function with validation
 *
 * @references
 * Carlsson, G. (2009). Topology and data. Bulletin of the American
 * Mathematical Society, 46(2), 255-308.
 *
 * Edelsbrunner, H., & Harer, J. (2010). Computational Topology:
 * An Introduction. American Mathematical Society.
 *
 * Lim, L. H. (2020). Hodge Laplacians on graphs. SIAM Review, 62(3), 685-715.
 *
 * @example
 * // From R:
 * // set.seed(123)
 * // n <- 100
 * // X <- cbind(rnorm(n), rnorm(n))
 * // y <- X[,1]^2 + X[,2]^2 + rnorm(n, sd = 0.1)
 * // dcx <- build.nerve.from.knn(X, y, k = 10, max_p = 2,
 * //                              use_counting_measure = TRUE,
 * //                              )
 * // summary(dcx)
 */
extern "C" SEXP S_build_nerve_from_knn(
    SEXP s_X,
    SEXP s_y,
    SEXP s_k,
    SEXP s_max_p,
    SEXP s_use_counting_measure,
    SEXP s_density_normalization
) {
    // ==================== Input Extraction ====================

    // Extract k
    if (!Rf_isInteger(s_k) || Rf_length(s_k) != 1) {
        Rf_error("k must be a single integer");
    }
    const int k = Rf_asInteger(s_k);

    // Extract max_p
    if (!Rf_isInteger(s_max_p) || Rf_length(s_max_p) != 1) {
        Rf_error("max_p must be a single integer");
    }
    const int max_p = Rf_asInteger(s_max_p);

    // Extract use_counting_measure
    if (!Rf_isLogical(s_use_counting_measure) || Rf_length(s_use_counting_measure) != 1) {
        Rf_error("use_counting_measure must be a single logical value");
    }
    const bool use_counting_measure = (bool)Rf_asLogical(s_use_counting_measure);

    // Extract density_normalization
    if (!Rf_isReal(s_density_normalization) || Rf_length(s_density_normalization) != 1) {
        Rf_error("density_normalization must be a single numeric value");
    }
    const double density_normalization = Rf_asReal(s_density_normalization);

    // ==================== Convert X to Eigen Sparse Matrix ====================
    Eigen::SparseMatrix<double> X_sparse;

    if (Rf_isMatrix(s_X) && TYPEOF(s_X) == REALSXP) {
        // Dense matrix input
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, R_DimSymbol));
        if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) != 2) {
            UNPROTECT(1);
            Rf_error("X must be a valid matrix with dim attribute");
        }

        const int n_rows = INTEGER(s_dim)[0];
        const int n_cols = INTEGER(s_dim)[1];
        UNPROTECT(1);

        const double* X_data = REAL(s_X);

        // Convert to sparse format
        X_sparse.resize(n_rows, n_cols);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(n_rows * n_cols / 2);

        for (int j = 0; j < n_cols; ++j) {
            for (int i = 0; i < n_rows; ++i) {
                double val = X_data[i + j * n_rows];
                if (val != 0.0) {
                    triplets.push_back(Eigen::Triplet<double>(i, j, val));
                }
            }
        }
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

    } else if (Rf_inherits(s_X, "dgCMatrix")) {
        // Sparse matrix input
        SEXP s_i = PROTECT(Rf_getAttrib(s_X, Rf_install("i")));
        SEXP s_p = PROTECT(Rf_getAttrib(s_X, Rf_install("p")));
        SEXP s_x = PROTECT(Rf_getAttrib(s_X, Rf_install("x")));
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, Rf_install("Dim")));

        if (s_i == R_NilValue || s_p == R_NilValue ||
            s_x == R_NilValue || s_dim == R_NilValue) {
            UNPROTECT(4);
            Rf_error("Invalid dgCMatrix: missing required slots");
        }

        const int* i_data = INTEGER(s_i);
        const int* p_data = INTEGER(s_p);
        const double* x_data = REAL(s_x);
        const int* dim_data = INTEGER(s_dim);

        const int n_rows = dim_data[0];
        const int n_cols = dim_data[1];
        const int nnz = Rf_length(s_x);

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(nnz);

        for (int j = 0; j < n_cols; ++j) {
            for (int idx = p_data[j]; idx < p_data[j + 1]; ++idx) {
                int i = i_data[idx];
                double val = x_data[idx];
                triplets.push_back(Eigen::Triplet<double>(i, j, val));
            }
        }

        X_sparse.resize(n_rows, n_cols);
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

        UNPROTECT(4);

    } else {
        Rf_error("X must be a numeric matrix or dgCMatrix");
    }

    const Eigen::Index n_points = X_sparse.rows();
    const Eigen::Index n_features = X_sparse.cols();

    // ==================== Convert y to Eigen Vector ====================
    Eigen::VectorXd y_vec;

    if (s_y == R_NilValue) {
        y_vec = Eigen::VectorXd::Zero(0);
    } else {
        if (!Rf_isReal(s_y)) {
            Rf_error("y must be a numeric vector or NULL");
        }
        const int y_len = Rf_length(s_y);
        if (y_len != n_points) {
            Rf_error("Length of y (%d) must equal number of rows in X (%ld)",
                     y_len, n_points);
        }

        const double* y_data = REAL(s_y);
        y_vec.resize(y_len);
        for (int i = 0; i < y_len; ++i) {
            y_vec[i] = y_data[i];
        }
    }

    // ==================== Validate Parameters ====================
    if (n_points <= 0 || n_features <= 0) {
        Rf_error("X must have positive dimensions (n_points=%ld, n_features=%ld)",
                 n_points, n_features);
    }

    if (k < 2) {
        Rf_error("k must be at least 2 (got k=%d)", k);
    }
    if (k >= n_points) {
        Rf_error("k must be less than n_points (got k=%d, n_points=%ld)",
                 k, n_points);
    }

    if (max_p < 1) {
        Rf_error("max_p must be at least 1 (got max_p=%d)", max_p);
    }
    if (max_p >= n_points) {
        Rf_error("max_p must be less than n_points (got max_p=%d, n_points=%ld)",
                 max_p, n_points);
    }
    if (max_p > 10) {
        Rf_warning("Large max_p=%d may result in excessive memory usage", max_p);
    }

    if (density_normalization < 0.0) {
        Rf_error("density_normalization must be non-negative (got %f)",
                 density_normalization);
    }

    // ==================== Build Complex ====================
    riem_dcx_t* dcx = nullptr;

    try {
        dcx = new riem_dcx_t();

        dcx->build_knn_riem_dcx(
            X_sparse,
            y_vec,
            static_cast<index_t>(k),
            static_cast<index_t>(max_p),
            use_counting_measure,
            density_normalization
        );

    } catch (const std::exception& e) {
        if (dcx != nullptr) delete dcx;
        Rf_error("Failed to build nerve complex: %s", e.what());
    } catch (...) {
        if (dcx != nullptr) delete dcx;
        Rf_error("Unknown error during nerve complex construction");
    }

    // ==================== Create External Pointer ====================
    // Create external pointer with finalizer
    SEXP ext_ptr = PROTECT(R_MakeExternalPtr(dcx, R_NilValue, R_NilValue));

    // Register finalizer
    R_RegisterCFinalizerEx(
        ext_ptr,
        [](SEXP ptr) {
            riem_dcx_t* obj = static_cast<riem_dcx_t*>(R_ExternalPtrAddr(ptr));
            if (obj != nullptr) {
                delete obj;
                R_ClearExternalPtr(ptr);
            }
        },
        TRUE
    );

    // Set class attribute
    SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(class_attr, 0, Rf_mkChar("Rcpp_riem_dcx"));
    Rf_setAttrib(ext_ptr, R_ClassSymbol, class_attr);

    UNPROTECT(2); // class_attr, ext_ptr

    return ext_ptr;
}


/**
 * @brief Extract summary from riem_dcx_t object
 *
 * Called via .Call() from R with manual registration in init.c
 *
 * @param s_dcx_ptr SEXP external pointer to riem_dcx_t
 * @return SEXP list with summary information
 */
extern "C" SEXP S_riem_dcx_summary(SEXP s_dcx_ptr) {
    // Validate pointer
    if (TYPEOF(s_dcx_ptr) != EXTPTRSXP) {
        Rf_error("Argument must be an external pointer");
    }

    riem_dcx_t* dcx = static_cast<riem_dcx_t*>(R_ExternalPtrAddr(s_dcx_ptr));
    if (dcx == nullptr) {
        Rf_error("Invalid or deleted riem_dcx_t pointer");
    }

    // Extract information
    const int max_dim = dcx->pmax;

    // Create n_simplices vector
    SEXP n_simplices = PROTECT(Rf_allocVector(INTSXP, max_dim + 1));
    int* n_simplices_data = INTEGER(n_simplices);

    for (int p = 0; p <= max_dim; ++p) {
        n_simplices_data[p] = static_cast<int>(dcx->S[p].size());
    }

    // Extract response information
    const bool has_response = (dcx->sig.y.size() > 0);
    const int response_length = static_cast<int>(dcx->sig.y.size());

    // Build result list
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 7));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));

    SET_STRING_ELT(names, 0, Rf_mkChar("max_dimension"));
    SET_STRING_ELT(names, 1, Rf_mkChar("n_simplices"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_vertices"));
    SET_STRING_ELT(names, 3, Rf_mkChar("n_edges"));
    SET_STRING_ELT(names, 4, Rf_mkChar("n_triangles"));
    SET_STRING_ELT(names, 5, Rf_mkChar("has_response"));
    SET_STRING_ELT(names, 6, Rf_mkChar("response_length"));

    SET_VECTOR_ELT(result, 0, Rf_ScalarInteger(max_dim));
    SET_VECTOR_ELT(result, 1, n_simplices);
    SET_VECTOR_ELT(result, 2, Rf_ScalarInteger(n_simplices_data[0]));
    SET_VECTOR_ELT(result, 3, Rf_ScalarInteger(max_dim >= 1 ? n_simplices_data[1] : 0));
    SET_VECTOR_ELT(result, 4, Rf_ScalarInteger(max_dim >= 2 ? n_simplices_data[2] : 0));
    SET_VECTOR_ELT(result, 5, Rf_ScalarLogical(has_response));
    SET_VECTOR_ELT(result, 6, Rf_ScalarInteger(response_length));

    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(3); // names, result, n_simplices

    return result;
}

/**
 * @brief Extract simplex vertex lists for a given dimension
 */
extern "C" SEXP S_get_simplices(SEXP s_dcx_ptr, SEXP s_dim) {
    if (TYPEOF(s_dcx_ptr) != EXTPTRSXP) {
        Rf_error("Argument must be an external pointer");
    }

    riem_dcx_t* dcx = static_cast<riem_dcx_t*>(R_ExternalPtrAddr(s_dcx_ptr));
    if (dcx == nullptr) {
        Rf_error("Invalid or deleted riem_dcx_t pointer");
    }

    int dim = Rf_asInteger(s_dim);
    if (dim < 0 || dim > dcx->pmax) {
        Rf_error("Invalid dimension %d (must be 0 to %d)", dim, dcx->pmax);
    }

    const auto& simplices = dcx->S[dim].simplex_verts;
    const index_t n_simplices = simplices.size();

    // Create list of integer vectors
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_simplices));

    for (index_t i = 0; i < n_simplices; ++i) {
        const auto& simplex = simplices[i];
        SEXP simplex_vec = PROTECT(Rf_allocVector(INTSXP, simplex.size()));
        int* p_simplex = INTEGER(simplex_vec);

        for (size_t j = 0; j < simplex.size(); ++j) {
            p_simplex[j] = static_cast<int>(simplex[j]) + 1; // R uses 1-indexing
        }

        SET_VECTOR_ELT(result, i, simplex_vec);
        UNPROTECT(1); // simplex_vec
    }

    UNPROTECT(1); // result
    return result;
}

/**
 * @brief Extract metric diagonal values for a given dimension
 */
extern "C" SEXP S_get_metric_diagonal(SEXP s_dcx_ptr, SEXP s_dim) {
    if (TYPEOF(s_dcx_ptr) != EXTPTRSXP) {
        Rf_error("Argument must be an external pointer");
    }

    riem_dcx_t* dcx = static_cast<riem_dcx_t*>(R_ExternalPtrAddr(s_dcx_ptr));
    if (dcx == nullptr) {
        Rf_error("Invalid or deleted riem_dcx_t pointer");
    }

    int dim = Rf_asInteger(s_dim);
    if (dim < 0 || dim > dcx->pmax) {
        Rf_error("Invalid dimension %d (must be 0 to %d)", dim, dcx->pmax);
    }

    const spmat_t& M = dcx->g.M[dim];
    const Eigen::Index n = M.rows();

    SEXP result = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_result = REAL(result);

    for (Eigen::Index i = 0; i < n; ++i) {
        p_result[i] = M.coeff(i, i);
    }

    UNPROTECT(1);
    return result;
}
