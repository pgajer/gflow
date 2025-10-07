#include "riem_dcx.hpp"
#include <R.h>
#include <Rinternals.h>
#include <cstring>

// Forward declare conversion utilities from SEXP_cpp_conversion_utils.cpp
// These should be available if that file is compiled into the package
namespace sexp_utils {
    Eigen::SparseMatrix<double> sexp_to_eigen_sparse(SEXP s_X);
    // Add other utility declarations as needed
}


// ================================================================
// HELPER FUNCTIONS TO BUILD NESTED COMPONENTS
// ================================================================

extern "C" SEXP create_graph_component(const riem_dcx_t& dcx) {
    const int n_fields = 6;
    SEXP graph = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // n.vertices
    SET_STRING_ELT(names, idx, Rf_mkChar("n.vertices"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.S[0].size()));

    // n.edges
    SET_STRING_ELT(names, idx, Rf_mkChar("n.edges"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.S[1].size()));

    // edge.lengths
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.lengths"));
    const size_t n_edges = dcx.edge_lengths.size();
    SEXP s_edge_lengths = PROTECT(Rf_allocVector(REALSXP, n_edges));
    for (size_t e = 0; e < n_edges; ++e) {
        REAL(s_edge_lengths)[e] = dcx.edge_lengths[e];
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_lengths);
    UNPROTECT(1);

    // vertex.densities
    SET_STRING_ELT(names, idx, Rf_mkChar("vertex.densities"));
    const Eigen::Index n_verts = dcx.rho.rho[0].size();
    SEXP s_vert_dens = PROTECT(Rf_allocVector(REALSXP, n_verts));
    for (Eigen::Index i = 0; i < n_verts; ++i) {
        REAL(s_vert_dens)[i] = dcx.rho.rho[0][i];
    }
    SET_VECTOR_ELT(graph, idx++, s_vert_dens);
    UNPROTECT(1);

    // edge.densities
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.densities"));
    const Eigen::Index n_edges_dens = dcx.rho.rho[1].size();
    SEXP s_edge_dens = PROTECT(Rf_allocVector(REALSXP, n_edges_dens));
    for (Eigen::Index e = 0; e < n_edges_dens; ++e) {
        REAL(s_edge_dens)[e] = dcx.rho.rho[1][e];
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_dens);
    UNPROTECT(1);

    // edge.list (n_edges × 2 matrix)
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.list"));
    SEXP s_edge_list = PROTECT(Rf_allocMatrix(INTSXP, n_edges, 2));
    int* edge_data = INTEGER(s_edge_list);
    for (size_t e = 0; e < n_edges; ++e) {
        const auto& verts = dcx.S[1].simplex_verts[e];
        edge_data[e] = verts[0] + 1;  // R uses 1-based indexing
        edge_data[e + n_edges] = verts[1] + 1;
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_list);
    UNPROTECT(1);

    Rf_setAttrib(graph, R_NamesSymbol, names);
    UNPROTECT(2); // names, graph
    return graph;
}

extern "C" SEXP create_iteration_component(const riem_dcx_t& dcx) {
    const int n_fields = 5;
    SEXP iteration = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // converged
    SET_STRING_ELT(names, idx, Rf_mkChar("converged"));
    SET_VECTOR_ELT(iteration, idx++, Rf_ScalarLogical(dcx.converged));

    // n.iterations
    SET_STRING_ELT(names, idx, Rf_mkChar("n.iterations"));
    SET_VECTOR_ELT(iteration, idx++, Rf_ScalarInteger(dcx.n_iterations));

    // response.changes
    SET_STRING_ELT(names, idx, Rf_mkChar("response.changes"));
    const size_t n_changes = dcx.response_changes.size();
    SEXP s_resp_changes = PROTECT(Rf_allocVector(REALSXP, n_changes));
    for (size_t i = 0; i < n_changes; ++i) {
        REAL(s_resp_changes)[i] = dcx.response_changes[i];
    }
    SET_VECTOR_ELT(iteration, idx++, s_resp_changes);
    UNPROTECT(1);

    // density.changes
    SET_STRING_ELT(names, idx, Rf_mkChar("density.changes"));
    const size_t n_dens_changes = dcx.density_changes.size();
    SEXP s_dens_changes = PROTECT(Rf_allocVector(REALSXP, n_dens_changes));
    for (size_t i = 0; i < n_dens_changes; ++i) {
        REAL(s_dens_changes)[i] = dcx.density_changes[i];
    }
    SET_VECTOR_ELT(iteration, idx++, s_dens_changes);
    UNPROTECT(1);

    // fitted.history (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("fitted.history"));
    const size_t n_hist = dcx.sig.y_hat_hist.size();
    SEXP s_history = PROTECT(Rf_allocVector(VECSXP, n_hist));
    for (size_t iter = 0; iter < n_hist; ++iter) {
        const vec_t& y_hat_iter = dcx.sig.y_hat_hist[iter];
        const Eigen::Index n = y_hat_iter.size();
        SEXP s_y_hat = PROTECT(Rf_allocVector(REALSXP, n));
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_hat)[i] = y_hat_iter[i];
        }
        SET_VECTOR_ELT(s_history, iter, s_y_hat);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(iteration, idx++, s_history);
    UNPROTECT(1);

    Rf_setAttrib(iteration, R_NamesSymbol, names);
    UNPROTECT(2); // names, iteration
    return iteration;
}

extern "C" SEXP create_parameters_component(
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double t_diffusion,
    double beta_damping,
    double gamma_modulation,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations
) {
    const int n_fields = 11;
    SEXP params = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    SET_STRING_ELT(names, idx, Rf_mkChar("k"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(k));

    SET_STRING_ELT(names, idx, Rf_mkChar("use.counting.measure"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarLogical(use_counting_measure));

    SET_STRING_ELT(names, idx, Rf_mkChar("density.normalization"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_normalization));

    SET_STRING_ELT(names, idx, Rf_mkChar("t.diffusion"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(t_diffusion));

    SET_STRING_ELT(names, idx, Rf_mkChar("beta.damping"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(beta_damping));

    SET_STRING_ELT(names, idx, Rf_mkChar("gamma.modulation"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(gamma_modulation));

    SET_STRING_ELT(names, idx, Rf_mkChar("n.eigenpairs"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(n_eigenpairs));

    SET_STRING_ELT(names, idx, Rf_mkChar("filter.type"));
    const char* filter_str = (filter_type == rdcx_filter_type_t::HEAT_KERNEL) ? "heat_kernel" :
                             (filter_type == rdcx_filter_type_t::TIKHONOV) ? "tikhonov" :
                             "cubic_spline";
    SET_VECTOR_ELT(params, idx++, Rf_mkString(filter_str));

    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.y"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_y));

    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.rho"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_rho));

    SET_STRING_ELT(names, idx, Rf_mkChar("max.iterations"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(max_iterations));

    Rf_setAttrib(params, R_NamesSymbol, names);
    UNPROTECT(2); // names, params
    return params;
}

/**
 * @brief R interface for kNN Riemannian graph regression
 *
 * Constructs 1-skeleton kNN complex and iteratively refines geometry to
 * reflect response structure for conditional expectation estimation.
 *
 * @param s_X Feature matrix (dense REALSXP or sparse dgCMatrix)
 * @param s_y Response vector (REALSXP)
 * @param s_k Number of nearest neighbors (INTSXP)
 * @param s_use_counting_measure Logical: use counting measure? (LGLSXP)
 * @param s_density_normalization Target density sum (REALSXP)
 * @param s_t_diffusion Heat diffusion time, 0=auto (REALSXP)
 * @param s_beta_damping Damping parameter, 0=auto (REALSXP)
 * @param s_gamma_modulation Response coherence exponent (REALSXP)
 * @param s_n_eigenpairs Number of eigenpairs for filtering (INTSXP)
 * @param s_filter_type Filter type string (STRSXP)
 * @param s_epsilon_y Response convergence threshold (REALSXP)
 * @param s_epsilon_rho Density convergence threshold (REALSXP)
 * @param s_max_iterations Maximum iteration count (INTSXP)
 *
 * @return External pointer to fitted riem_dcx_t object with class attribute
 *
 * @note This function is called from R via .Call(). Input validation is
 *       performed on the R side in fit.knn.riem.graph.regression().
 *       Additional defensive checks are included here for robustness.
 */
extern "C" SEXP S_fit_knn_riem_graph_regression(
    SEXP s_X,
    SEXP s_y,
    SEXP s_k,
    SEXP s_use_counting_measure,
    SEXP s_density_normalization,
    SEXP s_t_diffusion,
    SEXP s_beta_damping,
    SEXP s_gamma_modulation,
    SEXP s_n_eigenpairs,
    SEXP s_filter_type,
    SEXP s_epsilon_y,
    SEXP s_epsilon_rho,
    SEXP s_max_iterations,
    SEXP s_max_ratio_threshold,
    SEXP s_threshold_percentile,
    SEXP s_test_stage
) {
    // ================================================================
    // PART I: INPUT EXTRACTION (same as before)
    // ================================================================


    // -------------------- Feature Matrix X --------------------

    Eigen::SparseMatrix<double> X_sparse;
    Eigen::Index n_points = 0;
    Eigen::Index n_features = 0;

    // Check if dense matrix
    bool is_dense = (Rf_isMatrix(s_X) && TYPEOF(s_X) == REALSXP);
    bool is_sparse = Rf_inherits(s_X, "dgCMatrix");

    if (is_dense) {
        // Dense matrix: convert to sparse
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, R_DimSymbol));

        if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) != 2) {
            UNPROTECT(1);
            Rf_error("X must be a valid matrix with dim attribute");
        }

        n_points = INTEGER(s_dim)[0];
        n_features = INTEGER(s_dim)[1];
        UNPROTECT(1);

        if (n_points < 1 || n_features < 1) {
            Rf_error("X has invalid dimensions: %ld × %ld",
                     (long)n_points, (long)n_features);
        }

        const double* X_data = REAL(s_X);

        // Convert to sparse format
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(n_points * n_features / 2); // Heuristic reserve

        for (Eigen::Index j = 0; j < n_features; ++j) {
            for (Eigen::Index i = 0; i < n_points; ++i) {
                double val = X_data[i + j * n_points];
                if (val != 0.0) {
                    triplets.emplace_back(i, j, val);
                }
            }
        }

        X_sparse.resize(n_points, n_features);
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

    } else if (is_sparse) {
        // Sparse matrix (dgCMatrix): extract components
        SEXP s_i = PROTECT(Rf_getAttrib(s_X, Rf_install("i")));
        SEXP s_p = PROTECT(Rf_getAttrib(s_X, Rf_install("p")));
        SEXP s_x = PROTECT(Rf_getAttrib(s_X, Rf_install("x")));
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, Rf_install("Dim")));

        if (s_i == R_NilValue || s_p == R_NilValue ||
            s_x == R_NilValue || s_dim == R_NilValue) {
            UNPROTECT(4);
            Rf_error("Invalid dgCMatrix: missing required slots (i, p, x, or Dim)");
        }

        if (TYPEOF(s_i) != INTSXP || TYPEOF(s_p) != INTSXP ||
            TYPEOF(s_x) != REALSXP || TYPEOF(s_dim) != INTSXP) {
            UNPROTECT(4);
            Rf_error("Invalid dgCMatrix: slots have incorrect types");
        }

        const int* i_data = INTEGER(s_i);
        const int* p_data = INTEGER(s_p);
        const double* x_data = REAL(s_x);
        const int* dim_data = INTEGER(s_dim);

        n_points = dim_data[0];
        n_features = dim_data[1];

        if (n_points < 1 || n_features < 1) {
            UNPROTECT(4);
            Rf_error("dgCMatrix has invalid dimensions: %ld × %ld",
                     (long)n_points, (long)n_features);
        }

        if (Rf_length(s_dim) != 2) {
            UNPROTECT(4);
            Rf_error("dgCMatrix Dim slot must have length 2");
        }

        if (Rf_length(s_p) != n_features + 1) {
            UNPROTECT(4);
            Rf_error("dgCMatrix p slot has incorrect length");
        }

        const int nnz = Rf_length(s_x);

        // Build triplet list
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(nnz);

        for (int j = 0; j < n_features; ++j) {
            int col_start = p_data[j];
            int col_end = p_data[j + 1];

            if (col_start < 0 || col_end > nnz || col_start > col_end) {
                UNPROTECT(4);
                Rf_error("dgCMatrix has invalid p slot (column pointers)");
            }

            for (int idx = col_start; idx < col_end; ++idx) {
                int i = i_data[idx];

                if (i < 0 || i >= n_points) {
                    UNPROTECT(4);
                    Rf_error("dgCMatrix has invalid row index: %d (must be in [0, %ld))",
                             i, (long)n_points);
                }

                double val = x_data[idx];
                triplets.emplace_back(i, j, val);
            }
        }

        X_sparse.resize(n_points, n_features);
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

        UNPROTECT(4); // s_i, s_p, s_x, s_dim

    } else {
        Rf_error("X must be a numeric matrix or dgCMatrix (from Matrix package)");
    }

    // -------------------- Response Vector y --------------------

    if (TYPEOF(s_y) != REALSXP) {
        Rf_error("y must be a numeric (REAL) vector");
    }

    const R_xlen_t y_len = Rf_xlength(s_y);

    if (y_len != n_points) {
        Rf_error("Length of y (%ld) must equal number of rows in X (%ld)",
                 (long)y_len, (long)n_points);
    }

    vec_t y(y_len);
    const double* y_data = REAL(s_y);

    for (R_xlen_t i = 0; i < y_len; ++i) {
        y[i] = y_data[i];
    }

    // -------------------- Parameter k --------------------

    if (TYPEOF(s_k) != INTSXP || Rf_length(s_k) != 1) {
        Rf_error("k must be a single integer");
    }

    const int k_raw = INTEGER(s_k)[0];

    if (k_raw == NA_INTEGER) {
        Rf_error("k cannot be NA");
    }

    if (k_raw < 2 || k_raw >= n_points) {
        Rf_error("k must satisfy 2 <= k < n (got k=%d, n=%ld)",
                 k_raw, (long)n_points);
    }

    const index_t k = static_cast<index_t>(k_raw);

    // -------------------- use.counting.measure --------------------

    if (TYPEOF(s_use_counting_measure) != LGLSXP ||
        Rf_length(s_use_counting_measure) != 1) {
        Rf_error("use.counting.measure must be a single logical value");
    }

    const int ucm_raw = LOGICAL(s_use_counting_measure)[0];

    if (ucm_raw == NA_LOGICAL) {
        Rf_error("use.counting.measure cannot be NA");
    }

    const bool use_counting_measure = (ucm_raw != 0);

    // -------------------- density.normalization --------------------

    if (TYPEOF(s_density_normalization) != REALSXP ||
        Rf_length(s_density_normalization) != 1) {
        Rf_error("density.normalization must be a single numeric value");
    }

    const double density_normalization = REAL(s_density_normalization)[0];

    if (!R_FINITE(density_normalization) || density_normalization < 0.0) {
        Rf_error("density.normalization must be a finite non-negative number (got %.3f)",
                 density_normalization);
    }

    // -------------------- t.diffusion --------------------

    if (TYPEOF(s_t_diffusion) != REALSXP || Rf_length(s_t_diffusion) != 1) {
        Rf_error("t.diffusion must be a single numeric value");
    }

    const double t_diffusion = REAL(s_t_diffusion)[0];

    if (!R_FINITE(t_diffusion) || t_diffusion < 0.0) {
        Rf_error("t.diffusion must be a finite non-negative number (got %.3f)",
                 t_diffusion);
    }

    // -------------------- beta.damping --------------------

    if (TYPEOF(s_beta_damping) != REALSXP || Rf_length(s_beta_damping) != 1) {
        Rf_error("beta.damping must be a single numeric value");
    }

    const double beta_damping = REAL(s_beta_damping)[0];

    if (!R_FINITE(beta_damping) || beta_damping < 0.0) {
        Rf_error("beta.damping must be a finite non-negative number (got %.3f)",
                 beta_damping);
    }

    // -------------------- gamma.modulation --------------------

    if (TYPEOF(s_gamma_modulation) != REALSXP ||
        Rf_length(s_gamma_modulation) != 1) {
        Rf_error("gamma.modulation must be a single numeric value");
    }

    const double gamma_modulation = REAL(s_gamma_modulation)[0];

    if (!R_FINITE(gamma_modulation) || gamma_modulation <= 0.0) {
        Rf_error("gamma.modulation must be a finite positive number (got %.3f)",
                 gamma_modulation);
    }

    // -------------------- n.eigenpairs --------------------

    if (TYPEOF(s_n_eigenpairs) != INTSXP || Rf_length(s_n_eigenpairs) != 1) {
        Rf_error("n.eigenpairs must be a single integer");
    }

    const int n_eigenpairs_raw = INTEGER(s_n_eigenpairs)[0];

    if (n_eigenpairs_raw == NA_INTEGER) {
        Rf_error("n.eigenpairs cannot be NA");
    }

    if (n_eigenpairs_raw < 10 || n_eigenpairs_raw > n_points) {
        Rf_error("n.eigenpairs must satisfy 10 <= n.eigenpairs <= n (got %d, n=%ld)",
                 n_eigenpairs_raw, (long)n_points);
    }

    const int n_eigenpairs = n_eigenpairs_raw;

    // -------------------- filter.type --------------------

    if (TYPEOF(s_filter_type) != STRSXP || Rf_length(s_filter_type) != 1) {
        Rf_error("filter.type must be a single string");
    }

    const char* filter_str = CHAR(STRING_ELT(s_filter_type, 0));

    if (filter_str == nullptr || strlen(filter_str) == 0) {
        Rf_error("filter.type cannot be empty string");
    }

    rdcx_filter_type_t filter_type;

    if (strcmp(filter_str, "heat_kernel") == 0) {
        filter_type = rdcx_filter_type_t::HEAT_KERNEL;
    } else if (strcmp(filter_str, "tikhonov") == 0) {
        filter_type = rdcx_filter_type_t::TIKHONOV;
    } else if (strcmp(filter_str, "cubic_spline") == 0) {
        filter_type = rdcx_filter_type_t::CUBIC_SPLINE;
    } else {
        Rf_error("filter.type must be 'heat_kernel', 'tikhonov', or 'cubic_spline' (got '%s')",
                 filter_str);
    }

    // -------------------- epsilon.y --------------------

    if (TYPEOF(s_epsilon_y) != REALSXP || Rf_length(s_epsilon_y) != 1) {
        Rf_error("epsilon.y must be a single numeric value");
    }

    const double epsilon_y = REAL(s_epsilon_y)[0];

    if (!R_FINITE(epsilon_y) || epsilon_y <= 0.0) {
        Rf_error("epsilon.y must be a finite positive number (got %.3e)",
                 epsilon_y);
    }

    // -------------------- epsilon.rho --------------------

    if (TYPEOF(s_epsilon_rho) != REALSXP || Rf_length(s_epsilon_rho) != 1) {
        Rf_error("epsilon.rho must be a single numeric value");
    }

    const double epsilon_rho = REAL(s_epsilon_rho)[0];

    if (!R_FINITE(epsilon_rho) || epsilon_rho <= 0.0) {
        Rf_error("epsilon.rho must be a finite positive number (got %.3e)",
                 epsilon_rho);
    }

    // -------------------- max.iterations --------------------

    if (TYPEOF(s_max_iterations) != INTSXP || Rf_length(s_max_iterations) != 1) {
        Rf_error("max.iterations must be a single integer");
    }

    const int max_iterations = INTEGER(s_max_iterations)[0];

    if (max_iterations == NA_INTEGER) {
        Rf_error("max.iterations cannot be NA");
    }

    if (max_iterations < 1) {
        Rf_error("max.iterations must be at least 1 (got %d)", max_iterations);
    }

    // -------------------- max.ratio.threshold --------------------

    if (TYPEOF(s_max_ratio_threshold) != REALSXP || Rf_length(s_max_ratio_threshold) != 1) {
        Rf_error("max.ratio.threshold must be a single numeric value");
    }

    const double max_ratio_threshold = REAL(s_max_ratio_threshold)[0];

    if (!R_FINITE(max_ratio_threshold) || max_ratio_threshold <= 0.0) {
        Rf_error("max.ratio.threshold must be a finite positive number (got %.3e)",
                 max_ratio_threshold);
    }

    // -------------------- threshold.percentile --------------------

    if (TYPEOF(s_threshold_percentile) != REALSXP || Rf_length(s_threshold_percentile) != 1) {
        Rf_error("threshold.percentile must be a single numeric value");
    }

    const double threshold_percentile = REAL(s_threshold_percentile)[0];

    if (!R_FINITE(threshold_percentile) || threshold_percentile <= 0.0) {
        Rf_error("threshold.percentile must be a finite positive number (got %.3e)",
                 threshold_percentile);
    }

    // -------------------- s_test_stage --------------------

    if (TYPEOF(s_test_stage) != INTSXP || Rf_length(s_test_stage) != 1) {
        Rf_error("test_stage must be a single integer");
    }

    const int test_stage = INTEGER(s_test_stage)[0];

    if (test_stage == NA_INTEGER) {
        Rf_error("test_stage cannot be NA");
    }

    if (test_stage < -2) {
        Rf_error("test_stage must be at least -1 (got %d)", test_stage);
    }

    // ================================================================
    // PART II: CALL MEMBER FUNCTION
    // ================================================================

    riem_dcx_t dcx;  // Stack allocation now! No need for new/delete

    try {
        dcx.fit_knn_riem_graph_regression(
            X_sparse, y, k,
            use_counting_measure, density_normalization,
            t_diffusion, beta_damping, gamma_modulation,
            n_eigenpairs, filter_type,
            epsilon_y, epsilon_rho, max_iterations,
            max_ratio_threshold, threshold_percentile,
            test_stage
            );

    } catch (const std::exception& e) {
        Rf_error("Regression fitting failed: %s", e.what());
    }

    // ================================================================
    // PART III: BUILD RESULT LIST
    // ================================================================

    const Eigen::Index n = y.size();

    // ---------- Handle test stages with early termination ----------
    // When test_stage >= 0, the fit may have terminated early before
    // computing fitted values. We handle this by providing appropriate
    // placeholder values and diagnostic information.

    bool has_fitted_values = !dcx.sig.y_hat_hist.empty();
    vec_t y_hat_final;

    if (has_fitted_values) {
        y_hat_final = dcx.sig.y_hat_hist.back();
    } else {
        // Early termination: use observed y as placeholder
        // This allows the result structure to be consistent
        y_hat_final = y;

        // Warn user if this wasn't intentional
        if (test_stage < 0) {
            Rf_warning("No fitted values computed (unexpected early termination). "
                       "Returning observed y values as placeholder.");
        }
    }

    // Get final fitted values
    // if (dcx.sig.y_hat_hist.empty()) {
    //     Rf_error("No fitted values computed (internal error)");
    // }
    // const vec_t& y_hat_final = dcx.sig.y_hat_hist.back();

    // ---------- Main result list ----------

    const int n_components = 6;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_components));
    int component_idx = 0;

    // Component 1: fitted.values
    SET_STRING_ELT(names, component_idx, Rf_mkChar("fitted.values"));
    SEXP s_fitted = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_fitted)[i] = y_hat_final[i];
    }
    SET_VECTOR_ELT(result, component_idx++, s_fitted);
    UNPROTECT(1);

    // Component 2: residuals
    SET_STRING_ELT(names, component_idx, Rf_mkChar("residuals"));
    SEXP s_resid = PROTECT(Rf_allocVector(REALSXP, n));
    // Safe access: both y_hat_final and y are guaranteed to have size n
    for (Eigen::Index i = 0; i < n; ++i) {
        double residual;
        if (dcx.sig.y.size() == n) {
            residual = dcx.sig.y[i] - y_hat_final[i];
        } else {
            // Early termination: use input y
            residual = y[i] - y_hat_final[i];
        }
        REAL(s_resid)[i] = residual;
    }
    SET_VECTOR_ELT(result, component_idx++, s_resid);
    UNPROTECT(1);

    // Component 3: graph (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("graph"));
    SEXP s_graph = PROTECT(create_graph_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_graph);
    UNPROTECT(1);

    // Component 4: iteration (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("iteration"));
    SEXP s_iteration = PROTECT(create_iteration_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_iteration);
    UNPROTECT(1);

    // Component 5: parameters (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("parameters"));
    SEXP s_params = PROTECT(create_parameters_component(
        k, use_counting_measure, density_normalization,
        t_diffusion, beta_damping, gamma_modulation,
        n_eigenpairs, filter_type, epsilon_y, epsilon_rho, max_iterations
    ));
    SET_VECTOR_ELT(result, component_idx++, s_params);
    UNPROTECT(1);

    // Component 6: y (original response)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("y"));
    SEXP s_y_copy = PROTECT(Rf_allocVector(REALSXP, n));

    // CRITICAL FIX: Check if sig.y was populated during fitting
    if (dcx.sig.y.size() == n) {
        // Use the stored response from dcx
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_copy)[i] = dcx.sig.y[i];
        }
    } else {
        // Early termination: sig.y not populated, use input y
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_copy)[i] = y[i];
        }
    }
    SET_VECTOR_ELT(result, component_idx++, s_y_copy);
    UNPROTECT(1);

    // Set names attribute
    Rf_setAttrib(result, R_NamesSymbol, names);

    // Set class attribute
    SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(class_attr, 0, Rf_mkChar("knn.riem.fit"));
    SET_STRING_ELT(class_attr, 1, Rf_mkChar("list"));
    Rf_setAttrib(result, R_ClassSymbol, class_attr);

    UNPROTECT(3); // class_attr, names, result

    return result;
}
