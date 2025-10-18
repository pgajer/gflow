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

extern "C" SEXP create_density_history_component(const riem_dcx_t& dcx) {
    const size_t n_iters = dcx.density_history.rho_vertex.size();

    if (n_iters == 0) {
        // Return empty list if no history
        SEXP empty = PROTECT(Rf_allocVector(VECSXP, 0));
        UNPROTECT(1);
        return empty;
    }

    const int n_fields = 2;
    SEXP density = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // Field 1: rho.vertex (list of vectors, one per iteration)
    SET_STRING_ELT(names, idx, Rf_mkChar("rho.vertex"));
    SEXP s_rho_list = PROTECT(Rf_allocVector(VECSXP, n_iters));

    for (size_t i = 0; i < n_iters; ++i) {
        const vec_t& rho = dcx.density_history.rho_vertex[i];
        const Eigen::Index n_vertices = rho.size();

        SEXP s_rho_vec = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        for (Eigen::Index j = 0; j < n_vertices; ++j) {
            REAL(s_rho_vec)[j] = rho[j];
        }
        SET_VECTOR_ELT(s_rho_list, i, s_rho_vec);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(density, idx++, s_rho_list);
    UNPROTECT(1);

    // Field 2: n.iterations (scalar, for convenience)
    SET_STRING_ELT(names, idx, Rf_mkChar("n.iterations"));
    SET_VECTOR_ELT(density, idx++, Rf_ScalarInteger(n_iters));

    Rf_setAttrib(density, R_NamesSymbol, names);
    UNPROTECT(2); // names, density
    return density;
}

extern "C" SEXP create_gcv_component(const riem_dcx_t& dcx) {
    const size_t n_iters = dcx.gcv_history.iterations.size();

    const int n_fields = 4;
    SEXP gcv = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // eta.optimal (vector)
    SET_STRING_ELT(names, idx, Rf_mkChar("eta.optimal"));
    SEXP s_eta_opt = PROTECT(Rf_allocVector(REALSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        REAL(s_eta_opt)[i] = dcx.gcv_history.iterations[i].eta_optimal;
    }
    SET_VECTOR_ELT(gcv, idx++, s_eta_opt);
    UNPROTECT(1);

    // gcv.optimal (vector)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.optimal"));
    SEXP s_gcv_opt = PROTECT(Rf_allocVector(REALSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        REAL(s_gcv_opt)[i] = dcx.gcv_history.iterations[i].gcv_optimal;
    }
    SET_VECTOR_ELT(gcv, idx++, s_gcv_opt);
    UNPROTECT(1);

    // eta.grid (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("eta.grid"));
    SEXP s_eta_grid_list = PROTECT(Rf_allocVector(VECSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        const auto& grid = dcx.gcv_history.iterations[i].eta_grid;
        SEXP s_grid = PROTECT(Rf_allocVector(REALSXP, grid.size()));
        for (size_t j = 0; j < grid.size(); ++j) {
            REAL(s_grid)[j] = grid[j];
        }
        SET_VECTOR_ELT(s_eta_grid_list, i, s_grid);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(gcv, idx++, s_eta_grid_list);
    UNPROTECT(1);

    // gcv.scores (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.scores"));
    SEXP s_gcv_scores_list = PROTECT(Rf_allocVector(VECSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        const auto& scores = dcx.gcv_history.iterations[i].gcv_scores;
        SEXP s_scores = PROTECT(Rf_allocVector(REALSXP, scores.size()));
        for (size_t j = 0; j < scores.size(); ++j) {
            REAL(s_scores)[j] = scores[j];
        }
        SET_VECTOR_ELT(s_gcv_scores_list, i, s_scores);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(gcv, idx++, s_gcv_scores_list);
    UNPROTECT(1);

    Rf_setAttrib(gcv, R_NamesSymbol, names);
    UNPROTECT(2); // names, gcv
    return gcv;
}

extern "C" SEXP create_graph_component(const riem_dcx_t& dcx) {
    const int n_fields = 8;  // Expanded from 6 to 8
    SEXP graph = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // n.vertices
    SET_STRING_ELT(names, idx, Rf_mkChar("n.vertices"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.vertex_cofaces.size()));

    // n.edges
    SET_STRING_ELT(names, idx, Rf_mkChar("n.edges"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.edge_registry.size()));

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
    const Eigen::Index n_verts = dcx.vertex_cofaces.size();
    SEXP s_vert_dens = PROTECT(Rf_allocVector(REALSXP, n_verts));
    for (Eigen::Index i = 0; i < n_verts; ++i) {
        REAL(s_vert_dens)[i] = dcx.vertex_cofaces[i][0].density;
    }
    SET_VECTOR_ELT(graph, idx++, s_vert_dens);
    UNPROTECT(1);

    // edge.densities
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.densities"));
    const Eigen::Index n_edges_dens = dcx.edge_registry.size();
    SEXP s_edge_dens = PROTECT(Rf_allocVector(REALSXP, n_edges_dens));
    for (Eigen::Index e = 0; e < n_edges_dens; ++e) {
        // Find edge density from first vertex of the edge
        const auto [v0, v1] = dcx.edge_registry[e];
        double edge_density = 0.0;
        for (size_t k = 1; k < dcx.vertex_cofaces[v0].size(); ++k) {
            if (dcx.vertex_cofaces[v0][k].vertex_index == v1) {
                edge_density = dcx.vertex_cofaces[v0][k].density;
                break;
            }
        }
        REAL(s_edge_dens)[e] = edge_density;
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_dens);
    UNPROTECT(1);

    // edge.list (n_edges × 2 matrix)
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.list"));
    SEXP s_edge_list = PROTECT(Rf_allocMatrix(INTSXP, n_edges, 2));
    int* edge_data = INTEGER(s_edge_list);
    for (size_t e = 0; e < n_edges; ++e) {
        const auto [v0, v1] = dcx.edge_registry[e];
        edge_data[e] = v0 + 1;  // R uses 1-based indexing
        edge_data[e + n_edges] = v1 + 1;
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_list);
    UNPROTECT(1);

    // ---- adjacency list and edge lengths by vertex ----

    // adj.list (list of integer vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("adj.list"));
    SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, static_cast<size_t>(n_verts)));

    for (Eigen::Index i = 0; i < n_verts; ++i) {
        const auto& nbhrs = dcx.vertex_cofaces[i];
        // Note: nbhrs[0] corresponds to vertex i itself, so we skip it
        const size_t n_neighbors = nbhrs.size() - 1;

        SEXP RA = PROTECT(Rf_allocVector(INTSXP, (R_len_t)n_neighbors));
        int* A  = INTEGER(RA);
        for (size_t j = 1; j < nbhrs.size(); ++j) {
            *A++ = (int)nbhrs[j].vertex_index + 1;  // +1 for R's 1-based indexing
        }
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1); // RA
    }
    SET_VECTOR_ELT(graph, idx++, adj_list);
    UNPROTECT(1); // adj_list

    // edge.length.list (list of numeric vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.length.list"));
    SEXP edge_length_list = PROTECT(Rf_allocVector(VECSXP, static_cast<size_t>(n_verts)));

    for (Eigen::Index i = 0; i < n_verts; ++i) {
        const auto& nbhrs = dcx.vertex_cofaces[i];
        const size_t n_neighbors = nbhrs.size() - 1;

        SEXP RD = PROTECT(Rf_allocVector(REALSXP, (R_len_t)n_neighbors));
        double* D = REAL(RD);
        for (size_t j = 1; j < nbhrs.size(); ++j) {
            *D++ = nbhrs[j].dist;
        }
        SET_VECTOR_ELT(edge_length_list, i, RD);
        UNPROTECT(1); // RD
    }
    SET_VECTOR_ELT(graph, idx++, edge_length_list);
    UNPROTECT(1); // edge_length_list

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
    double t_scale_factor,
    double beta_coefficient_factor,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations,
    double density_alpha,
    double density_epsilon
    ) {

    const int n_fields = 15;
    SEXP params = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // k
    SET_STRING_ELT(names, idx, Rf_mkChar("k"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(k));

    // use.counting.measure
    SET_STRING_ELT(names, idx, Rf_mkChar("use.counting.measure"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarLogical(use_counting_measure));

    // density.normalization
    SET_STRING_ELT(names, idx, Rf_mkChar("density.normalization"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_normalization));

    // t.diffusion
    SET_STRING_ELT(names, idx, Rf_mkChar("t.diffusion"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(t_diffusion));

    // beta.damping (formerly density.uniform.pull)
    SET_STRING_ELT(names, idx, Rf_mkChar("beta.damping"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(beta_damping));

    // response.penalty.exp (formerly gamma.modulation)
    SET_STRING_ELT(names, idx, Rf_mkChar("response.penalty.exp"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(gamma_modulation));

    // NEW: t.scale.factor
    SET_STRING_ELT(names, idx, Rf_mkChar("t.scale.factor"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(t_scale_factor));

    // NEW: beta.coefficient.factor
    SET_STRING_ELT(names, idx, Rf_mkChar("beta.coefficient.factor"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(beta_coefficient_factor));

    // n.eigenpairs
    SET_STRING_ELT(names, idx, Rf_mkChar("n.eigenpairs"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(n_eigenpairs));

    // filter.type
    SET_STRING_ELT(names, idx, Rf_mkChar("filter.type"));

    const char* filter_str;
    switch (filter_type) {
    case rdcx_filter_type_t::HEAT_KERNEL:
        filter_str = "heat_kernel";
        break;
    case rdcx_filter_type_t::TIKHONOV:
        filter_str = "tikhonov";
        break;
    case rdcx_filter_type_t::CUBIC_SPLINE:
        filter_str = "cubic_spline";
        break;
    case rdcx_filter_type_t::GAUSSIAN:
        filter_str = "gaussian";
        break;
    case rdcx_filter_type_t::EXPONENTIAL:
        filter_str = "exponential";
        break;
    case rdcx_filter_type_t::BUTTERWORTH:
        filter_str = "butterworth";
        break;
    default:
        filter_str = "unknown";
        break;
    }

    SET_VECTOR_ELT(params, idx++, Rf_mkString(filter_str));

    // epsilon.y
    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.y"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_y));

    // epsilon.rho
    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.rho"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_rho));

    // max.iterations
    SET_STRING_ELT(names, idx, Rf_mkChar("max.iterations"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(max_iterations));

    // density.alpha
    SET_STRING_ELT(names, idx, Rf_mkChar("density.alpha"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_alpha));

    // density.epsilon
    SET_STRING_ELT(names, idx, Rf_mkChar("density.epsilon"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_epsilon));


    Rf_setAttrib(params, R_NamesSymbol, names);
    UNPROTECT(2); // names, params
    return params;
}

extern "C" SEXP create_gamma_selection_component(const riem_dcx_t& dcx) {
    // If gamma selection was not performed, return NULL
    if (!dcx.gamma_was_auto_selected) {
        return R_NilValue;
    }

    const auto& gsr = dcx.gamma_selection_result;

    const int n_fields = 5;
    SEXP gamma_sel = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // Field 1: gamma.optimal
    SET_STRING_ELT(names, idx, Rf_mkChar("gamma.optimal"));
    SET_VECTOR_ELT(gamma_sel, idx++, Rf_ScalarReal(gsr.gamma_optimal));

    // Field 2: gcv.optimal (FIRST-ITERATION GCV, not global optimal)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.first.iter"));
    SET_VECTOR_ELT(gamma_sel, idx++, Rf_ScalarReal(gsr.gcv_optimal));

    // Field 3: gamma.grid
    SET_STRING_ELT(names, idx, Rf_mkChar("gamma.grid"));
    SEXP s_gamma_grid = PROTECT(Rf_allocVector(REALSXP, gsr.gamma_grid.size()));
    for (size_t i = 0; i < gsr.gamma_grid.size(); ++i) {
        REAL(s_gamma_grid)[i] = gsr.gamma_grid[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_gamma_grid);
    UNPROTECT(1);

    // Field 4: gcv.scores (GCV at first iteration for each gamma)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.scores"));
    SEXP s_gcv_scores = PROTECT(Rf_allocVector(REALSXP, gsr.gcv_scores.size()));
    for (size_t i = 0; i < gsr.gcv_scores.size(); ++i) {
        REAL(s_gcv_scores)[i] = gsr.gcv_scores[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_gcv_scores);
    UNPROTECT(1);

    // Field 5: y.hat.optimal (fitted values at first iteration with optimal gamma)
    SET_STRING_ELT(names, idx, Rf_mkChar("y.hat.first.iter"));
    const Eigen::Index n = gsr.y_hat_optimal.size();
    SEXP s_y_hat_opt = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_y_hat_opt)[i] = gsr.y_hat_optimal[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_y_hat_opt);
    UNPROTECT(1);

    Rf_setAttrib(gamma_sel, R_NamesSymbol, names);
    UNPROTECT(2); // names, gamma_sel
    return gamma_sel;
}

extern "C" SEXP create_spectral_component(const riem_dcx_t& dcx,
                                          int optimal_iter,
                                          rdcx_filter_type_t filter_type) {
    const int n_fields = 2;
    SEXP spectral = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    if (!dcx.spectral_cache.is_valid) {
        Rf_warning("Spectral cache not valid");
        UNPROTECT(2);
        return R_NilValue;
    }

    const size_t n_eigen = dcx.spectral_cache.eigenvalues.size();
    const Eigen::Index n_verts = dcx.spectral_cache.eigenvectors.rows();

    // Get the optimal eta from GCV history
    double eta_opt = dcx.gcv_history.iterations[optimal_iter].eta_optimal;

    // 1. filtered.eigenvalues: F_η(λ)
    SET_STRING_ELT(names, idx, Rf_mkChar("filtered.eigenvalues"));
    SEXP s_f_lambda = PROTECT(Rf_allocVector(REALSXP, n_eigen));

    // Apply filter based on type
    for (size_t i = 0; i < n_eigen; ++i) {
        double lambda = dcx.spectral_cache.eigenvalues[i];
        double filtered_value;

        switch (filter_type) {
            case rdcx_filter_type_t::HEAT_KERNEL:
                filtered_value = std::exp(-eta_opt * lambda);
                break;
            case rdcx_filter_type_t::TIKHONOV:
                filtered_value = 1.0 / (1.0 + eta_opt * lambda);
                break;
            case rdcx_filter_type_t::CUBIC_SPLINE:
                filtered_value = 1.0 / (1.0 + eta_opt * lambda * lambda);
                break;
            case rdcx_filter_type_t::GAUSSIAN:
                filtered_value = std::exp(-eta_opt * lambda * lambda);
                break;
            case rdcx_filter_type_t::EXPONENTIAL:
                filtered_value = std::exp(-eta_opt * std::sqrt(std::max(lambda, 0.0)));
                break;
            case rdcx_filter_type_t::BUTTERWORTH: {
                // Default to n=2 for Butterworth
                double ratio = lambda / std::max(eta_opt, 1e-15);
                filtered_value = 1.0 / (1.0 + ratio * ratio * ratio * ratio);
                break;
            }
            default:
                filtered_value = 1.0;
                break;
        }

        REAL(s_f_lambda)[i] = filtered_value;
    }
    SET_VECTOR_ELT(spectral, idx++, s_f_lambda);
    UNPROTECT(1);

    // 2. eigenvectors: V
    SET_STRING_ELT(names, idx, Rf_mkChar("eigenvectors"));
    SEXP s_V = PROTECT(Rf_allocMatrix(REALSXP, n_verts, n_eigen));
    double* V_data = REAL(s_V);
    for (size_t j = 0; j < n_eigen; ++j) {
        for (Eigen::Index i = 0; i < n_verts; ++i) {
            V_data[i + j * n_verts] = dcx.spectral_cache.eigenvectors(i, j);
        }
    }
    SET_VECTOR_ELT(spectral, idx++, s_V);
    UNPROTECT(1);

    Rf_setAttrib(spectral, R_NamesSymbol, names);
    UNPROTECT(2);
    return spectral;
}

extern "C" SEXP create_extremality_component(
    const riem_dcx_t& dcx,
    const vec_t& y_hat
) {
    const Eigen::Index n = y_hat.size();

    // Compute extremality scores using full non-diagonal metric
    const bool use_iterative = false;
    vec_t extremality = dcx.compute_all_extremality_scores_full(y_hat,
                                                                use_iterative);


    // Convert to R vector
    SEXP s_extremality = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_extremality)[i] = extremality[i];
    }

    UNPROTECT(1);
    return s_extremality;
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
 * @param s_threshold_percentile
 * @param s_test_stage
 * @param s_verbose SEXP object (logical) controlling progress reporting during computation
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
    SEXP s_t_scale_factor,
    SEXP s_beta_coefficient_factor,
    SEXP s_n_eigenpairs,
    SEXP s_filter_type,
    SEXP s_epsilon_y,
    SEXP s_epsilon_rho,
    SEXP s_max_iterations,
    SEXP s_max_ratio_threshold,
    SEXP s_threshold_percentile,
    SEXP s_density_alpha,
    SEXP s_density_epsilon,
    SEXP s_test_stage,
    SEXP s_verbose
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

    if (!R_FINITE(gamma_modulation)) {
        Rf_error("gamma.modulation must be a finite number (got %.3f)",
                 gamma_modulation);
    }

    // -------------------- scale.factor --------------------
    if (TYPEOF(s_t_scale_factor) != REALSXP || Rf_length(s_t_scale_factor) != 1) {
        Rf_error("t.scale.factor must be a single numeric");
    }

    const double t_scale_factor = REAL(s_t_scale_factor)[0];

    if (!R_FINITE(t_scale_factor)) {
        Rf_error("t.scale.factor must be a finite number (got %.3f)",
                 t_scale_factor);
    }

    // -------------------- beta.coefficient.factor --------------------
    if (TYPEOF(s_beta_coefficient_factor) != REALSXP || Rf_length(s_beta_coefficient_factor) != 1) {
        Rf_error("beta.coefficient.factor must be a single numeric");
    }

    const double beta_coefficient_factor = REAL(s_beta_coefficient_factor)[0];

    if (!R_FINITE(beta_coefficient_factor)) {
        Rf_error("beta.coefficient.factor must be a finite number (got %.3f)",
                 beta_coefficient_factor);
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
    } else if (strcmp(filter_str, "gaussian") == 0) {
        filter_type = rdcx_filter_type_t::GAUSSIAN;
    } else if (strcmp(filter_str, "exponential") == 0) {
        filter_type = rdcx_filter_type_t::EXPONENTIAL;
    } else if (strcmp(filter_str, "butterworth") == 0) {
        filter_type = rdcx_filter_type_t::BUTTERWORTH;
    } else {
        Rf_error("filter.type must be 'heat_kernel', 'tikhonov', 'cubic_spline', "
                 "'gaussian', 'exponential', or 'butterworth' (got '%s')",
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

    // -------------------- density.alpha --------------------

    if (TYPEOF(s_density_alpha) != REALSXP || Rf_length(s_density_alpha) != 1) {
        Rf_error("density.alpha must be a single numeric value");
    }

    const double density_alpha = REAL(s_density_alpha)[0];

    if (!R_FINITE(density_alpha) || density_alpha < 1.0 || density_alpha > 2.0) {
        Rf_error("density.alpha must be finite and in [1, 2] (got %.3f)", density_alpha);
    }

    // -------------------- density.epsilon --------------------

    if (TYPEOF(s_density_epsilon) != REALSXP || Rf_length(s_density_epsilon) != 1) {
        Rf_error("density.epsilon must be a single numeric value");
    }

    const double density_epsilon = REAL(s_density_epsilon)[0];

    if (!R_FINITE(density_epsilon) || density_epsilon <= 0.0) {
        Rf_error("density.epsilon must be a finite positive number (got %.3e)",
                 density_epsilon);
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

    // -------------------- s_verbose --------------------
    if (TYPEOF(s_verbose) != LGLSXP ||
        Rf_length(s_verbose) != 1) {
        Rf_error("verbose must be a single logical value");
    }

    const int verbose_int = LOGICAL(s_verbose)[0];

    if (verbose_int == NA_LOGICAL) {
        Rf_error("verbose cannot be NA");
    }

    const bool verbose = (verbose_int != 0);


    // ================================================================
    // PART II: CALL MEMBER FUNCTION
    // ================================================================

    riem_dcx_t dcx;  // Stack allocation now! No need for new/delete

    try {
        dcx.fit_knn_riem_graph_regression(
            X_sparse,
            y,
            k,
            use_counting_measure,
            density_normalization,
            t_diffusion,
            beta_damping,
            gamma_modulation,
            t_scale_factor,
            beta_coefficient_factor,
            n_eigenpairs,
            filter_type,
            epsilon_y,
            epsilon_rho,
            max_iterations,
            max_ratio_threshold,
            threshold_percentile,
            density_alpha,
            density_epsilon,
            test_stage,
            verbose
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

    // ---------- Select fitted values based on minimum GCV ----------
    // Instead of simply using the last iteration's fitted values, we
    // select the iteration that minimizes the GCV score across all
    // iterations. This implements proper model selection and can
    // prevent overfitting when the iterative refinement continues
    // beyond the optimal point.

    bool has_fitted_values = !dcx.sig.y_hat_hist.empty();
    bool has_gcv_history = !dcx.gcv_history.iterations.empty();
    vec_t y_hat_final;
    int optimal_iteration = -1;

    if (has_fitted_values && has_gcv_history) {
        // Verify consistency between fitted values history and GCV history
        if (dcx.sig.y_hat_hist.size() != dcx.gcv_history.iterations.size()) {
            Rf_warning("Inconsistent history sizes: y_hat_hist has %d entries, "
                       "gcv_history has %d entries. Using last iteration.",
                       (int)dcx.sig.y_hat_hist.size(),
                       (int)dcx.gcv_history.iterations.size());
            y_hat_final = dcx.sig.y_hat_hist.back();
            optimal_iteration = (int)dcx.sig.y_hat_hist.size() - 1;
        } else {
            // Find iteration with minimum GCV score
            double min_gcv = std::numeric_limits<double>::max();
            size_t min_idx = 0;

            for (size_t i = 0; i < dcx.gcv_history.iterations.size(); ++i) {
                double gcv_score = dcx.gcv_history.iterations[i].gcv_optimal;
                if (gcv_score < min_gcv) {
                    min_gcv = gcv_score;
                    min_idx = i;
                }
            }

            // Use fitted values from the iteration with minimum GCV
            y_hat_final = dcx.sig.y_hat_hist[min_idx];
            optimal_iteration = (int)min_idx;

            // Inform user if optimal iteration differs from final iteration
            if (test_stage < 0 && min_idx != dcx.sig.y_hat_hist.size() - 1) {
                Rprintf("Selected iteration %d (GCV = %.6e) over final iteration %d (GCV = %.6e)\n",
                        (int)min_idx, min_gcv,
                        (int)(dcx.sig.y_hat_hist.size() - 1),
                        dcx.gcv_history.iterations.back().gcv_optimal);
            }
        }
    } else if (has_fitted_values) {
        // GCV history unavailable (should not happen in normal operation)
        // Fall back to last iteration
        Rf_warning("GCV history unavailable. Using last iteration's fitted values.");
        y_hat_final = dcx.sig.y_hat_hist.back();
        optimal_iteration = (int)dcx.sig.y_hat_hist.size() - 1;
    } else {
        // Early termination: use observed y as placeholder
        // This allows the result structure to be consistent
        y_hat_final = y;
        optimal_iteration = -1;

        // Warn user if this wasn't intentional
        if (test_stage < 0) {
            Rf_warning("No fitted values computed (unexpected early termination). "
                       "Returning observed y values as placeholder.");
        }
    }

    // ---------- Main result list ----------

    const int n_components = 12;
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

    // Component 3: optimal.iteration
    SET_STRING_ELT(names, component_idx, Rf_mkChar("optimal.iteration"));
    SET_VECTOR_ELT(result, component_idx++, Rf_ScalarInteger(optimal_iteration + 1));  // R uses 1-based indexing

    // Component 4: graph (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("graph"));
    SEXP s_graph = PROTECT(create_graph_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_graph);
    UNPROTECT(1);

    // Component 5: iteration (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("iteration"));
    SEXP s_iteration = PROTECT(create_iteration_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_iteration);
    UNPROTECT(1);

    // Component 6: parameters (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("parameters"));
    SEXP s_params = PROTECT(create_parameters_component(
                                k, use_counting_measure, density_normalization,
                                t_diffusion, beta_damping, gamma_modulation,
                                t_scale_factor, beta_coefficient_factor,
                                n_eigenpairs, filter_type, epsilon_y, epsilon_rho, max_iterations,
                                density_alpha, density_epsilon
                                ));
    SET_VECTOR_ELT(result, component_idx++, s_params);
    UNPROTECT(1);

    // Component 7: y (original response)
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

    // Component 8: gcv (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("gcv"));
    SEXP s_gcv = PROTECT(create_gcv_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_gcv);
    UNPROTECT(1);

    // Component 9: density
    SET_STRING_ELT(names, component_idx, Rf_mkChar("density"));
    SEXP s_density = PROTECT(create_density_history_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_density);
    UNPROTECT(1);

    // Component 10: gamma.selection
    SET_STRING_ELT(names, component_idx, Rf_mkChar("gamma.selection"));
    SEXP s_gamma_sel = PROTECT(create_gamma_selection_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_gamma_sel);
    UNPROTECT(1);

    // Component 11: spectral
    SET_STRING_ELT(names, component_idx, Rf_mkChar("spectral"));
    SEXP s_spectral = PROTECT(create_spectral_component(dcx, optimal_iteration, filter_type));
    SET_VECTOR_ELT(result, component_idx++, s_spectral);
    UNPROTECT(1);

    // Component 12: extremality.scores
    SET_STRING_ELT(names, component_idx, Rf_mkChar("extremality.scores"));
    SEXP s_extremality = PROTECT(create_extremality_component(dcx, y_hat_final));
    SET_VECTOR_ELT(result, component_idx++, s_extremality);
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
