#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Print.h>

// Undefine conflicting macros after including R headers
#undef length

#include "set_wgraph.hpp"
#include "deg0_lowess_graph_smoothing.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "error_utils.h" // for REPORT_ERROR()

extern "C" {
    SEXP S_deg0_lowess_graph_smoothing(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_X,
        SEXP s_max_iterations,
        SEXP s_convergence_threshold,
        SEXP s_convergence_type,
        SEXP s_k,
        SEXP s_pruning_thld,
        SEXP s_n_bws,
        SEXP s_log_grid,
        SEXP s_min_bw_factor,
        SEXP s_max_bw_factor,
        SEXP s_dist_normalization_factor,
        SEXP s_kernel_type,
        SEXP s_n_folds,
        SEXP s_use_uniform_weights,
		SEXP s_outlier_thld,
		SEXP s_with_bw_predictions,
        SEXP s_switch_to_residuals_after,
        SEXP s_verbose
    );
}

/**
 * @brief R interface for the deg0_lowess_graph_smoothing function
 */
SEXP S_deg0_lowess_graph_smoothing(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_X,
    SEXP s_max_iterations,
    SEXP s_convergence_threshold,
    SEXP s_convergence_type,
    SEXP s_k,
    SEXP s_pruning_thld,
    SEXP s_n_bws,
    SEXP s_log_grid,
    SEXP s_min_bw_factor,
    SEXP s_max_bw_factor,
    SEXP s_dist_normalization_factor,
    SEXP s_kernel_type,
    SEXP s_n_folds,
    SEXP s_use_uniform_weights,
	SEXP s_outlier_thld,
    SEXP s_with_bw_predictions,
    SEXP s_switch_to_residuals_after,
    SEXP s_verbose
) {
    int nprot = 0;

    // Extract and validate parameters
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Check and coerce X to matrix if necessary
    PROTECT(s_X = coerceVector(s_X, REALSXP)); nprot++;
    if (!isMatrix(s_X)) {
        UNPROTECT(nprot);
        REPORT_ERROR("X must be a numeric matrix");
    }

    // Get dimensions
    int *X_dims = INTEGER(getAttrib(s_X, R_DimSymbol));
    int n_samples = X_dims[0];
    int n_features = X_dims[1];

    // Extract scalar parameters
    size_t max_iterations = static_cast<size_t>(asInteger(s_max_iterations));
    double convergence_threshold = asReal(s_convergence_threshold);
    int convergence_type_int = asInteger(s_convergence_type);
    ConvergenceCriteriaType convergence_type;

    // Convert convergence type
    switch (convergence_type_int) {
        case 1:
            convergence_type = ConvergenceCriteriaType::MAX_ABSOLUTE_DIFF;
            break;
        case 2:
            convergence_type = ConvergenceCriteriaType::MEAN_ABSOLUTE_DIFF;
            break;
        case 3:
            convergence_type = ConvergenceCriteriaType::RELATIVE_CHANGE;
            break;
        default:
            UNPROTECT(nprot);
            REPORT_ERROR("Invalid convergence type: %d", convergence_type_int);
    }

    size_t k = static_cast<size_t>(asInteger(s_k));
    double pruning_thld = asReal(s_pruning_thld);
    size_t n_bws = static_cast<size_t>(asInteger(s_n_bws));
    bool log_grid = asLogical(s_log_grid);
    double min_bw_factor = asReal(s_min_bw_factor);
    double max_bw_factor = asReal(s_max_bw_factor);
    double dist_normalization_factor = asReal(s_dist_normalization_factor);
    size_t kernel_type = static_cast<size_t>(asInteger(s_kernel_type));
    size_t n_folds = static_cast<size_t>(asInteger(s_n_folds));
    bool use_uniform_weights = asLogical(s_use_uniform_weights);
	double outlier_thld = asReal(s_outlier_thld);
    bool with_bw_predictions = asLogical(s_with_bw_predictions);
    size_t switch_to_residuals_after = static_cast<size_t>(asInteger(s_switch_to_residuals_after));
    bool verbose = asLogical(s_verbose);

    // Construct initial graph
    set_wgraph_t initial_graph(adj_list, weight_list);

    // Find diameter endpoints
    auto [end1, diam] = initial_graph.get_vertex_eccentricity(0);  // Start from vertex 0
    auto [end2, diameter] = initial_graph.get_vertex_eccentricity(end1);
    initial_graph.graph_diameter = diameter;

    // Convert R matrix to C++ vectors (column-major to row-major)
    std::vector<std::vector<double>> X(n_features, std::vector<double>(n_samples));
    double *X_data = REAL(s_X);

    for (int j = 0; j < n_features; ++j) {
        for (int i = 0; i < n_samples; ++i) {
            X[j][i] = X_data[i + j * n_samples];
        }
    }

    // Call the C++ function
    deg0_lowess_graph_smoothing_t result;
    try {
        result = deg0_lowess_graph_smoothing(
            initial_graph,
            X,
            max_iterations,
            convergence_threshold,
            convergence_type,
            k,
            pruning_thld,
            n_bws,
            log_grid,
            min_bw_factor,
            max_bw_factor,
            dist_normalization_factor,
            kernel_type,
            n_folds,
            use_uniform_weights,
			outlier_thld,
            with_bw_predictions,
            switch_to_residuals_after,
            verbose
        );
    } catch (const std::exception &e) {
        UNPROTECT(nprot);
        REPORT_ERROR("Error in deg0_lowess_graph_smoothing: %s", e.what());
    }

    // Create R objects for the results
    SEXP s_result = PROTECT(allocVector(VECSXP, 7)); nprot++;  // 7 components in the result

    // 1. Smoothed graphs
    int n_iterations = result.iterations_performed;
    SEXP s_smoothed_graphs = PROTECT(allocVector(VECSXP, n_iterations)); nprot++;

    for (int iter = 0; iter < n_iterations; ++iter) {
        const set_wgraph_t& graph = result.smoothed_graphs[iter];
        size_t n_vertices = graph.adjacency_list.size();

        // Create lists for adjacency and weights
        SEXP s_graph_adj_list = PROTECT(allocVector(VECSXP, n_vertices));
        SEXP s_graph_weight_list = PROTECT(allocVector(VECSXP, n_vertices));

        for (size_t i = 0; i < n_vertices; ++i) {
            // Get neighbors and weights
            const auto& edges = graph.adjacency_list[i];

            // Create vectors for this vertex
            SEXP s_neighbors = PROTECT(allocVector(INTSXP, edges.size()));
            SEXP s_weights = PROTECT(allocVector(REALSXP, edges.size()));

            int* neighbors = INTEGER(s_neighbors);
            double* weights = REAL(s_weights);

            // Fill vectors
            size_t j = 0;
            for (const auto& edge : edges) {
                neighbors[j] = edge.vertex + 1;  // Convert to 1-based indexing
                weights[j] = edge.weight;
                j++;
            }

            // Set vectors in lists
            SET_VECTOR_ELT(s_graph_adj_list, i, s_neighbors);
            SET_VECTOR_ELT(s_graph_weight_list, i, s_weights);

            UNPROTECT(2);  // s_neighbors, s_weights
        }

        // Create graph list
        SEXP s_graph_list = PROTECT(allocVector(VECSXP, 2));
        SET_VECTOR_ELT(s_graph_list, 0, s_graph_adj_list);
        SET_VECTOR_ELT(s_graph_list, 1, s_graph_weight_list);

        // Set names
        SEXP s_graph_names = PROTECT(allocVector(STRSXP, 2));
        SET_STRING_ELT(s_graph_names, 0, mkChar("adj_list"));
        SET_STRING_ELT(s_graph_names, 1, mkChar("weight_list"));
        setAttrib(s_graph_list, R_NamesSymbol, s_graph_names);

        // Set in list of graphs
        SET_VECTOR_ELT(s_smoothed_graphs, iter, s_graph_list);

        UNPROTECT(4);  // s_graph_adj_list, s_graph_weight_list, s_graph_list, s_graph_names
    }

    // 2. Smoothed X matrices
    SEXP s_smoothed_X = PROTECT(allocVector(VECSXP, n_iterations + 1)); nprot++;  // +1 because we have one more X than graphs

    for (int iter = 0; iter <= n_iterations; ++iter) {
        const std::vector<std::vector<double>>& X_iter = result.smoothed_X[iter];

        // Create matrix
        SEXP s_X_matrix = PROTECT(allocMatrix(REALSXP, n_samples, n_features));
        double* X_data = REAL(s_X_matrix);

        // Fill matrix (convert from row-major to column-major)
        for (int j = 0; j < n_features; ++j) {
            for (int i = 0; i < n_samples; ++i) {
                X_data[i + j * n_samples] = X_iter[j][i];
            }
        }

        // Set in list
        SET_VECTOR_ELT(s_smoothed_X, iter, s_X_matrix);

        UNPROTECT(1);  // s_X_matrix
    }

    // 3. Convergence metrics
    SEXP s_convergence_metrics = PROTECT(allocVector(REALSXP, result.convergence_metrics.size())); nprot++;
    double* metrics_data = REAL(s_convergence_metrics);

    for (size_t i = 0; i < result.convergence_metrics.size(); ++i) {
        metrics_data[i] = result.convergence_metrics[i];
    }

    // 4. Iterations performed
    SEXP s_iterations_performed = PROTECT(ScalarInteger(result.iterations_performed)); nprot++;

    // 5. Used boosting flag
    SEXP s_used_boosting = PROTECT(ScalarLogical(result.used_boosting)); nprot++;

    // 6. Vertex mappings (list of integer vectors)
    SEXP s_vertex_mappings = PROTECT(allocVector(VECSXP, result.vertex_mappings.size())); nprot++;

    for (size_t i = 0; i < result.vertex_mappings.size(); ++i) {
        const std::vector<size_t>& mapping = result.vertex_mappings[i];

        // Create integer vector for this mapping
        SEXP s_mapping = PROTECT(allocVector(INTSXP, mapping.size()));
        int* mapping_data = INTEGER(s_mapping);

        // Fill vector with 1-based indexing for R
        for (size_t j = 0; j < mapping.size(); ++j) {
            mapping_data[j] = mapping[j] + 1;  // Convert to 1-based indexing for R
        }

        // Set in list
        SET_VECTOR_ELT(s_vertex_mappings, i, s_mapping);

        UNPROTECT(1);  // s_mapping
    }

    // 7. Outlier indices (integer vector)
    SEXP s_outlier_indices = PROTECT(allocVector(INTSXP, result.outlier_indices.size())); nprot++;
    int* outlier_data = INTEGER(s_outlier_indices);

    // Fill vector with 1-based indexing for R
    for (size_t i = 0; i < result.outlier_indices.size(); ++i) {
        outlier_data[i] = result.outlier_indices[i] + 1;  // Convert to 1-based indexing for R
    }

    // Set results in output list
    SET_VECTOR_ELT(s_result, 0, s_smoothed_graphs);
    SET_VECTOR_ELT(s_result, 1, s_smoothed_X);
    SET_VECTOR_ELT(s_result, 2, s_convergence_metrics);
    SET_VECTOR_ELT(s_result, 3, s_iterations_performed);
    SET_VECTOR_ELT(s_result, 4, s_used_boosting);
    SET_VECTOR_ELT(s_result, 5, s_vertex_mappings);
    SET_VECTOR_ELT(s_result, 6, s_outlier_indices);

    // Set names
    SEXP s_result_names = PROTECT(allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(s_result_names, 0, mkChar("smoothed_graphs"));
    SET_STRING_ELT(s_result_names, 1, mkChar("smoothed_X"));
    SET_STRING_ELT(s_result_names, 2, mkChar("convergence_metrics"));
    SET_STRING_ELT(s_result_names, 3, mkChar("iterations_performed"));
    SET_STRING_ELT(s_result_names, 4, mkChar("used_boosting"));
    SET_STRING_ELT(s_result_names, 5, mkChar("vertex_mappings"));
    SET_STRING_ELT(s_result_names, 6, mkChar("outlier_indices"));
    setAttrib(s_result, R_NamesSymbol, s_result_names);

    UNPROTECT(nprot);
    return s_result;
}
