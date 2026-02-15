#include "set_wgraph.hpp"               // set_wgraph_t
#include "graph_spectral_filter.hpp"    // graph_spectral_filter_t
#include "SEXP_cpp_conversion_utils.hpp"// convert_adj_list_from_R, convert_weight_list_from_R, Rvect_to_CppVect_double

#include <R.h>
#include <Rinternals.h>

extern "C" SEXP S_graph_spectral_filter(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_laplacian_type,
    SEXP s_filter_type,
    SEXP s_laplacian_power,
    SEXP s_kernel_tau_factor,
    SEXP s_kernel_radius_factor,
    SEXP s_kernel_type,
    SEXP s_kernel_adaptive,
    SEXP s_min_radius_factor,
    SEXP s_max_radius_factor,
    SEXP s_domain_min_size,
    SEXP s_precision,
    SEXP s_n_evectors_to_compute,
    SEXP s_n_candidates,
    SEXP s_log_grid,
    SEXP s_with_t_predictions,
    SEXP s_verbose
) {
    // 1) Convert inputs (LENGTH-first policy; avoid long vectors)
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    const R_len_t ny = Rf_length(s_y);
    std::vector<double> y(REAL(s_y), REAL(s_y) + (ny > 0 ? ny : 0));

    // Convert enum/scalar parameters
    laplacian_type_t laplacian_type = static_cast<laplacian_type_t>(Rf_asInteger(s_laplacian_type));
    filter_type_t    filter_type    = static_cast<filter_type_t>(Rf_asInteger(s_filter_type));
    size_t           laplacian_power = static_cast<size_t>(Rf_asInteger(s_laplacian_power));

    // Kernel params
    kernel_params_t kernel_params;
    kernel_params.tau_factor        = Rf_asReal(s_kernel_tau_factor);
    kernel_params.radius_factor     = Rf_asReal(s_kernel_radius_factor);
    kernel_params.kernel_type       = static_cast<kernel_type_t>(Rf_asInteger(s_kernel_type));
    kernel_params.adaptive          = (Rf_asLogical(s_kernel_adaptive) == TRUE);
    kernel_params.min_radius_factor = Rf_asReal(s_min_radius_factor);
    kernel_params.max_radius_factor = Rf_asReal(s_max_radius_factor);
    kernel_params.domain_min_size   = static_cast<size_t>(Rf_asInteger(s_domain_min_size));
    kernel_params.precision         = Rf_asReal(s_precision);

    // Remaining parameters
    size_t n_evectors_to_compute = static_cast<size_t>(Rf_asInteger(s_n_evectors_to_compute));
    size_t n_candidates          = static_cast<size_t>(Rf_asInteger(s_n_candidates));
    const bool log_grid          = (Rf_asLogical(s_log_grid) == TRUE);
    const bool with_t_predictions= (Rf_asLogical(s_with_t_predictions) == TRUE);
    const bool verbose           = (Rf_asLogical(s_verbose) == TRUE);

    // 2) Call the C++ graph filter
    set_wgraph_t graph(adj_list, weight_list);
    graph.compute_graph_diameter();

    graph_spectral_filter_t result = graph.graph_spectral_filter(
        y,
        laplacian_type,
        filter_type,
        laplacian_power,
        kernel_params,
        n_evectors_to_compute,
        n_candidates,
        log_grid,
        with_t_predictions,
        verbose
    );

    // 3) Build the R list to return
    const int n_fields = with_t_predictions ? 13 : 12;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_fields));                 // [1]

    // names
    {
        SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));                // [2]
        SET_STRING_ELT(names, 0,  Rf_mkChar("evalues"));
        SET_STRING_ELT(names, 1,  Rf_mkChar("evectors"));
        SET_STRING_ELT(names, 2,  Rf_mkChar("candidate_ts"));
        SET_STRING_ELT(names, 3,  Rf_mkChar("gcv_scores"));
        SET_STRING_ELT(names, 4,  Rf_mkChar("opt_t_idx"));
        SET_STRING_ELT(names, 5,  Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 6,  Rf_mkChar("t_predictions"));
        SET_STRING_ELT(names, 7,  Rf_mkChar("laplacian_type"));
        SET_STRING_ELT(names, 8,  Rf_mkChar("filter_type"));
        SET_STRING_ELT(names, 9,  Rf_mkChar("laplacian_power"));
        SET_STRING_ELT(names, 10, Rf_mkChar("kernel_params"));
        SET_STRING_ELT(names, 11, Rf_mkChar("compute_time_ms"));
        if (with_t_predictions) {
            SET_STRING_ELT(names, 12, Rf_mkChar("gcv_min_score"));
        }
        Rf_setAttrib(r_result, R_NamesSymbol, names);
        UNPROTECT(1); // names                                                  // [-]
    }

    // evalues
    {
        const int n = (int) result.evalues.size();
        SEXP Revalues = PROTECT(Rf_allocVector(REALSXP, n));                    // [3]
        double* p = REAL(Revalues);
        for (int i = 0; i < n; ++i) p[i] = result.evalues[i];
        SET_VECTOR_ELT(r_result, 0, Revalues);
        UNPROTECT(1); // Revalues                                               // [-]
    }

    // evectors (column-major)
    {
        const int nrow = (int) result.evectors.rows();
        const int ncol = (int) result.evectors.cols();
        SEXP Revectors = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol));          // [4]
        double* p = REAL(Revectors);
        for (int c = 0; c < ncol; ++c)
            for (int r = 0; r < nrow; ++r)
                p[r + nrow * c] = result.evectors(r, c);
        SET_VECTOR_ELT(r_result, 1, Revectors);
        UNPROTECT(1); // Revectors                                              // [-]
    }

    // candidate_ts
    {
        const int n = (int) result.candidate_ts.size();
        SEXP Rts = PROTECT(Rf_allocVector(REALSXP, n));                         // [5]
        double* p = REAL(Rts);
        for (int i = 0; i < n; ++i) p[i] = result.candidate_ts[i];
        SET_VECTOR_ELT(r_result, 2, Rts);
        UNPROTECT(1); // Rts                                                    // [-]
    }

    // gcv_scores
    {
        const int n = (int) result.gcv_scores.size();
        SEXP Rgcv = PROTECT(Rf_allocVector(REALSXP, n));                        // [6]
        double* p = REAL(Rgcv);
        for (int i = 0; i < n; ++i) p[i] = result.gcv_scores[i];
        SET_VECTOR_ELT(r_result, 3, Rgcv);
        UNPROTECT(1); // Rgcv                                                   // [-]
    }

    // opt_t_idx (1-based for R)
    {
        SEXP Ropt = PROTECT(Rf_allocVector(INTSXP, 1));                         // [7]
        INTEGER(Ropt)[0] = (int) (result.opt_t_idx + 1);
        SET_VECTOR_ELT(r_result, 4, Ropt);
        UNPROTECT(1); // Ropt                                                   // [-]
    }

    // predictions
    {
        const int n = (int) result.predictions.size();
        SEXP Rpred = PROTECT(Rf_allocVector(REALSXP, n));                       // [8]
        double* p = REAL(Rpred);
        for (int i = 0; i < n; ++i) p[i] = result.predictions[i];
        SET_VECTOR_ELT(r_result, 5, Rpred);
        UNPROTECT(1); // Rpred                                                  // [-]
    }

    // t_predictions (optional)
    {
        if (with_t_predictions && !result.t_predictions.empty()) {
            const int nvert = (int) result.t_predictions[0].size();
            const int nt    = (int) result.t_predictions.size();
            SEXP Rtp = PROTECT(Rf_allocMatrix(REALSXP, nvert, nt));             // [9]
            double* p = REAL(Rtp);
            for (int j = 0; j < nt; ++j) {
                const std::vector<double>& col = result.t_predictions[j];
                for (int i = 0; i < nvert; ++i) {
                    p[i + nvert * j] = col[i];
                }
            }
            SET_VECTOR_ELT(r_result, 6, Rtp);
            UNPROTECT(1); // Rtp                                                // [-]
        } else {
            SET_VECTOR_ELT(r_result, 6, R_NilValue);
        }
    }

    // laplacian_type
    {
        SEXP Rlap_type = PROTECT(Rf_allocVector(INTSXP, 1));                    // [10]
        INTEGER(Rlap_type)[0] = (int) result.laplacian_type;
        SET_VECTOR_ELT(r_result, 7, Rlap_type);
        UNPROTECT(1); // Rlap_type                                              // [-]
    }

    // filter_type
    {
        SEXP Rfilter_type = PROTECT(Rf_allocVector(INTSXP, 1));                 // [11]
        INTEGER(Rfilter_type)[0] = (int) result.filter_type;
        SET_VECTOR_ELT(r_result, 8, Rfilter_type);
        UNPROTECT(1); // Rfilter_type                                           // [-]
    }

    // laplacian_power
    {
        SEXP Rlap_power = PROTECT(Rf_allocVector(INTSXP, 1));                   // [12]
        INTEGER(Rlap_power)[0] = (int) result.laplacian_power;
        SET_VECTOR_ELT(r_result, 9, Rlap_power);
        UNPROTECT(1); // Rlap_power                                             // [-]
    }

    // kernel_params (list)
    {
        SEXP Rkernel_params = PROTECT(Rf_allocVector(VECSXP, 8));               // [13]

        // 0: tau_factor
        {
            SEXP x = PROTECT(Rf_allocVector(REALSXP, 1));                       // [14]
            REAL(x)[0] = result.kernel_params.tau_factor;
            SET_VECTOR_ELT(Rkernel_params, 0, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 1: radius_factor
        {
            SEXP x = PROTECT(Rf_allocVector(REALSXP, 1));                       // [15]
            REAL(x)[0] = result.kernel_params.radius_factor;
            SET_VECTOR_ELT(Rkernel_params, 1, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 2: kernel_type
        {
            SEXP x = PROTECT(Rf_allocVector(INTSXP, 1));                        // [16]
            INTEGER(x)[0] = (int) result.kernel_params.kernel_type;
            SET_VECTOR_ELT(Rkernel_params, 2, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 3: adaptive
        {
            SEXP x = PROTECT(Rf_allocVector(LGLSXP, 1));                        // [17]
            LOGICAL(x)[0] = result.kernel_params.adaptive ? 1 : 0;
            SET_VECTOR_ELT(Rkernel_params, 3, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 4: min_radius_factor
        {
            SEXP x = PROTECT(Rf_allocVector(REALSXP, 1));                       // [18]
            REAL(x)[0] = result.kernel_params.min_radius_factor;
            SET_VECTOR_ELT(Rkernel_params, 4, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 5: max_radius_factor
        {
            SEXP x = PROTECT(Rf_allocVector(REALSXP, 1));                       // [19]
            REAL(x)[0] = result.kernel_params.max_radius_factor;
            SET_VECTOR_ELT(Rkernel_params, 5, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 6: domain_min_size
        {
            SEXP x = PROTECT(Rf_allocVector(INTSXP, 1));                        // [20]
            INTEGER(x)[0] = (int) result.kernel_params.domain_min_size;
            SET_VECTOR_ELT(Rkernel_params, 6, x);
            UNPROTECT(1); // x                                                  // [-]
        }
        // 7: precision
        {
            SEXP x = PROTECT(Rf_allocVector(REALSXP, 1));                       // [21]
            REAL(x)[0] = result.kernel_params.precision;
            SET_VECTOR_ELT(Rkernel_params, 7, x);
            UNPROTECT(1); // x                                                  // [-]
        }

        // Set names
        {
            SEXP nms = PROTECT(Rf_allocVector(STRSXP, 8));                      // [22]
            SET_STRING_ELT(nms, 0, Rf_mkChar("tau_factor"));
            SET_STRING_ELT(nms, 1, Rf_mkChar("radius_factor"));
            SET_STRING_ELT(nms, 2, Rf_mkChar("kernel_type"));
            SET_STRING_ELT(nms, 3, Rf_mkChar("adaptive"));
            SET_STRING_ELT(nms, 4, Rf_mkChar("min_radius_factor"));
            SET_STRING_ELT(nms, 5, Rf_mkChar("max_radius_factor"));
            SET_STRING_ELT(nms, 6, Rf_mkChar("domain_min_size"));
            SET_STRING_ELT(nms, 7, Rf_mkChar("precision"));
            Rf_setAttrib(Rkernel_params, R_NamesSymbol, nms);
            UNPROTECT(1); // nms                                                // [-]
        }

        SET_VECTOR_ELT(r_result, 10, Rkernel_params);
        UNPROTECT(1); // Rkernel_params                                         // [-]
    }

    // compute_time_ms
    {
        SEXP Rcompute_time = PROTECT(Rf_allocVector(REALSXP, 1));               // [23]
        REAL(Rcompute_time)[0] = result.compute_time_ms;
        SET_VECTOR_ELT(r_result, 11, Rcompute_time);
        UNPROTECT(1); // Rcompute_time                                          // [-]
    }

    // gcv_min_score (if with_t_predictions)
    if (with_t_predictions) {
        SEXP Rgcv_min = PROTECT(Rf_allocVector(REALSXP, 1));                    // [24]
        REAL(Rgcv_min)[0] = result.gcv_min_score;
        SET_VECTOR_ELT(r_result, 12, Rgcv_min);
        UNPROTECT(1); // Rgcv_min                                               // [-]
    }

    UNPROTECT(1); // r_result                                                   // [-]
    return r_result;
}
