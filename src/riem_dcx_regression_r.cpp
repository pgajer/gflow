#include "riem_dcx.hpp"
#include <R.h>
#include <Rinternals.h>

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
    SEXP s_max_iterations
    ) {
    // ==================== Input Extraction ====================

    // Convert X to Eigen sparse matrix
    Eigen::SparseMatrix<double> X_sparse;
    // ... (same logic as S_build_nerve_from_knn)

    // Extract y
    const int n = Rf_length(s_y);
    vec_t y(n);
    // ... copy data

    // Extract k
    const int k = Rf_asInteger(s_k);

    // Extract parameters from list
    // ... extract all parameters

    // ==================== Call Member Function ====================

    riem_dcx_t* dcx = new riem_dcx_t();

    try {
        dcx->fit_knn_riem_graph_regression(
            X_sparse, y, k,
            use_counting_measure, density_normalization,
            t_diffusion, beta_damping, gamma_modulation,
            n_eigenpairs, filter_type,
            epsilon_y, epsilon_rho, max_iterations
        );
    } catch (const std::exception& e) {
        delete dcx;
        Rf_error("Regression fitting failed: %s", e.what());
    }

    // ==================== Create External Pointer ====================

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

    // Set class
    SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(class_attr, 0, Rf_mkChar("Rcpp_riem_dcx"));
    Rf_setAttrib(ext_ptr, R_ClassSymbol, class_attr);

    UNPROTECT(2);
    return ext_ptr;
}
