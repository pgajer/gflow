// centered_paths_hardened.cpp — rchk-hardened drop-in for S_ugg_get_path_data
//
// Hardenings vs provided version:
// - Defensive scalar extraction via Rf_as* + NA/range checks.
// - Coercion block for coordinates with PROTECT_WITH_INDEX/REPROTECT; copy to STL then UNPROTECT immediately.
// - Validate coordinates is a numeric matrix with n_points>0, n_dims>0.
// - Bounds checks for start/end vertices after 0-based conversion.
// - Result assembly: only `result` + `names` protected until tail (UNPROTECT(2)); all others are local protects.
//
// Assumptions: helper converters and core algorithm are available with same signatures.

#include <vector>
#include <stdexcept>
#include <R.h>
#include <Rinternals.h>

// Helpers (provided elsewhere)
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);
extern SEXP convert_vector_double_to_R(const std::vector<double>&);
extern SEXP convert_vector_int_to_R(const std::vector<int>&);
extern SEXP convert_matrix_double_to_R(const std::vector<std::vector<double>>&);

// Core structures
struct path_data_t {
    std::vector<std::vector<double>> path_coordinates;
    std::vector<int>    path_indices;
    std::vector<double> path_distances;
    std::vector<double> path_curvatures;
    double total_length;
    int    n_points;
    bool   is_closed;
};

// Core algorithm (provided elsewhere)
path_data_t ugg_get_path_data(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<std::vector<double>>& coordinates,
    int start_vertex,
    int end_vertex,
    double smoothing_bandwidth,
    bool compute_curvature,
    bool verbose);

extern "C" {

SEXP S_ugg_get_path_data(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP coordinates_sexp,
    SEXP start_vertex_sexp,
    SEXP end_vertex_sexp,
    SEXP smoothing_bandwidth_sexp,
    SEXP compute_curvature_sexp,
    SEXP verbose_sexp
) {
    // --- Coercion block for coordinates ---
    std::vector<std::vector<double>> coordinates;
    int n_points = 0;
    int n_dims   = 0;
    {
        SEXP coords = coordinates_sexp;
        PROTECT_INDEX pcoords;
        PROTECT_WITH_INDEX(coords, &pcoords);
        if (TYPEOF(coords) != REALSXP) {
            REPROTECT(coords = Rf_coerceVector(coords, REALSXP), pcoords);
        }

        // Expect a matrix
        SEXP dims = Rf_getAttrib(coords, R_DimSymbol);
        if (Rf_isNull(dims) || LENGTH(dims) != 2) {
            UNPROTECT(1);
            Rf_error("S_ugg_get_path_data(): 'coordinates' must be a numeric matrix.");
        }
        if (TYPEOF(dims) != INTSXP) {
            UNPROTECT(1);
            Rf_error("S_ugg_get_path_data(): 'coordinates' dim attribute must be integer.");
        }
        int* dimptr = INTEGER(dims);
        n_points = dimptr[0];
        n_dims   = dimptr[1];
        if (n_points <= 0 || n_dims <= 0) {
            UNPROTECT(1);
            Rf_error("S_ugg_get_path_data(): 'coordinates' must have positive dimensions.");
        }

        // Copy to STL
        coordinates.assign(n_points, std::vector<double>(n_dims));
        double* cptr = REAL(coords);
        for (int i = 0; i < n_points; ++i) {
            for (int j = 0; j < n_dims; ++j) {
                coordinates[(size_t)i][(size_t)j] = cptr[i + j * (R_xlen_t)n_points];
            }
        }
        UNPROTECT(1); // coords
    }

    // --- Container-first for graph inputs ---
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);

    // --- Defensive scalars ---
    int    start_vertex_i = Rf_asInteger(start_vertex_sexp);
    int    end_vertex_i   = Rf_asInteger(end_vertex_sexp);
    double smoothing_bw   = Rf_asReal(smoothing_bandwidth_sexp);
    int    compute_curv_i = Rf_asLogical(compute_curvature_sexp);
    int    verbose_i      = Rf_asLogical(verbose_sexp);

    // --- NA / range checks ---
    if (start_vertex_i == NA_INTEGER || end_vertex_i == NA_INTEGER) {
        Rf_error("S_ugg_get_path_data(): start_vertex and end_vertex cannot be NA.");
    }
    if (ISNAN(smoothing_bw) || smoothing_bw < 0.0) {
        Rf_error("S_ugg_get_path_data(): smoothing_bandwidth must be >= 0.");
    }
    if (compute_curv_i == NA_LOGICAL) {
        Rf_error("S_ugg_get_path_data(): compute_curvature must be TRUE/FALSE.");
    }
    if (verbose_i == NA_LOGICAL) {
        Rf_error("S_ugg_get_path_data(): verbose must be TRUE/FALSE.");
    }

    // Convert to 0-based; validate bounds vs n_points
    const int start_vertex = start_vertex_i - 1;
    const int end_vertex   = end_vertex_i   - 1;
    if (start_vertex < 0 || start_vertex >= n_points) {
        Rf_error("S_ugg_get_path_data(): start_vertex out of bounds (1..%d).", n_points);
    }
    if (end_vertex < 0 || end_vertex >= n_points) {
        Rf_error("S_ugg_get_path_data(): end_vertex out of bounds (1..%d).", n_points);
    }

    const bool compute_curvature = (compute_curv_i == TRUE);
    const bool verbose           = (verbose_i == TRUE);

    // --- Core computation ---
    path_data_t path_data = ugg_get_path_data(
        adj_list,
        weight_list,
        coordinates,
        start_vertex,
        end_vertex,
        smoothing_bw,
        compute_curvature,
        verbose
    );

    // --- Result assembly ---
    const int N = 7;
    SEXP result = PROTECT(Rf_allocVector(VECSXP, N));
    SEXP names  = PROTECT(Rf_allocVector(STRSXP, N));
    SET_STRING_ELT(names, 0, Rf_mkChar("path_coordinates"));
    SET_STRING_ELT(names, 1, Rf_mkChar("path_indices"));
    SET_STRING_ELT(names, 2, Rf_mkChar("path_distances"));
    SET_STRING_ELT(names, 3, Rf_mkChar("path_curvatures"));
    SET_STRING_ELT(names, 4, Rf_mkChar("total_length"));
    SET_STRING_ELT(names, 5, Rf_mkChar("n_points"));
    SET_STRING_ELT(names, 6, Rf_mkChar("is_closed"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    // 0: path_coordinates (matrix)
    {
        SEXP m = PROTECT(convert_matrix_double_to_R(path_data.path_coordinates));
        SET_VECTOR_ELT(result, 0, m);
        UNPROTECT(1);
    }
    // 1: path_indices (1-based)
    {
        std::vector<int> idx = path_data.path_indices;
        for (size_t k = 0; k < idx.size(); ++k) idx[k] += 1;
        SEXP v = PROTECT(convert_vector_int_to_R(idx));
        SET_VECTOR_ELT(result, 1, v);
        UNPROTECT(1);
    }
    // 2: path_distances
    {
        SEXP v = PROTECT(convert_vector_double_to_R(path_data.path_distances));
        SET_VECTOR_ELT(result, 2, v);
        UNPROTECT(1);
    }
    // 3: path_curvatures
    {
        SEXP v = PROTECT(convert_vector_double_to_R(path_data.path_curvatures));
        SET_VECTOR_ELT(result, 3, v);
        UNPROTECT(1);
    }
    // 4: total_length
    {
        SEXP s = PROTECT(Rf_ScalarReal(path_data.total_length));
        SET_VECTOR_ELT(result, 4, s);
        UNPROTECT(1);
    }
    // 5: n_points
    {
        SEXP s = PROTECT(Rf_ScalarInteger(path_data.n_points));
        SET_VECTOR_ELT(result, 5, s);
        UNPROTECT(1);
    }
    // 6: is_closed
    {
        SEXP s = PROTECT(Rf_ScalarLogical(path_data.is_closed ? TRUE : FALSE));
        SET_VECTOR_ELT(result, 6, s);
        UNPROTECT(1);
    }

    UNPROTECT(2); // result, names
    return result;
}

} // extern "C"
