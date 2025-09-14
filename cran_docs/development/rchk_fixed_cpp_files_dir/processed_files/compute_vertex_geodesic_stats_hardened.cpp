// compute_vertex_geodesic_stats_hardened.cpp — rchk-hardened drop-in
//
// Hardenings applied:
// - Defensive scalar extraction via Rf_as* with full NA/range checks.
// - 1-based → 0-based conversion for grid_vertex with bounds checking.
// - Long-vector safety: use R_xlen_t for sizes when allocating/filling R objects.
// - Result assembly uses local PROTECT scopes and keeps only `result` + `names` protected until tail → final UNPROTECT(2).
// - Avoids variable UNPROTECT counters and protects all temporaries before insertion.
//
// Assumptions: the following helpers and types exist in your build:
//   convert_adj_list_from_R, convert_weight_list_from_R
//   uniform_grid_graph_t (type), create_maximal_packing(...)
//   compute_vertex_geodesic_stats(uniform_grid_graph_t,...)
//   uniform_grid_graph_t exposes: graph_diameter (double) and get_vertex_eccentricity(size_t).
//
// Keep the function signature unchanged.

#include <vector>
#include <numeric>
#include <R.h>
#include <Rinternals.h>

// Forward decls for helpers (assumed available at link time)
extern std::vector<std::vector<int>>    convert_adj_list_from_R(SEXP);
extern std::vector<std::vector<double>> convert_weight_list_from_R(SEXP);

// Core API (assumed available)
class uniform_grid_graph_t {
public:
    double graph_diameter;
    std::pair<size_t,double> get_vertex_eccentricity(size_t v);
    uniform_grid_graph_t() = default;
    uniform_grid_graph_t(const std::vector<std::vector<int>>& adj,
                         const std::vector<std::vector<double>>& w,
                         const std::vector<size_t>& packing);
};
extern uniform_grid_graph_t create_maximal_packing(const std::vector<std::vector<int>>& adj_list,
                                                   const std::vector<std::vector<double>>& weight_list,
                                                   size_t n_packing_vertices,
                                                   size_t max_packing_iterations,
                                                   double packing_precision);

// compute_vertex_geodesic_stats signature is intentionally omitted; we rely on 'auto' return with structured bindings.

extern "C" SEXP S_compute_vertex_geodesic_stats(
    SEXP adj_list_sexp,
    SEXP weight_list_sexp,
    SEXP grid_vertex_sexp,
    SEXP min_radius_sexp,
    SEXP max_radius_sexp,
    SEXP n_steps_sexp,
    SEXP n_packing_vertices_sexp,
    SEXP packing_precision_sexp
) {
    R_CheckUserInterrupt();

    // Graph inputs (container-first)
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(adj_list_sexp);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(weight_list_sexp);
    const size_t n_vertices = adj_list.size();

    // Defensive scalars
    const int    grid_vertex_i       = Rf_asInteger(grid_vertex_sexp);
    const double min_radius          = Rf_asReal(min_radius_sexp);
    const double max_radius          = Rf_asReal(max_radius_sexp);
    const int    n_steps_i           = Rf_asInteger(n_steps_sexp);
    const int    n_pack_i            = Rf_asInteger(n_packing_vertices_sexp);
    const double packing_precision   = Rf_asReal(packing_precision_sexp);

    // NA / range checks
    if (grid_vertex_i == NA_INTEGER) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'grid_vertex' cannot be NA.");
    }
    if (n_steps_i == NA_INTEGER || n_steps_i <= 0) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'n_steps' must be a positive integer.");
    }
    if (ISNAN(min_radius) || ISNAN(max_radius)) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'min_radius' and 'max_radius' cannot be NA.");
    }
    if (min_radius < 0.0) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'min_radius' must be >= 0.");
    }
    if (max_radius <= 0.0) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'max_radius' must be > 0.");
    }
    if (max_radius < min_radius) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'max_radius' must be >= 'min_radius'.");
    }
    if (n_pack_i == NA_INTEGER || n_pack_i < 0) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'n_packing_vertices' must be >= 0.");
    }
    if (ISNAN(packing_precision) || packing_precision <= 0.0) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'packing_precision' must be > 0.");
    }

    // 1-based (R) → 0-based (C++) with bounds check
    const size_t grid_vertex = static_cast<size_t>(grid_vertex_i - 1);
    if (grid_vertex_i <= 0 || grid_vertex >= n_vertices) {
        Rf_error("S_compute_vertex_geodesic_stats(): 'grid_vertex' out of bounds (1..%zu).", n_vertices);
    }

    // Build uniform grid graph (packing or full)
    uniform_grid_graph_t grid_graph;
    const size_t n_packing_vertices = static_cast<size_t>(n_pack_i);

    if (n_packing_vertices > 0 && n_packing_vertices < n_vertices) {
        const size_t max_packing_iterations = 20;
        grid_graph = create_maximal_packing(adj_list, weight_list,
                                            n_packing_vertices,
                                            max_packing_iterations,
                                            packing_precision);
    } else {
        // Use all vertices as packing
        std::vector<size_t> packing(n_vertices);
        std::iota(packing.begin(), packing.end(), static_cast<size_t>(0));
        grid_graph = uniform_grid_graph_t(adj_list, weight_list, packing);

        // Initialize diameter info (two-sweep)
        auto ecc0         = grid_graph.get_vertex_eccentricity(static_cast<size_t>(0));
        auto ecc_far      = grid_graph.get_vertex_eccentricity(ecc0.first);
        grid_graph.graph_diameter = ecc_far.second;
    }

    // Compute vertex-specific geodesic statistics
    auto stats = compute_vertex_geodesic_stats(
        grid_graph,
        grid_vertex,
        min_radius,
        max_radius,
        n_steps_i
    );

    const R_xlen_t n_radii = static_cast<R_xlen_t>(stats.size());

    // Result list (2 elements): data matrix [n_radii x 10], vertex (1-based)
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));

    // 0: data matrix
    {
        SEXP mat = PROTECT(Rf_allocMatrix(REALSXP, n_radii, 10));
        double* ptr = REAL(mat);

        for (R_xlen_t i = 0; i < n_radii; ++i) {
            auto [radius, geodesic_rays, composite_geodesics, overlap_stats] = stats[static_cast<size_t>(i)];

            // Column-major fill
            ptr[i + 0 * n_radii] = radius;
            ptr[i + 1 * n_radii] = static_cast<double>(geodesic_rays);
            ptr[i + 2 * n_radii] = static_cast<double>(composite_geodesics);

            ptr[i + 3 * n_radii] = overlap_stats.min;
            ptr[i + 4 * n_radii] = overlap_stats.p05;
            ptr[i + 5 * n_radii] = overlap_stats.p25;
            ptr[i + 6 * n_radii] = overlap_stats.median;
            ptr[i + 7 * n_radii] = overlap_stats.p75;
            ptr[i + 8 * n_radii] = overlap_stats.p95;
            ptr[i + 9 * n_radii] = overlap_stats.max;
        }

        // Column names
        SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 10));
        SET_STRING_ELT(colnames, 0, Rf_mkChar("radius"));
        SET_STRING_ELT(colnames, 1, Rf_mkChar("rays"));
        SET_STRING_ELT(colnames, 2, Rf_mkChar("composite_geodesics"));
        SET_STRING_ELT(colnames, 3, Rf_mkChar("overlap_min"));
        SET_STRING_ELT(colnames, 4, Rf_mkChar("overlap_p05"));
        SET_STRING_ELT(colnames, 5, Rf_mkChar("overlap_p25"));
        SET_STRING_ELT(colnames, 6, Rf_mkChar("overlap_median"));
        SET_STRING_ELT(colnames, 7, Rf_mkChar("overlap_p75"));
        SET_STRING_ELT(colnames, 8, Rf_mkChar("overlap_p95"));
        SET_STRING_ELT(colnames, 9, Rf_mkChar("overlap_max"));

        // dimnames: list(NULL, colnames)
        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, R_NilValue);
        SET_VECTOR_ELT(dimnames, 1, colnames);
        Rf_setAttrib(mat, R_DimNamesSymbol, dimnames);
        UNPROTECT(2); // colnames, dimnames

        SET_VECTOR_ELT(result, 0, mat);
        UNPROTECT(1); // mat
    }

    // 1: vertex (1-based echo)
    {
        SEXP v = PROTECT(Rf_ScalarInteger(static_cast<int>(grid_vertex + 1)));
        SET_VECTOR_ELT(result, 1, v);
        UNPROTECT(1);
    }

    // names
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, Rf_mkChar("data"));
    SET_STRING_ELT(names, 1, Rf_mkChar("vertex"));
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(2); // result, names
    return result;
}
