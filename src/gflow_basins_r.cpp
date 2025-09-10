#include <R.h>
#include <Rinternals.h>
 // Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

#include "set_wgraph.hpp"  // for set_wgraph_t
#include "SEXP_cpp_conversion_utils.hpp"
#include "cpp_utils.hpp"

extern "C" {
    SEXP S_find_gflow_basins(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_min_basin_size,
        SEXP s_min_path_size,
        SEXP s_q_edge_thld
        );

    SEXP S_find_local_extrema(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_min_basin_size
        );

    SEXP S_create_basin_cx(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y
        );
}

/**
 * @brief R interface to set_wgraph_t::find_gflow_basins().
 *
 * @param s_y             Numeric vector of responses y (length n).
 * @param s_min_basin_size Integer scalar: minimum basin size.
 * @param s_min_path_size  Integer scalar: minimum path size for monotonicity.
 * @param s_q_edge_thld    Numeric scalar: edge-length quantile threshold.
 *
 * @return An R list with two components:
 *   \item lmin_basins A list of minima basins, each a list with:
 *     \item vertex 1-based index of extremum.
 *     \item value  Numeric value at that extremum.
 *     \item basin  Numeric matrix [m x 2] of {vertex, distance}.
 *   \item lmax_basins Similar list for maxima basins.
 */
SEXP S_find_gflow_basins(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_basin_size,
    SEXP s_min_path_size,
    SEXP s_q_edge_thld
    ) {
    // ---- Input validation and conversion ----
    std::vector<std::vector<int>> adj_list       = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    // Convert numeric vector directly
    double* y_ptr = REAL(s_y);
    std::vector<double> y(y_ptr, y_ptr + LENGTH(s_y));

    int min_basin_size = INTEGER(s_min_basin_size)[0];
    int min_path_size  = INTEGER(s_min_path_size)[0];
    double q_edge_thld = REAL(s_q_edge_thld)[0];

    // ---- Call C++ implementation ----
    set_wgraph_t graph(adj_list, weight_list);

    // Ensure graph diameter is available
    graph.compute_graph_diameter();

    auto basins = graph.find_gflow_basins(
        y,
        static_cast<size_t>(min_basin_size),
        static_cast<size_t>(min_path_size),
        q_edge_thld
        );
    const auto &lmin = basins.first;
    const auto &lmax = basins.second;

    // ---- Build R objects ----
    SEXP r_lmin_list = PROTECT(Rf_allocVector(VECSXP, lmin.size()));
    SEXP r_lmax_list = PROTECT(Rf_allocVector(VECSXP, lmax.size()));

    // Helper macro to populate one basin list
    auto fill_basin = [&](const basin_t &b, SEXP r_blist) {
        // Names: vertex, value, basin
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(names, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(names, 1, Rf_mkChar("value"));
        SET_STRING_ELT(names, 2, Rf_mkChar("basin"));
        Rf_setAttrib(r_blist, R_NamesSymbol, names);

        // vertex (1-based)
        SEXP r_vertex = PROTECT(Rf_ScalarInteger(
                                    static_cast<int>(b.reachability_map.ref_vertex) + 1
                                    ));
        // value
        SEXP r_value  = PROTECT(Rf_ScalarReal(b.value));

        // basin matrix
        size_t m = b.reachability_map.sorted_vertices.size();
        SEXP r_basin = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
        double *pr = REAL(r_basin);
        for (size_t i = 0; i < m; ++i) {
            const auto &vi = b.reachability_map.sorted_vertices[i];
            pr[i]     = vi.vertex + 1;   // row i, col 1
            pr[i + m] = vi.distance;     // row i, col 2
        }

        // Set list entries
        SET_VECTOR_ELT(r_blist, 0, r_vertex);
        SET_VECTOR_ELT(r_blist, 1, r_value);
        SET_VECTOR_ELT(r_blist, 2, r_basin);

        UNPROTECT(4); // names, r_vertex, r_value, r_basin
    };

    // Fill minima basins
    for (size_t i = 0; i < lmin.size(); ++i) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 3));
        fill_basin(lmin[i], r_blist);
        SET_VECTOR_ELT(r_lmin_list, i, r_blist);
        UNPROTECT(1); // r_blist
    }

    // Fill maxima basins
    for (size_t i = 0; i < lmax.size(); ++i) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 3));
        fill_basin(lmax[i], r_blist);
        SET_VECTOR_ELT(r_lmax_list, i, r_blist);
        UNPROTECT(1);
    }

    // Name the two top‐level components
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_result_names, 0, Rf_mkChar("lmin_basins"));
    SET_STRING_ELT(r_result_names, 1, Rf_mkChar("lmax_basins"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
    SET_VECTOR_ELT(r_result, 0, r_lmin_list);
    SET_VECTOR_ELT(r_result, 1, r_lmax_list);

    UNPROTECT(4); // r_lmin_list, r_lmax_list, result, r_result_names

    return r_result;
}

/**
 * @brief R interface to find local extrema (minima and maxima) on a weighted graph.
 *
 * This function constructs a weighted graph from R adjacency and weight lists, evaluates
 * a scalar function defined on the graph's vertices, and identifies local minima and maxima
 * basins based on the structure of the graph and function values. Each extremum is reported
 * together with its associated basin — a set of vertices and distances reachable under a
 * gradient flow regime.
 *
 * Each returned basin is a list containing:
 * - `vertex`: the 1-based index of the extremum vertex,
 * - `value`: the function value at the extremum,
 * - `basin`: a matrix with two columns. Each row corresponds to a vertex in the basin:
 *     - column 1: vertex index (1-based),
 *     - column 2: geodesic distance from the extremum.
 *
 * The final return value is a list with two components:
 * - `lmin_basins`: list of local minima basins.
 * - `lmax_basins`: list of local maxima basins.
 *
 * @param s_adj_list R list of integer vectors representing the adjacency list of the graph.
 * @param s_weight_list R list of numeric vectors representing the corresponding edge weights.
 * @param s_y R numeric vector of function values defined at each vertex.
 * @param s_min_basin_size R integer scalar specifying the minimum number of vertices in a basin.
 *
 * @return A named R list with components `lmin_basins` and `lmax_basins`, each a list of basin descriptors.
 */
SEXP S_find_local_extrema(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_min_basin_size
    ) {
    // Convert input from R format
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);
    std::vector<double> y(REAL(s_y), REAL(s_y) + Rf_length(s_y));
    size_t min_basin_size = static_cast<size_t>(INTEGER(s_min_basin_size)[0]);

    // Create graph
    set_wgraph_t graph(adj_list, weight_list);

    // Ensure graph diameter is available
    graph.compute_graph_diameter();

    // Find local extrema
    auto basins = graph.find_local_extrema(
        y,
        min_basin_size
        );

    const auto &lmin = basins.first;
    const auto &lmax = basins.second;

    // ---- Build R objects ----
    SEXP r_lmin_list = PROTECT(Rf_allocVector(VECSXP, lmin.size()));
    SEXP r_lmax_list = PROTECT(Rf_allocVector(VECSXP, lmax.size()));

    // Helper macro to populate one basin list
    auto fill_basin = [&](const basin_t &b, SEXP r_blist) {
        // Names: vertex, value, basin
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(names, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(names, 1, Rf_mkChar("value"));
        SET_STRING_ELT(names, 2, Rf_mkChar("basin"));
        Rf_setAttrib(r_blist, R_NamesSymbol, names);

        // vertex (1-based)
        SEXP r_vertex = PROTECT(Rf_ScalarInteger(
                                    static_cast<int>(b.reachability_map.ref_vertex) + 1
                                    ));
        // value
        SEXP r_value  = PROTECT(Rf_ScalarReal(b.value));

        // basin matrix
        size_t m = b.reachability_map.sorted_vertices.size();
        SEXP r_basin = PROTECT(Rf_allocMatrix(REALSXP, m, 2));
        double *pr = REAL(r_basin);
        for (size_t i = 0; i < m; ++i) {
            const auto &vi = b.reachability_map.sorted_vertices[i];
            pr[i]     = vi.vertex + 1;   // row i, col 1
            pr[i + m] = vi.distance;     // row i, col 2
        }

        // Set list entries
        SET_VECTOR_ELT(r_blist, 0, r_vertex);
        SET_VECTOR_ELT(r_blist, 1, r_value);
        SET_VECTOR_ELT(r_blist, 2, r_basin);

        UNPROTECT(4); // names, r_vertex, r_value, r_basin
    };

    // Fill minima basins
    for (size_t i = 0; i < lmin.size(); ++i) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 3));
        fill_basin(lmin[i], r_blist);
        SET_VECTOR_ELT(r_lmin_list, i, r_blist);
        UNPROTECT(1); // r_blist
    }

    // Fill maxima basins
    for (size_t i = 0; i < lmax.size(); ++i) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 3));
        fill_basin(lmax[i], r_blist);
        SET_VECTOR_ELT(r_lmax_list, i, r_blist);
        UNPROTECT(1);
    }

    // Name the two top‐level components
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_result_names, 0, Rf_mkChar("lmin_basins"));
    SET_STRING_ELT(r_result_names, 1, Rf_mkChar("lmax_basins"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
    SET_VECTOR_ELT(r_result, 0, r_lmin_list);
    SET_VECTOR_ELT(r_result, 1, r_lmax_list);

    UNPROTECT(4); // r_lmin_list, r_lmax_list, result, r_result_names

    return r_result;
}

/**
 * @brief R interface to create a gradient flow basin complex on a weighted graph.
 *
 * This function constructs a topological denoising of a scalar function defined on a graph through
 * the identification and simplification of gradient flow basins. It implements the basin complex
 * algorithm described in the theoretical framework, which:
 *
 * 1. Identifies local extrema (minima and maxima) on the graph
 * 2. Constructs monotonic basins around each extremum
 * 3. Calculates the minimum monotonicity span for each basin
 * 4. Identifies cancellation pairs of basins
 * 5. Absorbs less significant basins based on the rel_min_monotonicity_span_thld parameter
 * 6. Performs harmonic repair of function values in absorbed regions
 *
 * The returned basin complex includes:
 * - A matrix summarizing all basins (basins_matrix)
 * - Harmonically repaired function values (topologically simplified)
 * - Remaining local minima basins after simplification
 * - Remaining local maxima basins after simplification
 * - Initial relative monotonicity spans
 * - Pairwise geodesic distance matrices between local minima and between local maxima
 *
 * @param s_adj_list R list of integer vectors representing the graph's adjacency list
 * @param s_weight_list R list of numeric vectors with corresponding edge weights
 * @param s_y R numeric vector of function values at each vertex
 * @param s_rel_min_monotonicity_span_thld R numeric scalar threshold - basins with relative
 *        monotonicity span below this threshold will be simplified
 *
 * @return A named R list with components:
 *   - basins_matrix: Matrix summarizing all basins with columns:
 *       * extremum_vertex: 1-based index of the extremum vertex
 *       * value: Function value at the extremum
 *       * is_maximum: 1 for maxima, 0 for minima
 *       * rel_min_span: Relative minimum monotonicity span
 *       * rel_max_span: Relative maximum monotonicity span
 *       * delta_rel_span: Difference between max and min spans
 *       * size: Number of vertices in the basin
 *       * rel_size: Size relative to total graph size
 *   - harmonic_predictions: Numeric vector of repaired function values
 *   - lmin_basins: List of local minima basins after simplification
 *   - lmax_basins: List of local maxima basins after simplification
 *   - init_rel_min_monotonicity_spans: Initial relative monotonicity spans
 *   - lmin_dist_mat: Matrix of pairwise geodesic distances between all local minima
 *   - lmax_dist_mat: Matrix of pairwise geodesic distances between all local maxima
 *
 * Each basin in lmin_basins and lmax_basins is a list with components:
 *   - vertex: 1-based index of the extremum vertex
 *   - value: Function value at the extremum
 *   - basin: Matrix [m x 3] with vertex index, vertex distance from extremum and the value of y at the vertex
 *   - basin_bd: Matrix [b x 3] of boundary vertices, distance from extremum, and their monotonicity spans
 *
 * The distance matrices have row and column names corresponding to the basin labels
 * (e.g., "m1", "m2", ... for minima and "M1", "M2", ... for maxima).
 *
 * @see create_basin_cx
 * @see find_local_extremum_basin
 * @see basin_t
 * @see basin_cx_t
 * @see compute_shortest_path_distances
 */
SEXP S_create_basin_cx(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y
    ) {
    // Convert input from R format
    std::vector<std::vector<int>>    adj_list    = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    std::vector<double> y(REAL(s_y), REAL(s_y) + Rf_length(s_y));

    // Create graph
    set_wgraph_t graph(adj_list, weight_list);

    // Ensure graph diameter is available
    graph.compute_graph_diameter();

    // Create basins
    basin_cx_t basin_cx = graph.create_basin_cx(y);

    // ---- Build R objects ----
    const char* names[] = {
        "basins_matrix",
        "harmonic_predictions",
        "lmin_basins",
        "lmax_basins",
        "init_rel_min_monotonicity_spans",
        "lmin_dist_mat",
        "lmax_dist_mat",
        "graph_diameter",
        NULL
    };

    const auto& lmin_map = basin_cx.lmin_basins_map;
    const auto& lmax_map = basin_cx.lmax_basins_map;

    int n_elements = 0;
    while (names[n_elements] != NULL) n_elements++;

    size_t n_protected = 0;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, n_elements)); n_protected++;
    SEXP r_result_names = PROTECT(Rf_allocVector(STRSXP, n_elements));
    // Set names
    for (int i = 0; i < n_elements; i++) {
        SET_STRING_ELT(r_result_names, i, Rf_mkChar(names[i]));
    }
    Rf_setAttrib(r_result, R_NamesSymbol, r_result_names);
    UNPROTECT(1); // for r_result_names

    auto vec_real = [&](const std::vector<double>& v) {
        SEXP ans = PROTECT(Rf_allocVector(REALSXP, v.size()));
        std::copy(v.begin(), v.end(), REAL(ans));
        return ans;
    };

    SEXP r_lmin_list = PROTECT(Rf_allocVector(VECSXP, lmin_map.size())); n_protected++;
    SEXP r_lmax_list = PROTECT(Rf_allocVector(VECSXP, lmax_map.size())); n_protected++;

    // Helper lambda to populate one basin list - MODIFIED to include basin_bd
    auto fill_basin = [&](const basin_t &b, SEXP r_blist) {
        // Names: vertex, value, basin, basin_bd (MODIFIED)
        SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
        SET_STRING_ELT(names, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(names, 1, Rf_mkChar("value"));
        SET_STRING_ELT(names, 2, Rf_mkChar("basin"));
        SET_STRING_ELT(names, 3, Rf_mkChar("basin_bd"));
        Rf_setAttrib(r_blist, R_NamesSymbol, names);

        // vertex (1-based)
        SEXP r_vertex = PROTECT(Rf_ScalarInteger(
                                    static_cast<int>(b.reachability_map.ref_vertex) + 1
                                    ));
        // value
        SEXP r_value  = PROTECT(Rf_ScalarReal(b.value));


        // basin matrix
        size_t m = b.reachability_map.sorted_vertices.size();

        // Create a copy of the sorted_vertices that we can sort by distance
        std::vector<vertex_info_t> sorted_by_distance = b.reachability_map.sorted_vertices;

        // Sort by distance (ascending)
        std::sort(sorted_by_distance.begin(), sorted_by_distance.end(),
                  [](const vertex_info_t& a, const vertex_info_t& b) {
                      return a.distance < b.distance;
                  });

        SEXP r_basin = PROTECT(Rf_allocMatrix(REALSXP, m, 3));
        double *pr = REAL(r_basin);
        for (size_t i = 0; i < m; ++i) {
            const auto &vi = sorted_by_distance[i];
            pr[i]         = vi.vertex + 1;   // row i, col 1
            pr[i + m]     = vi.distance;     // row i, col 2
            pr[i + 2 * m] = y[vi.vertex];    // row i, col 3
        }

        // Add column names to basin matrix
        SEXP basin_colnames = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(basin_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_colnames, 1, Rf_mkChar("distance"));
        SET_STRING_ELT(basin_colnames, 2, Rf_mkChar("value"));
        Rf_setAttrib(r_basin, R_DimNamesSymbol,
                  Rf_list2(R_NilValue, basin_colnames)); // No row names
        UNPROTECT(1); // basin_colnames


        // Basin boundary matrix from boundary_monotonicity_spans_map
        size_t bd_size = b.boundary_monotonicity_spans_map.size();
        SEXP r_basin_bd = PROTECT(Rf_allocMatrix(REALSXP, bd_size, 3));
        double *bd_pr = REAL(r_basin_bd);
        size_t bd_idx = 0;
        for (const auto& [vertex, span] : b.boundary_monotonicity_spans_map) {
            bd_pr[bd_idx]               = vertex + 1;                  // row i, col 1 (1-based vertex index)
            bd_pr[bd_idx + bd_size]     = b.reachability_map.distances.at(vertex);  // row i, col 2 (distance)
            bd_pr[bd_idx + 2 * bd_size] = span;                        // row i, col 3 (monotonicity span)
            bd_idx++;
        }

        // Add column names to basin_bd matrix
        SEXP basin_bd_colnames = PROTECT(Rf_allocVector(STRSXP, 3));
        SET_STRING_ELT(basin_bd_colnames, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(basin_bd_colnames, 1, Rf_mkChar("distance"));
        SET_STRING_ELT(basin_bd_colnames, 2, Rf_mkChar("span"));
        Rf_setAttrib(r_basin_bd, R_DimNamesSymbol,
                  Rf_list2(R_NilValue, basin_bd_colnames)); // No row names
        UNPROTECT(1); // basin_bd_colnames

        // Set list entries
        SET_VECTOR_ELT(r_blist, 0, r_vertex);
        SET_VECTOR_ELT(r_blist, 1, r_value);
        SET_VECTOR_ELT(r_blist, 2, r_basin);
        SET_VECTOR_ELT(r_blist, 3, r_basin_bd);

        UNPROTECT(5); // names, r_vertex, r_value, r_basin, r_basin_bd
    };

    // Fill minima basins
    size_t counter = 0;
    for (const auto& [v, basin] : lmin_map) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 4)); // MODIFIED: now 4 elements
        fill_basin(basin, r_blist);
        SET_VECTOR_ELT(r_lmin_list, counter++, r_blist);
        UNPROTECT(1); // r_blist
    }

    // Fill maxima basins
    counter = 0;
    for (const auto& [v, basin] : lmax_map) {
        SEXP r_blist = PROTECT(Rf_allocVector(VECSXP, 4)); // MODIFIED: now 4 elements
        fill_basin(basin, r_blist);
        SET_VECTOR_ELT(r_lmax_list, counter++, r_blist);
        UNPROTECT(1);
    }

    // Create basins_matrix
    // Count total number of basins (minima + maxima)
    size_t total_basins = lmin_map.size() + lmax_map.size();

    // Column names for basins_matrix
    const char* basin_matrix_colnames[] = {
        "extremum_vertex",
        "hop_idx",
        "value",
        "is_maximum",
        "rel_min_span",
        "rel_max_span",
        "delta_rel_span",
        "size",
        "rel_size"
    };
    int n_cols = sizeof(basin_matrix_colnames) / sizeof(basin_matrix_colnames[0]);

    // Create the matrix
    SEXP basins_matrix = PROTECT(Rf_allocMatrix(REALSXP, total_basins, n_cols)); n_protected++;

    // Set column names
    SEXP colnames = PROTECT(Rf_allocVector(STRSXP, n_cols));
    for (int i = 0; i < n_cols; i++) {
        SET_STRING_ELT(colnames, i, Rf_mkChar(basin_matrix_colnames[i]));
    }
    Rf_setAttrib(basins_matrix, R_DimNamesSymbol,
              Rf_list2(R_NilValue, colnames)); // No row names
    UNPROTECT(1); // colnames

    // Fill the matrix with data from all basins
    double* matrix_data = REAL(basins_matrix);
    size_t row = 0;

    // Calculate total graph size for relative size calculations
    double total_graph_size = static_cast<double>(graph.num_vertices());

    // First process all minima
    for (const auto& [vertex, basin] : lmin_map) {
        size_t size           = basin.reachability_map.sorted_vertices.size();
        double rel_size       = total_graph_size > 0 ? static_cast<double>(size) / total_graph_size : 0.0;
        double rel_min_span   = basin.rel_min_monotonicity_span;
        double rel_max_span   = basin.rel_max_monotonicity_span;
        double delta_rel_span = rel_max_span - rel_min_span;

        // Fill row data
        matrix_data[row]                  = vertex + 1;         // extremum_vertex (1-based)
        matrix_data[row +   total_basins] = basin.extremum_hop_index;  // extremum_hop_index
        matrix_data[row + 2*total_basins] = basin.value;        // value
        matrix_data[row + 3*total_basins] = 0.0;                // is_maximum (0 = false)
        matrix_data[row + 4*total_basins] = rel_min_span;       // rel_min_span
        matrix_data[row + 5*total_basins] = rel_max_span;       // rel_max_span
        matrix_data[row + 6*total_basins] = delta_rel_span;     // delta_rel_span
        matrix_data[row + 7*total_basins] = static_cast<double>(size); // size
        matrix_data[row + 8*total_basins] = rel_size;           // rel_size

        row++;
    }

    // Then process all maxima
    for (const auto& [vertex, basin] : lmax_map) {
        size_t size           = basin.reachability_map.sorted_vertices.size();
        double rel_size       = total_graph_size > 0 ? static_cast<double>(size) / total_graph_size : 0.0;
        double rel_min_span   = basin.rel_min_monotonicity_span;
        double rel_max_span   = basin.rel_max_monotonicity_span;
        double delta_rel_span = rel_max_span - rel_min_span;

        // Fill row data
        matrix_data[row]                  = vertex + 1;         // extremum_vertex (1-based)
        matrix_data[row +   total_basins] = basin.extremum_hop_index;  // extremum_hop_index
        matrix_data[row + 2*total_basins] = basin.value;        // value
        matrix_data[row + 3*total_basins] = 1.0;                // is_maximum (1 = true)
        matrix_data[row + 4*total_basins] = rel_min_span;       // rel_min_span
        matrix_data[row + 5*total_basins] = rel_max_span;       // rel_max_span
        matrix_data[row + 6*total_basins] = delta_rel_span;     // delta_rel_span
        matrix_data[row + 7*total_basins] = static_cast<double>(size); // size
        matrix_data[row + 8*total_basins] = rel_size;           // rel_size

        row++;
    }

    // Calculate pairwise geodesic distances between local minima
    SEXP r_lmin_dist_mat = R_NilValue;
    if (lmin_map.size() > 0) {
        // Create a vector of local minima vertex indices
        std::vector<size_t> lmin_vertices;
        for (const auto& [vertex, _] : lmin_map) {
            lmin_vertices.push_back(vertex);
        }

        // Create distance matrix
        r_lmin_dist_mat = PROTECT(Rf_allocMatrix(REALSXP, lmin_vertices.size(), lmin_vertices.size())); n_protected++;
        double* lmin_dist_data = REAL(r_lmin_dist_mat);

        // Calculate distances for each pair
        for (size_t i = 0; i < lmin_vertices.size(); i++) {
            // Create a set of target vertices (all except current)
            std::unordered_set<size_t> targets;
            for (size_t j = 0; j < lmin_vertices.size(); j++) {
                if (i != j) {
                    targets.insert(lmin_vertices[j]);
                }
            }

            // Compute shortest path distances from this vertex to all targets
            auto distances = graph.compute_shortest_path_distances(lmin_vertices[i], targets);

            // Fill in matrix (including diagonal)
            for (size_t j = 0; j < lmin_vertices.size(); j++) {
                if (i == j) {
                    lmin_dist_data[i + j * lmin_vertices.size()] = 0.0; // Diagonal = 0
                } else {
                    lmin_dist_data[i + j * lmin_vertices.size()] = distances[lmin_vertices[j]];
                }
            }
        }

        // Add row and column names (using 1-based vertex indices for R)
        SEXP lmin_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP lmin_rownames = PROTECT(Rf_allocVector(STRSXP, lmin_vertices.size()));
        SEXP lmin_colnames = PROTECT(Rf_allocVector(STRSXP, lmin_vertices.size()));

        for (size_t i = 0; i < lmin_vertices.size(); i++) {
            char name[32];
            snprintf(name, sizeof(name), "m%zu", lmin_vertices[i] + 1); // 1-based index for R
            SET_STRING_ELT(lmin_rownames, i, Rf_mkChar(name));
            SET_STRING_ELT(lmin_colnames, i, Rf_mkChar(name));
        }

        SET_VECTOR_ELT(lmin_dimnames, 0, lmin_rownames);
        SET_VECTOR_ELT(lmin_dimnames, 1, lmin_colnames);
        Rf_setAttrib(r_lmin_dist_mat, R_DimNamesSymbol, lmin_dimnames);

        UNPROTECT(3); // lmin_dimnames, lmin_rownames, lmin_colnames
    }

    // Calculate pairwise geodesic distances between local maxima
    SEXP r_lmax_dist_mat = R_NilValue;
    if (lmax_map.size() > 0) {
        // Create a vector of local maxima vertex indices
        std::vector<size_t> lmax_vertices;
        for (const auto& [vertex, _] : lmax_map) {
            lmax_vertices.push_back(vertex);
        }

        // Create distance matrix
        r_lmax_dist_mat = PROTECT(Rf_allocMatrix(REALSXP, lmax_vertices.size(), lmax_vertices.size())); n_protected++;
        double* lmax_dist_data = REAL(r_lmax_dist_mat);

        // Calculate distances for each pair
        for (size_t i = 0; i < lmax_vertices.size(); i++) {
            // Create a set of target vertices (all except current)
            std::unordered_set<size_t> targets;
            for (size_t j = 0; j < lmax_vertices.size(); j++) {
                if (i != j) {
                    targets.insert(lmax_vertices[j]);
                }
            }

            // Compute shortest path distances from this vertex to all targets
            auto distances = graph.compute_shortest_path_distances(lmax_vertices[i], targets);

            // Fill in matrix (including diagonal)
            for (size_t j = 0; j < lmax_vertices.size(); j++) {
                if (i == j) {
                    lmax_dist_data[i + j * lmax_vertices.size()] = 0.0; // Diagonal = 0
                } else {
                    lmax_dist_data[i + j * lmax_vertices.size()] = distances[lmax_vertices[j]];
                }
            }
        }

        // Add row and column names (using 1-based vertex indices for R)
        SEXP lmax_dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SEXP lmax_rownames = PROTECT(Rf_allocVector(STRSXP, lmax_vertices.size()));
        SEXP lmax_colnames = PROTECT(Rf_allocVector(STRSXP, lmax_vertices.size()));

        for (size_t i = 0; i < lmax_vertices.size(); i++) {
            char name[32];
            snprintf(name, sizeof(name), "M%zu", lmax_vertices[i] + 1); // 1-based index for R
            SET_STRING_ELT(lmax_rownames, i, Rf_mkChar(name));
            SET_STRING_ELT(lmax_colnames, i, Rf_mkChar(name));
        }

        SET_VECTOR_ELT(lmax_dimnames, 0, lmax_rownames);
        SET_VECTOR_ELT(lmax_dimnames, 1, lmax_colnames);
        Rf_setAttrib(r_lmax_dist_mat, R_DimNamesSymbol, lmax_dimnames);

        UNPROTECT(3); // lmax_dimnames, lmax_rownames, lmax_colnames
    }

    SEXP r_graph_diameter = PROTECT(Rf_allocVector(REALSXP, 1)); n_protected++;
    REAL(r_graph_diameter)[0] = graph.graph_diameter;

    // Set results
    size_t i = 0;
    SET_VECTOR_ELT(r_result, i++, basins_matrix);
    SET_VECTOR_ELT(r_result, i++, vec_real(basin_cx.harmonic_predictions)); n_protected++;
    SET_VECTOR_ELT(r_result, i++, r_lmin_list);
    SET_VECTOR_ELT(r_result, i++, r_lmax_list);
    SET_VECTOR_ELT(r_result, i++, vec_real(basin_cx.init_rel_min_monotonicity_spans)); n_protected++;
    SET_VECTOR_ELT(r_result, i++, r_lmin_dist_mat);
    SET_VECTOR_ELT(r_result, i++, r_lmax_dist_mat);
    SET_VECTOR_ELT(r_result, i++, r_graph_diameter);

    UNPROTECT(n_protected);
    return r_result;
}


