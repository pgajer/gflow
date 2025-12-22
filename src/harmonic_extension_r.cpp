/**
 * @file harmonic_extension_r.cpp
 * @brief SEXP interface for harmonic extension computation
 */

#include <R.h>
#include <Rinternals.h>
#include "harmonic_extension.hpp"
#include "set_wgraph.hpp"

#include <vector>
#include <unordered_set>

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);
size_t select_max_density_trajectory(
    const std::vector<std::vector<size_t>>& trajectories,
    const std::vector<double>& density
	);

/**
 * @brief Compute harmonic extension from R
 *
 * @param s_adj_list Adjacency list (0-based from R conversion)
 * @param s_weight_list Edge weight list
 * @param s_trajectory Trajectory vertices (0-based)
 * @param s_tube_radius Tubular neighborhood radius in hops
 * @param s_use_edge_weights Use inverse edge length weights
 * @param s_max_iterations Maximum solver iterations
 * @param s_tolerance Convergence tolerance
 * @param s_basin_restriction Optional basin vertices to restrict neighborhood
 * @param s_verbose Print progress
 *
 * @return R list with harmonic extension results
 */
extern "C" SEXP S_compute_harmonic_extension(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_trajectory,
    SEXP s_tube_radius,
    SEXP s_use_edge_weights,
    SEXP s_max_iterations,
    SEXP s_tolerance,
    SEXP s_basin_restriction,
    SEXP s_verbose
) {
    // Convert graph structure
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    // Convert trajectory (already 0-based from R)
    const int n_traj = LENGTH(s_trajectory);
    std::vector<size_t> trajectory(n_traj);
    const int* p_traj = INTEGER(s_trajectory);
    for (int i = 0; i < n_traj; ++i) {
        trajectory[i] = static_cast<size_t>(p_traj[i]);
    }

    // Build parameters
    harmonic_extension_params_t params;
    params.tube_radius = Rf_asInteger(s_tube_radius);
    params.use_edge_weights = Rf_asLogical(s_use_edge_weights);
    params.max_iterations = Rf_asInteger(s_max_iterations);
    params.tolerance = Rf_asReal(s_tolerance);

    // Convert basin restriction if provided
    if (!Rf_isNull(s_basin_restriction) && LENGTH(s_basin_restriction) > 0) {
        const int* p_basin = INTEGER(s_basin_restriction);
        const int n_basin = LENGTH(s_basin_restriction);
        for (int i = 0; i < n_basin; ++i) {
            params.basin_restriction.insert(static_cast<size_t>(p_basin[i]));
        }
    }

    bool verbose = Rf_asLogical(s_verbose);

    // Compute harmonic extension
    harmonic_extension_result_t result = graph.compute_harmonic_extension(
        trajectory, params, verbose
    );

    // Convert result to R list
    const int n_components = 10;
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // trajectory (1-based)
    const int n_traj_out = static_cast<int>(result.trajectory.size());
    SEXP s_traj_out = PROTECT(Rf_allocVector(INTSXP, n_traj_out));
    int* p_traj_out = INTEGER(s_traj_out);
    for (int i = 0; i < n_traj_out; ++i) {
        p_traj_out[i] = static_cast<int>(result.trajectory[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("trajectory"));
    SET_VECTOR_ELT(s_result, idx++, s_traj_out);
    UNPROTECT(1);

    // trajectory.coords
    const int n_tc = static_cast<int>(result.trajectory_coords.size());
    SEXP s_tc = PROTECT(Rf_allocVector(REALSXP, n_tc));
    double* p_tc = REAL(s_tc);
    for (int i = 0; i < n_tc; ++i) {
        p_tc[i] = result.trajectory_coords[i];
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("trajectory.coords"));
    SET_VECTOR_ELT(s_result, idx++, s_tc);
    UNPROTECT(1);

    // trajectory.length
    SET_STRING_ELT(s_names, idx, Rf_mkChar("trajectory.length"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarReal(result.trajectory_length));

    // tubular.vertices (1-based)
    const int n_tv = static_cast<int>(result.tubular_vertices.size());
    SEXP s_tv = PROTECT(Rf_allocVector(INTSXP, n_tv));
    int* p_tv = INTEGER(s_tv);
    for (int i = 0; i < n_tv; ++i) {
        p_tv[i] = static_cast<int>(result.tubular_vertices[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("tubular.vertices"));
    SET_VECTOR_ELT(s_result, idx++, s_tv);
    UNPROTECT(1);

    // hop.distances
    SEXP s_hd = PROTECT(Rf_allocVector(INTSXP, n_tv));
    int* p_hd = INTEGER(s_hd);
    for (int i = 0; i < n_tv; ++i) {
        p_hd[i] = result.hop_distances[i];
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("hop.distances"));
    SET_VECTOR_ELT(s_result, idx++, s_hd);
    UNPROTECT(1);

	// nearest.traj.idx (1-based trajectory index)
    SEXP s_nti = PROTECT(Rf_allocVector(INTSXP, n_tv));
    int* p_nti = INTEGER(s_nti);
    for (int i = 0; i < n_tv; ++i) {
        p_nti[i] = static_cast<int>(result.nearest_traj_idx[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("nearest.traj.idx"));
    SET_VECTOR_ELT(s_result, idx++, s_nti);
    UNPROTECT(1);

    // extended.coords
    SEXP s_ec = PROTECT(Rf_allocVector(REALSXP, n_tv));
    double* p_ec = REAL(s_ec);
    for (int i = 0; i < n_tv; ++i) {
        p_ec[i] = result.extended_coords[i];
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("extended.coords"));
    SET_VECTOR_ELT(s_result, idx++, s_ec);
    UNPROTECT(1);

    // n.iterations
    SET_STRING_ELT(s_names, idx, Rf_mkChar("n.iterations"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarInteger(result.n_iterations));

    // final.max.change
    SET_STRING_ELT(s_names, idx, Rf_mkChar("final.max.change"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarReal(result.final_max_change));

    // tube.radius (for reference)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("tube.radius"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarInteger(params.tube_radius));

    Rf_setAttrib(s_result, R_NamesSymbol, s_names);

    // Add class
    SEXP s_class = PROTECT(Rf_mkString("harmonic_extension"));
    Rf_setAttrib(s_result, R_ClassSymbol, s_class);

    UNPROTECT(3);  // s_result, s_names, s_class
    return s_result;
}

/**
 * @brief Select trajectory with maximal mean density
 *
 * @param s_trajectories List of trajectory vectors (0-based vertices)
 * @param s_density Density values for all vertices
 *
 * @return Index of best trajectory (1-based for R)
 */
extern "C" SEXP S_select_max_density_trajectory(
    SEXP s_trajectories,
    SEXP s_density
) {
    const int n_traj = LENGTH(s_trajectories);

    if (n_traj == 0) {
        return Rf_ScalarInteger(NA_INTEGER);
    }

    // Convert density
    const int n_dens = LENGTH(s_density);
    std::vector<double> density(n_dens);
    const double* p_dens = REAL(s_density);
    for (int i = 0; i < n_dens; ++i) {
        density[i] = p_dens[i];
    }

    // Convert trajectories
    std::vector<std::vector<size_t>> trajectories(n_traj);
    for (int i = 0; i < n_traj; ++i) {
        SEXP s_traj = VECTOR_ELT(s_trajectories, i);
        const int len = LENGTH(s_traj);
        const int* p_traj = INTEGER(s_traj);

        trajectories[i].resize(len);
        for (int j = 0; j < len; ++j) {
            trajectories[i][j] = static_cast<size_t>(p_traj[j]);
        }
    }

    size_t best_idx = select_max_density_trajectory(trajectories, density);

    // Return 1-based index
    return Rf_ScalarInteger(static_cast<int>(best_idx) + 1);
}
