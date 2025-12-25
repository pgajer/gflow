/**
 * @file madag_r.cpp
 * @brief R interface implementation for MADAG functions
 */

#include "madag_r.h"
#include "madag.hpp"
#include "set_wgraph.hpp"

#include <vector>
#include <string>

// Forward declarations for conversion utilities
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Convert MADAG structure to R list
 */
static SEXP madag_to_R(const madag_t& madag, const std::vector<double>& y) {
    // Count items for PROTECT
    int n_protect = 0;
    
    // Create output list with named elements
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 12));
    n_protect++;
    
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 12));
    n_protect++;
    
    // 1. source_vertex (convert to 1-based)
    SEXP s_source = PROTECT(Rf_ScalarInteger(static_cast<int>(madag.source_vertex) + 1));
    n_protect++;
    SET_VECTOR_ELT(result, 0, s_source);
    SET_STRING_ELT(names, 0, Rf_mkChar("source_vertex"));
    
    // 2. source_value
    SEXP s_source_val = PROTECT(Rf_ScalarReal(madag.source_value));
    n_protect++;
    SET_VECTOR_ELT(result, 1, s_source_val);
    SET_STRING_ELT(names, 1, Rf_mkChar("source_value"));
    
    // 3. reachable_vertices (convert to 1-based)
    size_t n_reach = madag.reachable_vertices.size();
    SEXP s_reachable = PROTECT(Rf_allocVector(INTSXP, n_reach));
    n_protect++;
    int* p_reachable = INTEGER(s_reachable);
    for (size_t i = 0; i < n_reach; ++i) {
        p_reachable[i] = static_cast<int>(madag.reachable_vertices[i]) + 1;
    }
    SET_VECTOR_ELT(result, 2, s_reachable);
    SET_STRING_ELT(names, 2, Rf_mkChar("reachable_vertices"));
    
    // 4. reachable_maxima (convert to 1-based)
    size_t n_max = madag.reachable_maxima.size();
    SEXP s_maxima = PROTECT(Rf_allocVector(INTSXP, n_max));
    n_protect++;
    int* p_maxima = INTEGER(s_maxima);
    for (size_t i = 0; i < n_max; ++i) {
        p_maxima[i] = static_cast<int>(madag.reachable_maxima[i]) + 1;
    }
    SET_VECTOR_ELT(result, 3, s_maxima);
    SET_STRING_ELT(names, 3, Rf_mkChar("reachable_maxima"));
    
    // 5. predecessors (list of integer vectors, 1-based)
    SEXP s_predecessors = PROTECT(Rf_allocVector(VECSXP, n_reach));
    n_protect++;
    SEXP s_pred_names = PROTECT(Rf_allocVector(STRSXP, n_reach));
    n_protect++;
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = madag.reachable_vertices[i];
        auto it = madag.predecessors.find(v);
        
        if (it != madag.predecessors.end() && !it->second.empty()) {
            SEXP pred_vec = PROTECT(Rf_allocVector(INTSXP, it->second.size()));
            int* p_pred = INTEGER(pred_vec);
            for (size_t j = 0; j < it->second.size(); ++j) {
                p_pred[j] = static_cast<int>(it->second[j]) + 1;
            }
            SET_VECTOR_ELT(s_predecessors, i, pred_vec);
            UNPROTECT(1);
        } else {
            SET_VECTOR_ELT(s_predecessors, i, Rf_allocVector(INTSXP, 0));
        }
        
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%d", static_cast<int>(v) + 1);
        SET_STRING_ELT(s_pred_names, i, Rf_mkChar(name_buf));
    }
    Rf_setAttrib(s_predecessors, R_NamesSymbol, s_pred_names);
    SET_VECTOR_ELT(result, 4, s_predecessors);
    SET_STRING_ELT(names, 4, Rf_mkChar("predecessors"));
    
    // 6. successors (list of integer vectors, 1-based)
    SEXP s_successors = PROTECT(Rf_allocVector(VECSXP, n_reach));
    n_protect++;
    SEXP s_succ_names = PROTECT(Rf_allocVector(STRSXP, n_reach));
    n_protect++;
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = madag.reachable_vertices[i];
        auto it = madag.successors.find(v);
        
        if (it != madag.successors.end() && !it->second.empty()) {
            SEXP succ_vec = PROTECT(Rf_allocVector(INTSXP, it->second.size()));
            int* p_succ = INTEGER(succ_vec);
            for (size_t j = 0; j < it->second.size(); ++j) {
                p_succ[j] = static_cast<int>(it->second[j]) + 1;
            }
            SET_VECTOR_ELT(s_successors, i, succ_vec);
            UNPROTECT(1);
        } else {
            SET_VECTOR_ELT(s_successors, i, Rf_allocVector(INTSXP, 0));
        }
        
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%d", static_cast<int>(v) + 1);
        SET_STRING_ELT(s_succ_names, i, Rf_mkChar(name_buf));
    }
    Rf_setAttrib(s_successors, R_NamesSymbol, s_succ_names);
    SET_VECTOR_ELT(result, 5, s_successors);
    SET_STRING_ELT(names, 5, Rf_mkChar("successors"));
    
    // 7. topological_order (convert to 1-based)
    size_t n_topo = madag.topological_order.size();
    SEXP s_topo = PROTECT(Rf_allocVector(INTSXP, n_topo));
    n_protect++;
    int* p_topo = INTEGER(s_topo);
    for (size_t i = 0; i < n_topo; ++i) {
        p_topo[i] = static_cast<int>(madag.topological_order[i]) + 1;
    }
    SET_VECTOR_ELT(result, 6, s_topo);
    SET_STRING_ELT(names, 6, Rf_mkChar("topological_order"));
    
    // 8. path_count_from_source
    SEXP s_path_from = PROTECT(Rf_allocVector(INTSXP, n_reach));
    n_protect++;
    SEXP s_path_from_names = PROTECT(Rf_allocVector(STRSXP, n_reach));
    n_protect++;
    int* p_path_from = INTEGER(s_path_from);
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = madag.reachable_vertices[i];
        auto it = madag.path_count_from_source.find(v);
        p_path_from[i] = (it != madag.path_count_from_source.end()) ? 
                         static_cast<int>(it->second) : 0;
        
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%d", static_cast<int>(v) + 1);
        SET_STRING_ELT(s_path_from_names, i, Rf_mkChar(name_buf));
    }
    Rf_setAttrib(s_path_from, R_NamesSymbol, s_path_from_names);
    SET_VECTOR_ELT(result, 7, s_path_from);
    SET_STRING_ELT(names, 7, Rf_mkChar("path_count_from_source"));
    
    // 9. path_count_to_sinks
    SEXP s_path_to = PROTECT(Rf_allocVector(INTSXP, n_reach));
    n_protect++;
    SEXP s_path_to_names = PROTECT(Rf_allocVector(STRSXP, n_reach));
    n_protect++;
    int* p_path_to = INTEGER(s_path_to);
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = madag.reachable_vertices[i];
        auto it = madag.path_count_to_sinks.find(v);
        p_path_to[i] = (it != madag.path_count_to_sinks.end()) ? 
                       static_cast<int>(it->second) : 0;
        
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%d", static_cast<int>(v) + 1);
        SET_STRING_ELT(s_path_to_names, i, Rf_mkChar(name_buf));
    }
    Rf_setAttrib(s_path_to, R_NamesSymbol, s_path_to_names);
    SET_VECTOR_ELT(result, 8, s_path_to);
    SET_STRING_ELT(names, 8, Rf_mkChar("path_count_to_sinks"));
    
    // 10. cells (list of cell structures)
    size_t n_cells = madag.cells.size();
    SEXP s_cells = PROTECT(Rf_allocVector(VECSXP, n_cells));
    n_protect++;
    for (size_t i = 0; i < n_cells; ++i) {
        const ms_cell_t& cell = madag.cells[i];
        
        SEXP cell_list = PROTECT(Rf_allocVector(VECSXP, 9));
        SEXP cell_names = PROTECT(Rf_allocVector(STRSXP, 9));
        
        // min_vertex (1-based)
        SET_VECTOR_ELT(cell_list, 0, Rf_ScalarInteger(static_cast<int>(cell.min_vertex) + 1));
        SET_STRING_ELT(cell_names, 0, Rf_mkChar("min_vertex"));
        
        // max_vertex (1-based)
        SET_VECTOR_ELT(cell_list, 1, Rf_ScalarInteger(static_cast<int>(cell.max_vertex) + 1));
        SET_STRING_ELT(cell_names, 1, Rf_mkChar("max_vertex"));
        
        // min_value
        SET_VECTOR_ELT(cell_list, 2, Rf_ScalarReal(cell.min_value));
        SET_STRING_ELT(cell_names, 2, Rf_mkChar("min_value"));
        
        // max_value
        SET_VECTOR_ELT(cell_list, 3, Rf_ScalarReal(cell.max_value));
        SET_STRING_ELT(cell_names, 3, Rf_mkChar("max_value"));
        
        // support (1-based)
        SEXP support = PROTECT(Rf_allocVector(INTSXP, cell.support.size()));
        int* p_support = INTEGER(support);
        for (size_t j = 0; j < cell.support.size(); ++j) {
            p_support[j] = static_cast<int>(cell.support[j]) + 1;
        }
        SET_VECTOR_ELT(cell_list, 4, support);
        SET_STRING_ELT(cell_names, 4, Rf_mkChar("support"));
        UNPROTECT(1);
        
        // n_trajectories
        SET_VECTOR_ELT(cell_list, 5, Rf_ScalarInteger(static_cast<int>(cell.n_trajectories)));
        SET_STRING_ELT(cell_names, 5, Rf_mkChar("n_trajectories"));
        
        // explicitly_enumerated
        SET_VECTOR_ELT(cell_list, 6, Rf_ScalarLogical(cell.explicitly_enumerated));
        SET_STRING_ELT(cell_names, 6, Rf_mkChar("explicitly_enumerated"));
        
        // bottlenecks (1-based)
        SEXP bottlenecks = PROTECT(Rf_allocVector(INTSXP, cell.bottlenecks.size()));
        int* p_bottlenecks = INTEGER(bottlenecks);
        for (size_t j = 0; j < cell.bottlenecks.size(); ++j) {
            p_bottlenecks[j] = static_cast<int>(cell.bottlenecks[j]) + 1;
        }
        SET_VECTOR_ELT(cell_list, 7, bottlenecks);
        SET_STRING_ELT(cell_names, 7, Rf_mkChar("bottlenecks"));
        UNPROTECT(1);
        
        // n_clusters
        SET_VECTOR_ELT(cell_list, 8, Rf_ScalarInteger(cell.n_clusters));
        SET_STRING_ELT(cell_names, 8, Rf_mkChar("n_clusters"));
        
        Rf_setAttrib(cell_list, R_NamesSymbol, cell_names);
        SET_VECTOR_ELT(s_cells, i, cell_list);
        UNPROTECT(2);  // cell_list, cell_names
    }
    SET_VECTOR_ELT(result, 9, s_cells);
    SET_STRING_ELT(names, 9, Rf_mkChar("cells"));
    
    // 11. n_vertices
    SET_VECTOR_ELT(result, 10, Rf_ScalarInteger(static_cast<int>(madag.n_vertices())));
    SET_STRING_ELT(names, 10, Rf_mkChar("n_vertices"));
    
    // 12. n_cells
    SET_VECTOR_ELT(result, 11, Rf_ScalarInteger(static_cast<int>(madag.n_cells())));
    SET_STRING_ELT(names, 11, Rf_mkChar("n_cells"));
    
    Rf_setAttrib(result, R_NamesSymbol, names);
    
    UNPROTECT(n_protect);
    return result;
}

/**
 * @brief Convert trajectory to R list
 */
static SEXP trajectory_to_R(const madag_trajectory_t& traj) {
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    
    // vertices (1-based)
    SEXP vertices = PROTECT(Rf_allocVector(INTSXP, traj.vertices.size()));
    int* p_vertices = INTEGER(vertices);
    for (size_t i = 0; i < traj.vertices.size(); ++i) {
        p_vertices[i] = static_cast<int>(traj.vertices[i]) + 1;
    }
    SET_VECTOR_ELT(result, 0, vertices);
    SET_STRING_ELT(names, 0, Rf_mkChar("vertices"));
    
    // source_min (1-based)
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(static_cast<int>(traj.source_min) + 1));
    SET_STRING_ELT(names, 1, Rf_mkChar("source_min"));
    
    // sink_max (1-based)
    SET_VECTOR_ELT(result, 2, Rf_ScalarInteger(static_cast<int>(traj.sink_max) + 1));
    SET_STRING_ELT(names, 2, Rf_mkChar("sink_max"));
    
    // total_ascent
    SET_VECTOR_ELT(result, 3, Rf_ScalarReal(traj.total_ascent));
    SET_STRING_ELT(names, 3, Rf_mkChar("total_ascent"));
    
    // cluster_id
    SET_VECTOR_ELT(result, 4, Rf_ScalarInteger(traj.cluster_id));
    SET_STRING_ELT(names, 4, Rf_mkChar("cluster_id"));
    
    Rf_setAttrib(result, R_NamesSymbol, names);
    
    UNPROTECT(3);
    return result;
}

/**
 * @brief Convert vector of trajectories to R list
 */
static SEXP trajectories_to_R(const std::vector<madag_trajectory_t>& trajectories) {
    size_t n = trajectories.size();
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n));
    
    for (size_t i = 0; i < n; ++i) {
        SEXP traj_r = PROTECT(trajectory_to_R(trajectories[i]));
        SET_VECTOR_ELT(result, i, traj_r);
        UNPROTECT(1);
    }
    
    UNPROTECT(1);
    return result;
}

// ============================================================================
// Main SEXP Functions
// ============================================================================

extern "C" {

SEXP S_construct_madag(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_source_vertex,
    SEXP s_params,
    SEXP s_verbose
) {
    // Convert inputs
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    
    // Source vertex (R is 1-based, convert to 0-based)
    int source_vertex_r = INTEGER(s_source_vertex)[0];
    size_t source_vertex = static_cast<size_t>(source_vertex_r - 1);
    
    bool verbose = LOGICAL(s_verbose)[0];
    
    // Extract parameters
    madag_params_t params;
    if (!Rf_isNull(s_params)) {
        SEXP s_max_traj = PROTECT(VECTOR_ELT(s_params, 0));
        SEXP s_compute_counts = PROTECT(VECTOR_ELT(s_params, 1));
        SEXP s_enumerate = PROTECT(VECTOR_ELT(s_params, 2));
        SEXP s_edge_thld = PROTECT(VECTOR_ELT(s_params, 3));
        SEXP s_min_support = PROTECT(VECTOR_ELT(s_params, 4));
        
        params.max_trajectories_per_cell = static_cast<size_t>(INTEGER(s_max_traj)[0]);
        params.compute_path_counts = LOGICAL(s_compute_counts)[0];
        params.enumerate_trajectories = LOGICAL(s_enumerate)[0];
        params.edge_length_quantile_thld = REAL(s_edge_thld)[0];
        params.min_cell_support = static_cast<size_t>(INTEGER(s_min_support)[0]);
        
        UNPROTECT(5);
    }
    
    // Build graph
    set_wgraph_t graph(adj_list, weight_list);
    
    // Construct MADAG
    madag_t madag = construct_madag(graph, y, source_vertex, params, verbose);
    
    // Convert to R
    return madag_to_R(madag, y);
}

SEXP S_enumerate_cell_trajectories(
    SEXP s_madag,
    SEXP s_y,
    SEXP s_max_vertex,
    SEXP s_max_trajectories
) {
    // This function requires reconstructing the MADAG from the R structure
    // For now, we'll implement a simpler version that works directly
    
    // Extract key information from s_madag
    SEXP s_source = VECTOR_ELT(s_madag, 0);  // source_vertex
    SEXP s_reachable = VECTOR_ELT(s_madag, 2);  // reachable_vertices
    SEXP s_predecessors = VECTOR_ELT(s_madag, 4);  // predecessors
    SEXP s_successors = VECTOR_ELT(s_madag, 5);  // successors
    
    size_t source_vertex = static_cast<size_t>(INTEGER(s_source)[0] - 1);
    size_t max_vertex = static_cast<size_t>(INTEGER(s_max_vertex)[0] - 1);
    size_t max_traj = static_cast<size_t>(INTEGER(s_max_trajectories)[0]);
    
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    
    // Reconstruct MADAG structure
    madag_t madag;
    madag.source_vertex = source_vertex;
    madag.source_value = y[source_vertex];
    
    // Reconstruct reachable vertices
    int* p_reachable = INTEGER(s_reachable);
    size_t n_reach = LENGTH(s_reachable);
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = static_cast<size_t>(p_reachable[i] - 1);
        madag.reachable_vertices.push_back(v);
        madag.reachable_set.insert(v);
    }
    
    // Reconstruct successors
    SEXP succ_names = Rf_getAttrib(s_successors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_successors); ++i) {
        const char* name = CHAR(STRING_ELT(succ_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP succ_vec = VECTOR_ELT(s_successors, i);
        std::vector<size_t> successors;
        int* p_succ = INTEGER(succ_vec);
        for (int j = 0; j < LENGTH(succ_vec); ++j) {
            successors.push_back(static_cast<size_t>(p_succ[j] - 1));
        }
        madag.successors[v] = successors;
    }
    
    // Reconstruct predecessors
    SEXP pred_names = Rf_getAttrib(s_predecessors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_predecessors); ++i) {
        const char* name = CHAR(STRING_ELT(pred_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP pred_vec = VECTOR_ELT(s_predecessors, i);
        std::vector<size_t> predecessors;
        int* p_pred = INTEGER(pred_vec);
        for (int j = 0; j < LENGTH(pred_vec); ++j) {
            predecessors.push_back(static_cast<size_t>(p_pred[j] - 1));
        }
        madag.predecessors[v] = predecessors;
    }
    
    // Enumerate trajectories
    std::vector<madag_trajectory_t> trajectories = 
        enumerate_cell_trajectories(madag, y, max_vertex, max_traj);
    
    // Convert to R
    return trajectories_to_R(trajectories);
}

SEXP S_sample_cell_trajectories(
    SEXP s_madag,
    SEXP s_y,
    SEXP s_max_vertex,
    SEXP s_n_samples,
    SEXP s_seed
) {
    // Similar reconstruction as above
    SEXP s_source = VECTOR_ELT(s_madag, 0);
    SEXP s_reachable = VECTOR_ELT(s_madag, 2);
    SEXP s_successors = VECTOR_ELT(s_madag, 5);
    SEXP s_predecessors = VECTOR_ELT(s_madag, 4);
    
    size_t source_vertex = static_cast<size_t>(INTEGER(s_source)[0] - 1);
    size_t max_vertex = static_cast<size_t>(INTEGER(s_max_vertex)[0] - 1);
    size_t n_samples = static_cast<size_t>(INTEGER(s_n_samples)[0]);
    unsigned int seed = static_cast<unsigned int>(INTEGER(s_seed)[0]);
    
    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    
    // Reconstruct MADAG
    madag_t madag;
    madag.source_vertex = source_vertex;
    madag.source_value = y[source_vertex];
    
    int* p_reachable = INTEGER(s_reachable);
    size_t n_reach = LENGTH(s_reachable);
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = static_cast<size_t>(p_reachable[i] - 1);
        madag.reachable_vertices.push_back(v);
        madag.reachable_set.insert(v);
    }
    
    SEXP succ_names = Rf_getAttrib(s_successors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_successors); ++i) {
        const char* name = CHAR(STRING_ELT(succ_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP succ_vec = VECTOR_ELT(s_successors, i);
        std::vector<size_t> successors;
        int* p_succ = INTEGER(succ_vec);
        for (int j = 0; j < LENGTH(succ_vec); ++j) {
            successors.push_back(static_cast<size_t>(p_succ[j] - 1));
        }
        madag.successors[v] = successors;
    }
    
    SEXP pred_names = Rf_getAttrib(s_predecessors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_predecessors); ++i) {
        const char* name = CHAR(STRING_ELT(pred_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP pred_vec = VECTOR_ELT(s_predecessors, i);
        std::vector<size_t> predecessors;
        int* p_pred = INTEGER(pred_vec);
        for (int j = 0; j < LENGTH(pred_vec); ++j) {
            predecessors.push_back(static_cast<size_t>(p_pred[j] - 1));
        }
        madag.predecessors[v] = predecessors;
    }
    
    // Sample trajectories
    std::vector<madag_trajectory_t> trajectories = 
        sample_cell_trajectories(madag, y, max_vertex, n_samples, seed);
    
    return trajectories_to_R(trajectories);
}

SEXP S_trajectory_similarity_matrix(
    SEXP s_trajectories,
    SEXP s_similarity_type
) {
    // Convert trajectories from R list
    std::vector<madag_trajectory_t> trajectories;
    size_t n_traj = LENGTH(s_trajectories);
    
    for (size_t i = 0; i < n_traj; ++i) {
        SEXP traj_r = VECTOR_ELT(s_trajectories, i);
        SEXP vertices_r = VECTOR_ELT(traj_r, 0);  // vertices
        
        madag_trajectory_t traj;
        int* p_v = INTEGER(vertices_r);
        for (int j = 0; j < LENGTH(vertices_r); ++j) {
            traj.vertices.push_back(static_cast<size_t>(p_v[j] - 1));
        }
        trajectories.push_back(traj);
    }
    
    std::string sim_type = CHAR(STRING_ELT(s_similarity_type, 0));
    
    // Compute similarity matrix
    std::vector<std::vector<double>> sim_matrix = 
        compute_trajectory_similarity_matrix(trajectories, sim_type);
    
    // Convert to R matrix
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, n_traj, n_traj));
    double* p_result = REAL(result);
    
    for (size_t i = 0; i < n_traj; ++i) {
        for (size_t j = 0; j < n_traj; ++j) {
            p_result[i + j * n_traj] = sim_matrix[i][j];
        }
    }
    
    UNPROTECT(1);
    return result;
}

SEXP S_identify_bottlenecks(
    SEXP s_madag,
    SEXP s_max_vertex,
    SEXP s_min_fraction
) {
    // Reconstruct minimal MADAG structure
    SEXP s_source = VECTOR_ELT(s_madag, 0);
    SEXP s_reachable = VECTOR_ELT(s_madag, 2);
    SEXP s_successors = VECTOR_ELT(s_madag, 5);
    SEXP s_predecessors = VECTOR_ELT(s_madag, 4);
    SEXP s_path_from = VECTOR_ELT(s_madag, 7);
    
    size_t source_vertex = static_cast<size_t>(INTEGER(s_source)[0] - 1);
    size_t max_vertex = static_cast<size_t>(INTEGER(s_max_vertex)[0] - 1);
    double min_fraction = REAL(s_min_fraction)[0];
    
    // Reconstruct MADAG
    madag_t madag;
    madag.source_vertex = source_vertex;
    
    int* p_reachable = INTEGER(s_reachable);
    size_t n_reach = LENGTH(s_reachable);
    for (size_t i = 0; i < n_reach; ++i) {
        size_t v = static_cast<size_t>(p_reachable[i] - 1);
        madag.reachable_vertices.push_back(v);
        madag.reachable_set.insert(v);
    }
    
    // Reconstruct successors
    SEXP succ_names = Rf_getAttrib(s_successors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_successors); ++i) {
        const char* name = CHAR(STRING_ELT(succ_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP succ_vec = VECTOR_ELT(s_successors, i);
        std::vector<size_t> successors;
        int* p_succ = INTEGER(succ_vec);
        for (int j = 0; j < LENGTH(succ_vec); ++j) {
            successors.push_back(static_cast<size_t>(p_succ[j] - 1));
        }
        madag.successors[v] = successors;
    }
    
    // Reconstruct predecessors
    SEXP pred_names = Rf_getAttrib(s_predecessors, R_NamesSymbol);
    for (int i = 0; i < LENGTH(s_predecessors); ++i) {
        const char* name = CHAR(STRING_ELT(pred_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        
        SEXP pred_vec = VECTOR_ELT(s_predecessors, i);
        std::vector<size_t> predecessors;
        int* p_pred = INTEGER(pred_vec);
        for (int j = 0; j < LENGTH(pred_vec); ++j) {
            predecessors.push_back(static_cast<size_t>(p_pred[j] - 1));
        }
        madag.predecessors[v] = predecessors;
    }
    
    // Reconstruct path counts from source
    SEXP path_from_names = Rf_getAttrib(s_path_from, R_NamesSymbol);
    int* p_path_from = INTEGER(s_path_from);
    for (int i = 0; i < LENGTH(s_path_from); ++i) {
        const char* name = CHAR(STRING_ELT(path_from_names, i));
        size_t v = static_cast<size_t>(atoi(name) - 1);
        madag.path_count_from_source[v] = static_cast<size_t>(p_path_from[i]);
    }
    
    // Compute topological order
    madag.topological_order = compute_topological_order(madag);
    madag.reverse_topological_order = madag.topological_order;
    std::reverse(madag.reverse_topological_order.begin(),
                 madag.reverse_topological_order.end());
    
    // Identify reachable maxima
    for (size_t v : madag.reachable_vertices) {
        if (madag.is_sink(v)) {
            madag.reachable_maxima.push_back(v);
        }
    }
    
    // Identify bottlenecks
    std::vector<size_t> bottlenecks = identify_bottlenecks(madag, max_vertex, min_fraction);
    
    // Convert to R (1-based)
    SEXP result = PROTECT(Rf_allocVector(INTSXP, bottlenecks.size()));
    int* p_result = INTEGER(result);
    for (size_t i = 0; i < bottlenecks.size(); ++i) {
        p_result[i] = static_cast<int>(bottlenecks[i]) + 1;
    }
    
    UNPROTECT(1);
    return result;
}

} // extern "C"
