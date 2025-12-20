/**
 * @file se_tree_r.cpp
 * @brief SEXP interface for SE tree computation
 */

#include <R.h>
#include <Rinternals.h>
#include "se_tree.hpp"
#include "set_wgraph.hpp"
#include "gfc.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>

// Forward declarations
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);
void compute_se_tree_hr_support(se_tree_t& tree, bool verbose);

/**
 * @brief Convert se_node_t to R list
 */
static SEXP se_node_to_R(const se_node_t& node) {
    const int n_components = 7;
    SEXP s_node = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    int idx = 0;

    // vertex (1-based)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("vertex"));
    SET_VECTOR_ELT(s_node, idx++, Rf_ScalarInteger(static_cast<int>(node.vertex) + 1));

    // is.maximum
    SET_STRING_ELT(s_names, idx, Rf_mkChar("is.maximum"));
    SET_VECTOR_ELT(s_node, idx++, Rf_ScalarLogical(node.is_maximum ? TRUE : FALSE));

    // is.spurious
    SET_STRING_ELT(s_names, idx, Rf_mkChar("is.spurious"));
    SET_VECTOR_ELT(s_node, idx++, Rf_ScalarLogical(node.is_spurious ? TRUE : FALSE));

    // parent (1-based, NA if root)
    SET_STRING_ELT(s_names, idx, Rf_mkChar("parent"));
    if (node.parent == SIZE_MAX) {
        SET_VECTOR_ELT(s_node, idx++, Rf_ScalarInteger(NA_INTEGER));
    } else {
        SET_VECTOR_ELT(s_node, idx++, Rf_ScalarInteger(static_cast<int>(node.parent) + 1));
    }

    // children (1-based)
    const int n_children = static_cast<int>(node.children.size());
    SEXP s_children = PROTECT(Rf_allocVector(INTSXP, n_children));
    int* p_children = INTEGER(s_children);
    for (int i = 0; i < n_children; ++i) {
        p_children[i] = static_cast<int>(node.children[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("children"));
    SET_VECTOR_ELT(s_node, idx++, s_children);
    UNPROTECT(1);

    // basin.vertices (1-based)
    const int n_basin = static_cast<int>(node.basin_vertices.size());
    SEXP s_basin = PROTECT(Rf_allocVector(INTSXP, n_basin));
    int* p_basin = INTEGER(s_basin);
    for (int i = 0; i < n_basin; ++i) {
        p_basin[i] = static_cast<int>(node.basin_vertices[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("basin.vertices"));
    SET_VECTOR_ELT(s_node, idx++, s_basin);
    UNPROTECT(1);

    // basin.boundary (1-based)
    const int n_boundary = static_cast<int>(node.basin_boundary.size());
    SEXP s_boundary = PROTECT(Rf_allocVector(INTSXP, n_boundary));
    int* p_boundary = INTEGER(s_boundary);
    for (int i = 0; i < n_boundary; ++i) {
        p_boundary[i] = static_cast<int>(node.basin_boundary[i]) + 1;
    }
    SET_STRING_ELT(s_names, idx, Rf_mkChar("basin.boundary"));
    SET_VECTOR_ELT(s_node, idx++, s_boundary);
    UNPROTECT(1);

    Rf_setAttrib(s_node, R_NamesSymbol, s_names);
    UNPROTECT(2);  // s_node, s_names

    return s_node;
}

extern "C" SEXP S_build_se_tree(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_root_vertex,
    SEXP s_spurious_min,
    SEXP s_spurious_max,
    SEXP s_basins,          // Pre-computed basins from compute.gfc.basins()
    SEXP s_compute_hr_support,
    SEXP s_verbose
) {
    // Convert graph structure
    auto adj_list = convert_adj_list_from_R(s_adj_list);
    auto weight_list = convert_weight_list_from_R(s_weight_list);
    set_wgraph_t graph(adj_list, weight_list);

    // Convert y vector
    const int n = LENGTH(s_y);
    std::vector<double> y(n);
    const double* p_y = REAL(s_y);
    for (int i = 0; i < n; ++i) {
        y[i] = p_y[i];
    }

    // Root vertex (convert from 1-based to 0-based)
    size_t root_vertex = static_cast<size_t>(Rf_asInteger(s_root_vertex) - 1);

    // Convert spurious sets (1-based to 0-based)
    std::unordered_set<size_t> spurious_min;
    const int* p_smin = INTEGER(s_spurious_min);
    for (int i = 0; i < LENGTH(s_spurious_min); ++i) {
        spurious_min.insert(static_cast<size_t>(p_smin[i] - 1));
    }

    std::unordered_set<size_t> spurious_max;
    const int* p_smax = INTEGER(s_spurious_max);
    for (int i = 0; i < LENGTH(s_spurious_max); ++i) {
        spurious_max.insert(static_cast<size_t>(p_smax[i] - 1));
    }

    // Convert basins from R list to C++ map
    std::unordered_map<size_t, bbasin_t> basins;
    const int n_basins = LENGTH(s_basins);
    for (int i = 0; i < n_basins; ++i) {
        SEXP s_basin = VECTOR_ELT(s_basins, i);

        bbasin_t basin;
        basin.extremum_vertex = static_cast<size_t>(
            Rf_asInteger(VECTOR_ELT(s_basin, 0)) - 1);  // 1-based to 0-based
        basin.value = Rf_asReal(VECTOR_ELT(s_basin, 1));
        basin.is_maximum = Rf_asLogical(VECTOR_ELT(s_basin, 2));

        // vertices
        SEXP s_verts = VECTOR_ELT(s_basin, 3);
        const int* p_verts = INTEGER(s_verts);
        for (int j = 0; j < LENGTH(s_verts); ++j) {
            basin.vertices.push_back(static_cast<size_t>(p_verts[j] - 1));
        }

        // boundary
        SEXP s_bd = VECTOR_ELT(s_basin, 4);
        const int* p_bd = INTEGER(s_bd);
        for (int j = 0; j < LENGTH(s_bd); ++j) {
            basin.boundary.push_back(static_cast<size_t>(p_bd[j] - 1));
        }

        basins[basin.extremum_vertex] = basin;
    }

    bool compute_hr_support = Rf_asLogical(s_compute_hr_support);
    bool verbose = Rf_asLogical(s_verbose);

    // Build SE tree
    se_tree_t tree = graph.build_se_tree(
        root_vertex, basins, spurious_min, spurious_max, verbose
    );

    // Compute HR support if requested
    if (compute_hr_support) {
        compute_se_tree_hr_support(tree, verbose);
    }

    // Convert to R list
    const int n_result_components = 8;
    SEXP s_result = PROTECT(Rf_allocVector(VECSXP, n_result_components));
    SEXP s_result_names = PROTECT(Rf_allocVector(STRSXP, n_result_components));

    int idx = 0;

    // root.vertex (1-based)
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("root.vertex"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarInteger(static_cast<int>(tree.root_vertex) + 1));

    // root.is.maximum
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("root.is.maximum"));
    SET_VECTOR_ELT(s_result, idx++, Rf_ScalarLogical(tree.root_is_maximum ? TRUE : FALSE));

    // classification
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("classification"));
    const char* class_str =
        tree.classification == se_tree_class_t::BOTH_TYPES ? "both_types" :
        tree.classification == se_tree_class_t::ONLY_MIN ? "only_min" :
        tree.classification == se_tree_class_t::ONLY_MAX ? "only_max" : "no_terminals";
    SET_VECTOR_ELT(s_result, idx++, Rf_mkString(class_str));

    // nodes (list of node structures)
    const int n_nodes = static_cast<int>(tree.nodes.size());
    SEXP s_nodes = PROTECT(Rf_allocVector(VECSXP, n_nodes));
    SEXP s_node_names = PROTECT(Rf_allocVector(STRSXP, n_nodes));
    int node_idx = 0;
    for (const auto& [vertex, node] : tree.nodes) {
        char name_buf[32];
        snprintf(name_buf, sizeof(name_buf), "%s_%d",
                 node.is_maximum ? "max" : "min",
                 static_cast<int>(vertex) + 1);
        SET_STRING_ELT(s_node_names, node_idx, Rf_mkChar(name_buf));

        SEXP s_node = PROTECT(se_node_to_R(node));
        SET_VECTOR_ELT(s_nodes, node_idx, s_node);
        UNPROTECT(1);
        ++node_idx;
    }
    Rf_setAttrib(s_nodes, R_NamesSymbol, s_node_names);
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("nodes"));
    SET_VECTOR_ELT(s_result, idx++, s_nodes);
    UNPROTECT(2);  // s_nodes, s_node_names

    // ns.min.terminals (1-based)
    const int n_ns_min = static_cast<int>(tree.ns_min_terminals.size());
    SEXP s_ns_min = PROTECT(Rf_allocVector(INTSXP, n_ns_min));
    int* p_ns_min = INTEGER(s_ns_min);
    for (int i = 0; i < n_ns_min; ++i) {
        p_ns_min[i] = static_cast<int>(tree.ns_min_terminals[i]) + 1;
    }
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("ns.min.terminals"));
    SET_VECTOR_ELT(s_result, idx++, s_ns_min);
    UNPROTECT(1);

    // ns.max.terminals (1-based)
    const int n_ns_max = static_cast<int>(tree.ns_max_terminals.size());
    SEXP s_ns_max = PROTECT(Rf_allocVector(INTSXP, n_ns_max));
    int* p_ns_max = INTEGER(s_ns_max);
    for (int i = 0; i < n_ns_max; ++i) {
        p_ns_max[i] = static_cast<int>(tree.ns_max_terminals[i]) + 1;
    }
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("ns.max.terminals"));
    SET_VECTOR_ELT(s_result, idx++, s_ns_max);
    UNPROTECT(1);

    // hr.support.vertices (1-based)
    const int n_hr_verts = static_cast<int>(tree.hr_support_vertices.size());
    SEXP s_hr_verts = PROTECT(Rf_allocVector(INTSXP, n_hr_verts));
    int* p_hr_verts = INTEGER(s_hr_verts);
    for (int i = 0; i < n_hr_verts; ++i) {
        p_hr_verts[i] = static_cast<int>(tree.hr_support_vertices[i]) + 1;
    }
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("hr.support.vertices"));
    SET_VECTOR_ELT(s_result, idx++, s_hr_verts);
    UNPROTECT(1);

    // hr.support.boundary (1-based)
    const int n_hr_bd = static_cast<int>(tree.hr_support_boundary.size());
    SEXP s_hr_bd = PROTECT(Rf_allocVector(INTSXP, n_hr_bd));
    int* p_hr_bd = INTEGER(s_hr_bd);
    for (int i = 0; i < n_hr_bd; ++i) {
        p_hr_bd[i] = static_cast<int>(tree.hr_support_boundary[i]) + 1;
    }
    SET_STRING_ELT(s_result_names, idx, Rf_mkChar("hr.support.boundary"));
    SET_VECTOR_ELT(s_result, idx++, s_hr_bd);
    UNPROTECT(1);

    Rf_setAttrib(s_result, R_NamesSymbol, s_result_names);

    // Add class
    SEXP s_class = PROTECT(Rf_mkString("se_tree"));
    Rf_setAttrib(s_result, R_ClassSymbol, s_class);

    UNPROTECT(3);  // s_result, s_result_names, s_class
    return s_result;
}
