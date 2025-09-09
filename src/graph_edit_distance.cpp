#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// Undefine conflicting macros after including R headers
#undef length
#undef eval

#include <vector>
#include <cmath>

#include "SEXP_cpp_conversion_utils.hpp"

extern "C" {
    SEXP S_graph_edit_distance(SEXP s_graph1_adj_list,
                               SEXP s_graph1_weights_list,
                               SEXP s_graph2_adj_list,
                               SEXP s_graph2_weights_list,
                               SEXP s_edge_cost,
                               SEXP s_weight_cost_factor);
}

/**
 * @brief Calculate the graph edit distance between two weighted graphs
 *
 * This function computes the graph edit distance between two weighted graphs
 * represented by their adjacency lists and corresponding weight lists.
 *
 * @param graph1_adj_list Adjacency list of the first graph
 * @param graph1_weights Weight list of the first graph
 * @param graph2_adj_list Adjacency list of the second graph
 * @param graph2_weights Weight list of the second graph
 * @param edge_cost Cost of adding or removing an edge
 * @param weight_cost_factor Factor to scale weight differences
 *
 * @return The graph edit distance between the two input graphs
 *
 * @throws std::runtime_error if the graphs have different numbers of vertices
 *
 * @note The adjacency lists and weight lists should have the same shape for each graph
 */
// Main computation function
double graph_edit_distance(const std::vector<std::vector<int>>& graph1_adj_list,
                           const std::vector<std::vector<double>>& graph1_weights,
                           const std::vector<std::vector<int>>& graph2_adj_list,
                           const std::vector<std::vector<double>>& graph2_weights,
                           double edge_cost,
                           double weight_cost_factor) {
    // Input validation
    if (graph1_adj_list.empty() || graph2_adj_list.empty()) {
        Rf_error("Empty graph provided");
    }

    int n = graph1_adj_list.size();
    if (n != graph2_adj_list.size()) {
        Rf_error("Graphs must have the same number of vertices");
    }

    if (graph1_adj_list.size() != graph1_weights.size() ||
        graph2_adj_list.size() != graph2_weights.size()) {
        Rf_error("Adjacency lists and weights must have matching dimensions");
    }

    // Validate vertex indices and weights
    for (int v = 0; v < n; ++v) {
        for (int adj_idx : graph1_adj_list[v]) {
            if (adj_idx < 0 || adj_idx >= n) {
                Rf_error("Invalid vertex index in graph1");
            }
        }
        for (int adj_idx : graph2_adj_list[v]) {
            if (adj_idx < 0 || adj_idx >= n) {
                Rf_error("Invalid vertex index in graph2");
            }
        }
        for (double weight : graph1_weights[v]) {
            if (std::isnan(weight) || std::isinf(weight)) {
                Rf_error("Invalid weight value in graph1");
            }
        }
        for (double weight : graph2_weights[v]) {
            if (std::isnan(weight) || std::isinf(weight)) {
                Rf_error("Invalid weight value in graph2");
            }
        }
    }

    double ged = 0.0;

    for (int v1 = 0; v1 < n; ++v1) {
        for (int v2 = v1 + 1; v2 < n; ++v2) {
            bool edge1_exists = false;
            bool edge2_exists = false;
            double weight1 = 0.0;
            double weight2 = 0.0;

            // Check if edge exists in graph1
            for (size_t i = 0; i < graph1_adj_list[v1].size(); ++i) {
                if (graph1_adj_list[v1][i] == v2) {
                    edge1_exists = true;
                    weight1 = graph1_weights[v1][i];
                    break;
                }
            }

            // Check if edge exists in graph2
            for (size_t i = 0; i < graph2_adj_list[v1].size(); ++i) {
                if (graph2_adj_list[v1][i] == v2) {
                    edge2_exists = true;
                    weight2 = graph2_weights[v1][i];
                    break;
                }
            }

            if (!edge1_exists && !edge2_exists) {
                continue;
            } else if (edge1_exists != edge2_exists) {
                ged += edge_cost;
            } else {
                double weight_diff = std::abs(weight1 - weight2);
                if (std::isnan(weight_diff) || std::isinf(weight_diff)) {
                    Rf_error("Invalid weight difference computation");
                }
                ged += weight_diff * weight_cost_factor;
            }
        }
    }

    return ged;
}

#if 0
double graph_edit_distance(const std::vector<std::vector<int>>& graph1_adj_list,
                           const std::vector<std::vector<double>>& graph1_weights,
                           const std::vector<std::vector<int>>& graph2_adj_list,
                           const std::vector<std::vector<double>>& graph2_weights,
                           double edge_cost,
                           double weight_cost_factor) {
    int n = graph1_adj_list.size();
    if (n != graph2_adj_list.size()) {
        Rf_error("Graphs must have the same number of vertices");
    }

    double ged = 0.0;

    for (int v1 = 0; v1 < n; ++v1) {
        for (int v2 = v1 + 1; v2 < n; ++v2) {
            bool edge1_exists = false;
            bool edge2_exists = false;
            double weight1 = 0.0;
            double weight2 = 0.0;

            // Check if edge exists in graph1
            for (size_t i = 0; i < graph1_adj_list[v1].size(); ++i) {
                if (graph1_adj_list[v1][i] == v2) {
                    edge1_exists = true;
                    weight1 = graph1_weights[v1][i];
                    break;
                }
            }

            // Check if edge exists in graph2
            for (size_t i = 0; i < graph2_adj_list[v1].size(); ++i) {
                if (graph2_adj_list[v1][i] == v2) {
                    edge2_exists = true;
                    weight2 = graph2_weights[v1][i];
                    break;
                }
            }

            if (!edge1_exists && !edge2_exists) {
                // Edge doesn't exist in both graphs, no change needed
                continue;
            } else if (edge1_exists != edge2_exists) {
                // Edge exists in one graph but not the other, count as insertion/deletion
                ged += edge_cost;
            } else {
                // Edge exists in both graphs, calculate weight difference
                ged += std::abs(weight1 - weight2) * weight_cost_factor;
            }
        }
    }

    return ged;
}
#endif


/**
 * @brief R-callable wrapper for the graph edit distance calculation
 *
 * This function serves as a wrapper to call the C++ implementation of graph edit distance
 * calculation from R. It handles the conversion of R data structures to C++ and vice versa.
 *
 * @param s_graph1_adj_list R list of integer vectors representing the adjacency list of the first graph
 * @param s_graph1_weights_list R list of double vectors representing the weights of the first graph
 * @param s_graph2_adj_list R list of integer vectors representing the adjacency list of the second graph
 * @param s_graph2_weights_list R list of double vectors representing the weights of the second graph
 * @param s_edge_cost R numeric scalar representing the cost of adding or removing an edge
 * @param s_weight_cost_factor R numeric scalar representing the factor to scale weight differences
 *
 * @return SEXP (R numeric) containing the calculated graph edit distance
 */
SEXP S_graph_edit_distance(SEXP s_graph1_adj_list,
                           SEXP s_graph1_weights_list,
                           SEXP s_graph2_adj_list,
                           SEXP s_graph2_weights_list,
                           SEXP s_edge_cost,
                           SEXP s_weight_cost_factor) {
    if (TYPEOF(s_edge_cost) != REALSXP || TYPEOF(s_weight_cost_factor) != REALSXP) {
        Rf_error("Cost parameters must be numeric");
    }

    if (Rf_length(s_edge_cost) != 1 || Rf_length(s_weight_cost_factor) != 1) {
        Rf_error("Cost parameters must be scalar values");
    }

    try {
        std::vector<std::vector<int>> graph1_adj_vv        = convert_adj_list_from_R(s_graph1_adj_list);
        std::vector<std::vector<double>> graph1_weights_vv = convert_weight_list_from_R(s_graph1_weights_list);
        std::vector<std::vector<int>> graph2_adj_vv        = convert_adj_list_from_R(s_graph2_adj_list);
        std::vector<std::vector<double>> graph2_weights_vv = convert_weight_list_from_R(s_graph2_weights_list);

        // Additional safety checks
        for (size_t i = 0; i < graph1_adj_vv.size(); ++i) {
            if (graph1_adj_vv[i].size() != graph1_weights_vv[i].size()) {
                Rf_error("Mismatch between adjacency list and weights dimensions in graph1");
            }
        }

        for (size_t i = 0; i < graph2_adj_vv.size(); ++i) {
            if (graph2_adj_vv[i].size() != graph2_weights_vv[i].size()) {
                Rf_error("Mismatch between adjacency list and weights dimensions in graph2");
            }
        }

        double edge_cost = REAL(s_edge_cost)[0];
        double weight_cost_factor = REAL(s_weight_cost_factor)[0];

        if (std::isnan(edge_cost) || std::isnan(weight_cost_factor)) {
            Rf_error("Cost parameters cannot be NA/NaN");
        }

        double res = graph_edit_distance(graph1_adj_vv,
                                       graph1_weights_vv,
                                       graph2_adj_vv,
                                       graph2_weights_vv,
                                       edge_cost,
                                       weight_cost_factor);

        if (std::isnan(res) || std::isinf(res)) {
            Rf_error("Computation resulted in invalid value (NA/NaN/Inf)");
        }

        SEXP s_res = PROTECT(Rf_allocVector(REALSXP, 1));
        REAL(s_res)[0] = res;
        UNPROTECT(1);
        return s_res;

    } catch (const std::exception& e) {
        Rf_error("C++ error: %s", e.what());
    } catch (...) {
        Rf_error("Unknown C++ error occurred");
    }
}



#if 0
SEXP S_graph_edit_distance(SEXP s_graph1_adj_list,
                           SEXP s_graph1_weights_list,
                           SEXP s_graph2_adj_list,
                           SEXP s_graph2_weights_list,
                           SEXP s_edge_cost,
                           SEXP s_weight_cost_factor) {

    std::vector<std::vector<int>> graph1_adj_vv        = convert_adj_list_from_R(s_graph1_adj_list);
    std::vector<std::vector<double>> graph1_weights_vv = convert_weight_list_from_R(s_graph1_weights_list);

    std::vector<std::vector<int>> graph2_adj_vv        = convert_adj_list_from_R(s_graph2_adj_list);
    std::vector<std::vector<double>> graph2_weights_vv = convert_weight_list_from_R(s_graph2_weights_list);

    double edge_cost = REAL(s_edge_cost)[0];
    double weight_cost_factor = REAL(s_weight_cost_factor)[0];

    double res = graph_edit_distance(graph1_adj_vv,
                                     graph1_weights_vv,
                                     graph2_adj_vv,
                                     graph2_weights_vv,
                                     edge_cost,
                                     weight_cost_factor);

    // Converting the double result to an SEXP object
    SEXP s_res = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(s_res)[0] = res;
    UNPROTECT(1);

    return s_res;
}
#endif
