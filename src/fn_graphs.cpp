#include "set_wgraph.hpp"
#include "error_utils.h"   // REPORT_ERROR()

/**
 * @brief Constructs a new graph with weights modified based on function values
 *
 * @details
 * This function creates a new graph with the same connectivity as the original,
 * but with weights modified according to a specified weighting function that
 * considers both the original weight and the difference in function values.
 * Optionally, edges with weights exceeding a threshold can be pruned.
 *
 * @param function_values Vector of function values at each vertex
 * @param weight_type Type of weight modification function to apply:
 *   0: Inverse relationship: w_new = w_old / (|f(i) - f(j)| + epsilon)
 *   1: Direct relationship: w_new = w_old * |f(i) - f(j)|
 *   2: Exponential decay: w_new = w_old * exp(-lambda * |f(i) - f(j)|)
 *   3: Power law: w_new = w_old * |f(i) - f(j)|^(-alpha)
 *   4: Sigmoid: w_new = w_old * 1/(1 + exp(beta * (|f(i) - f(j)| - tau)))
 *   5: L_p embedding: w_new = (w_old^p + alpha*|f(i) - f(j)|^q)^(1/r)
 * @param epsilon Small constant to avoid division by zero (for weight_type 0)
 * @param lambda Decay rate parameter (for weight_type 2)
 * @param alpha Power law exponent (for weight_type 3) or scaling factor (for weight_type 5)
 * @param beta Sigmoid steepness parameter (for weight_type 4)
 * @param tau Sigmoid threshold parameter (for weight_type 4)
 * @param p Power for feature distance term (for weight_type 5)
 * @param q Power for function difference term (for weight_type 5)
 * @param r Power for the overall normalization (for weight_type 5)
 * @param normalize Whether to normalize the weights after modification
 * @param weight_thld Threshold for pruning edges: edges with weight > weight_thld will be pruned
 *                   Set to a negative value to disable pruning
 * @return A new graph with modified weights
 */
set_wgraph_t set_wgraph_t::construct_function_aware_graph(
    const std::vector<double>& function_values,
    int weight_type,
    double epsilon,
    double lambda,
    double alpha,
    double beta,
    double tau,
    double p,
    double q,
    double r,
    bool normalize,
    double weight_thld
	) const {

    // Ensure we have function values for every vertex
    if (function_values.size() != adjacency_list.size()) {
        REPORT_ERROR("Function values vector must have the same length as the number of vertices");
    }

    // Create a new graph with the same structure
    set_wgraph_t new_graph(adjacency_list.size());

    // Copy over the properties that should be preserved
    new_graph.graph_diameter = -1;     // This needs to be recomputed
    new_graph.max_packing_radius = -1; // This needs to be recomputed

    // Process each vertex and collect normalization data first if needed
    std::vector<std::vector<std::tuple<size_t, size_t, double>>> edge_weights_by_vertex;
    if (normalize) {
        edge_weights_by_vertex.resize(adjacency_list.size());
    }

    // First pass: calculate all new weights
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        double weight_sum = 0.0;

        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;

            // Only process each edge once (for undirected graph)
            if (i < j) {
                double original_weight = edge.weight;
                double f_diff = std::abs(function_values[i] - function_values[j]);

                // Apply the selected weight modification function
                double new_weight;
                switch (weight_type) {
                case 0: // Inverse relationship
                    new_weight = original_weight / (f_diff + epsilon);
                    break;

                case 1: // Direct relationship
                    new_weight = original_weight * f_diff;
                    break;

                case 2: // Exponential decay
                    new_weight = original_weight * std::exp(-lambda * f_diff);
                    break;

                case 3: // Power law
                    new_weight = original_weight * std::pow(f_diff + epsilon, -alpha);
                    break;

                case 4: // Sigmoid
                    new_weight = original_weight * (1.0 / (1.0 + std::exp(beta * (f_diff - tau))));
                    break;

                case 5: // L_p embedding
                    {
                        // Small epsilon to avoid numerical issues with small powers
                        const double eps = 1e-10;
                        double feat_term = std::pow(original_weight + eps, p);
                        double func_term = alpha * std::pow(f_diff + eps, q);
                        new_weight = std::pow(feat_term + func_term, 1.0/r);
                    }
                    break;

                default:
                    REPORT_ERROR("Unsupported weight_type");
                }

                // Skip edges that exceed the weight threshold (if threshold is provided)
                if (weight_thld > 0 && new_weight > weight_thld) {
                    continue;
                }

                // For normalization, collect the edge and its weight
                if (normalize) {
                    edge_weights_by_vertex[i].push_back(std::make_tuple(i, j, new_weight));
                    edge_weights_by_vertex[j].push_back(std::make_tuple(j, i, new_weight));
                    weight_sum += new_weight;
                } else {
                    // If not normalizing, directly add the edge
                    new_graph.add_edge(i, j, new_weight);
                }
            }
        }

        // If normalizing, process all collected edges for this vertex
        if (normalize && weight_sum > 0) {
            for (const auto& [src, dst, weight] : edge_weights_by_vertex[i]) {
                // Only add edges from i to higher indices to avoid duplicates
                if (src == i && dst > i) {
                    double normalized_weight = weight / weight_sum;
                    new_graph.add_edge(src, dst, normalized_weight);
                }
            }
        }
    }

    return new_graph;
}


/**
 * @brief Returns all modified edge weights for different function-aware weighting schemes
 *
 * @details
 * This function calculates how different weighting schemes would modify the edge weights
 * based on function values, and returns the complete collection of weights for each scheme.
 * This allows for detailed distribution analysis in R to identify modes, thresholds,
 * or other patterns in the weight distribution.
 * The original, unmodified weights are included as the first element in the result.
 *
 * @param function_values Vector of function values at each vertex
 * @param weight_types Vector of weight_type values to analyze
 * @param epsilon Small constant to avoid division by zero (for weight_type 0)
 * @param lambda Decay rate parameter (for weight_type 2)
 * @param alpha Power law exponent (for weight_type 3) or scaling factor (for weight_type 5)
 * @param beta Sigmoid steepness parameter (for weight_type 4)
 * @param tau Sigmoid threshold parameter (for weight_type 4)
 * @param p Power for feature distance term (for weight_type 5)
 * @param q Power for function difference term (for weight_type 5)
 * @param r Power for the overall normalization (for weight_type 5)
 * @return A vector of weight vectors, where each inner vector contains all calculated
 *         weights for the corresponding weight_type, with the first vector containing
 *         the original, unmodified weights
 */
std::vector<std::vector<double>> set_wgraph_t::analyze_function_aware_weights(
    const std::vector<double>& function_values,
    const std::vector<int>& weight_types,
    double epsilon,
    double lambda,
    double alpha,
    double beta,
    double tau,
    double p,
    double q,
    double r
	) const {

    // Calculate exact number of edges in the graph
    size_t edge_count = 0;
    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        for (const auto& edge : adjacency_list[i]) {
            if (i < edge.vertex) { // Count each edge only once
                edge_count++;
            }
        }
    }

    // Result vector containing a weight vector for each weight type, plus one for original weights
    std::vector<std::vector<double>> results(weight_types.size() + 1);

    // Pre-allocate vectors with exact size
    for (auto& weights : results) {
        weights.resize(edge_count);
    }

    // First collect the original weights
    size_t orig_counter = 0;
    std::vector<double>& original_weights = results[0];

    for (size_t i = 0; i < adjacency_list.size(); ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            if (i < j) { // Process each edge only once
                original_weights[orig_counter++] = edge.weight;
            }
        }
    }

    // For each weight type
    for (size_t type_idx = 0; type_idx < weight_types.size(); ++type_idx) {
        int weight_type = weight_types[type_idx];
        std::vector<double>& modified_weights = results[type_idx + 1]; // +1 to account for original weights
        size_t counter = 0;

        // Collect all modified weights
        for (size_t i = 0; i < adjacency_list.size(); ++i) {
            for (const auto& edge : adjacency_list[i]) {
                size_t j = edge.vertex;
                if (i < j) { // Process each edge only once
                    double original_weight = edge.weight;
                    double f_diff = std::abs(function_values[i] - function_values[j]);

                    // Apply the selected weight modification function
                    double new_weight;
                    switch (weight_type) {
                        case 0: // Inverse relationship
                            new_weight = original_weight / (f_diff + epsilon);
                            break;
                        case 1: // Direct relationship
                            new_weight = original_weight * f_diff;
                            break;
                        case 2: // Exponential decay
                            new_weight = original_weight * std::exp(-lambda * f_diff);
                            break;
                        case 3: // Power law
                            new_weight = original_weight * std::pow(f_diff + epsilon, -alpha);
                            break;
                        case 4: // Sigmoid
                            new_weight = original_weight * (1.0 / (1.0 + std::exp(beta * (f_diff - tau))));
                            break;
                        case 5: // L_p embedding
                            {
                                // Small epsilon to avoid numerical issues with small powers
                                const double eps = 1e-10;
                                double feat_term = std::pow(original_weight + eps, p);
                                double func_term = alpha * std::pow(f_diff + eps, q);
                                new_weight = std::pow(feat_term + func_term, 1.0/r);
                            }
                            break;
                        default:
                            REPORT_ERROR("Unsupported weight_type");
                    }

                    modified_weights[counter++] = new_weight;
                }
            }
        }
    }

    return results;
}
