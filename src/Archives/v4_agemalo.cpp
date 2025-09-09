/**
 * @brief Implementation of Adaptive GEMALO algorithm
 */

#include <execution>
#include <thread>      // For std::thread
#include <random>      // For std::mt19937, std::random_device, std::gamma_distribution
#include <filesystem>  // for debugging
#include <fstream>

#include "agemalo.hpp"
#include "kernels.h"
#include "error_utils.h"         // For REPORT_ERROR()
#include "predictive_errors.hpp" // For struct bb_cri_t
#include "ext_ulm_priority_queue.hpp"

double calculate_directional_threshold(
    double bw,
    double bw_min,
    double graph_diameter,
    size_t min_path_size,
    double dir_thld_factor = 1.0/3.0);


/**
 * @brief Implements the Adaptive Graph Geodesic Model-Averaged LOcal Linear regression (Adaptive-GEMALO) algorithm
 *
 * @details This function implements an conditional expectation estimation of a
 * random variable defined over vertices of a grid graph using local geodesic
 * paths linear models with model averaging and adaptive bandwidth selection.
 *
 * @param adj_list Adjacency list representation of the graph, where adj_list[i] contains indices of vertices adjacent to vertex i
 * @param weight_list Corresponding edge weights, where weight_list[i][j] is the weight of edge from vertex i to adj_list[i][j]
 * @param y Response values at each vertex in the original graph
 * @param start_vertex Index of the vertex to start graph traversal
 * @param min_path_size Minimum number of vertices required in valid paths
 * @param n_grid_vertices Number of vertices in the uniform grid representation
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param snap_tolerance Maximum distance for snapping original vertices to grid vertices
 * @param dist_normalization_factor Factor for normalizing distances in the graph
 * @param kernel_type Type of kernel function used for weight calculation (e.g., Gaussian, triangular)
 * @param tolerance Convergence tolerance for linear model fitting
 * @param n_bb Number of Bayesian bootstrap iterations (0 for no bootstrap)
 * @param cri_probability Probability level for confidence intervals in bootstrap
 * @param n_perms Number of permutation test iterations (0 for no permutation testing)
 * @param blending_coef Number between 0 and 1 allowing for smooth interpolation between position-based weights and mean-error times position-based weights using the formula effective_weight = x.weight * pow(x.mean_error, blending_coef);
 * @param verbose Whether to print progress information
 *
 * @return adaptive_uggmalo_result_t structure containing:
 *         - graph_diameter: Computed diameter of the input graph
 *         - grid_opt_bw: Optimal bandwidth for each grid vertex
 *         - predictions: Model-averaged predictions for original graph vertices
 *         - grid_predictions: Model-averaged predictions for grid vertices
 *         - Bootstrap results (if n_bb > 0):
 *           - bb_predictions: Bootstrap replication predictions
 *           - cri_lower: Lower confidence interval bounds
 *           - cri_upper: Upper confidence interval bounds
 *         - Permutation test results (if n_perms > 0):
 *           - null_predictions: Predictions under null hypothesis
 *           - permutation_tests: Statistical test results including p-values and effect sizes
 *
 */
agemalo_result_t agemalo(
    const uniform_grid_graph_t& grid_graph,
    const std::vector<double>& y,
	// geodesic parameter
    size_t min_path_size,
	// bw parameters
	size_t n_bws,
    double max_bw_factor,
	// kernel parameters
    double dist_normalization_factor,
    size_t kernel_type,
	// model parameters
    double model_tolerance,
	double model_blending_coef,
	// Bayesian bootstrap parameters
    size_t n_bb,
    double cri_probability,
	// permutation parameters
    size_t n_perms,
    // verbose
    bool verbose
    ) {

    // std::for_each flag
    // auto exec = std::execution::seq;
    // auto exec = std::execution::par_unseq;

    // Initialize result structure
    agemalo_result_t result;
    result.graph_diameter = grid_graph.graph_diameter;
	result.max_packing_radius = grid_graph.max_packing_radius;

	// Initialize vectors to store predictions for all vertices
    const size_t n_vertices = grid_graph.adjacency_list.size();
    result.predictions.resize(n_vertices, INFINITY);
    result.errors.resize(n_vertices, INFINITY);
    result.scale.resize(n_vertices, INFINITY);

	bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});
    double max_bw = max_bw_factor * result.graph_diameter;

	initialize_kernel(kernel_type, 1.0);

    if (verbose) {
        Rprintf("In agemalo()\n");
		Rprintf("n_vertices(grid_graph): %zu\n", grid_graph.adjacency_list.size());
		// packing parameters
		Rprintf("graph_diameter: %f\n", result.graph_diameter);
		Rprintf("max_packing_radius: %f\n", result.max_packing_radius);
		Rprintf("grid_graph.grid_vertices.size(): %zu\n", grid_graph.grid_vertices.size());
		// bw parameters
        Rprintf("max_bw_factor: %f\n", max_bw_factor);
        Rprintf("max_bw: %f\n", max_bw);
		Rprintf("n_bws: %zu\n", n_bws);
    }

	// Before create_grid_vertex_models lambda
	std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data/";
	// Ensure directory exists
	std::filesystem::create_directories(debug_dir);

	std::ofstream grid_vertices_file(debug_dir + "_grid_vertices.csv");
	for (const auto& vertex : grid_graph.grid_vertices) {
		grid_vertices_file << (vertex + 1) << "\n";
	}

	// Export full graph structure
	{
		std::ofstream graph_file(debug_dir + "full_graph.csv");
		graph_file << "from,to,weight\n";

		for (size_t i = 0; i < grid_graph.adjacency_list.size(); ++i) {
			for (const auto& edge : grid_graph.adjacency_list[i]) {
				graph_file << (i + 1) << "," << (edge.vertex + 1) << "," << edge.weight << "\n";
			}
		}
	}

	// Lambda for creating models at grid vertices - supposet to be thread safe
    // Returns: Min Heap: Queue of models sorted in the ascending order by mean error
    auto create_grid_vertex_models = [&](
        size_t grid_vertex,
        const std::vector<double>& y,
        const std::optional<std::vector<double>>& weights
        ) {

		// Fit local models over classes of geodesic paths defined over a range of bandwidth balls within the given graph
        ext_ulm_priority_queue model_queue; // a min heap priority queue of models sorted in the increasing order of their mean predition error

		double radius = 0.2 * result.graph_diameter; // <--- this was set to this constant only for testing/debugging purposes
		shortest_paths_t shortest_paths = grid_graph.find_graph_paths_within_radius(grid_vertex, radius);

		// Create a vector to store the processed ray data
		std::vector<gray_xyw_t> ray_data;
		ray_data.reserve(shortest_paths.paths.size());

		// Process each path to extract x/y/w data
		for (auto& path : shortest_paths.paths) {
			ray_data.push_back(grid_graph.get_xyw_along_path(y, path, dist_normalization_factor));
		}

        return model_queue;
    }; // END OF create_grid_vertex_models()

    // weight/prediction/error/mean_error/be struct needed for mode averaging and local scale estimation; we use it to record weight/prediction/error/bw of the given vertex in each model where the vertext is in the support of the model
    struct wpe_t {
        double weight;
        double prediction;
        double error;
        double mean_error;
        double bw;

        // Constructor needed for emplace_back(x,y,z)
        wpe_t(double w, double p, double e, double me, double bw)
            : weight(w), prediction(p), error(e), mean_error(me), bw(bw) {}
    };

    struct grid_wpe_t {
        double weight;
        double prediction;
        double mean_error;
        double bw;

        // Constructor needed for emplace_back(x,y,z)
        grid_wpe_t(double w, double p, double me, double bw)
            : weight(w), prediction(p), mean_error(me), bw(bw) {}
    };


    auto average_models = [&model_blending_coef,&n_vertices,&verbose,&grid_graph](
        std::unordered_map<size_t, ext_ulm_priority_queue>& grid_vertex_models_map,
        const std::vector<double>& y,
        std::vector<double>& predictions,
        std::vector<double>& errors,
        std::vector<double>& scale,
        std::optional<std::unordered_map<size_t, double>>& grid_predictions_map
        ) {

        bool process_errors = errors.size() == predictions.size();
        bool process_scale  = scale.size() == predictions.size();

        std::vector<std::vector<wpe_t>> wpe(n_vertices); // wpe[i] stores a vector of {weight, prediction, error, mean_error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions
        std::unordered_map<size_t, std::vector<grid_wpe_t>> grid_wpe_map;

        // Rprintf("In average_models lambda grid_vertex_models_map.size(): %zu\n", grid_vertex_models_map.size());

        for (auto& [i, models_queue] : grid_vertex_models_map) {

            // auto model = models_queue.top();
            // Rprintf("model.grid_vertices.size(): %zu\n", model.grid_vertices.size());
            // if (models_queue.empty()) {
            //     Rprintf("models_queue is empty for i: %zu\n",i);
            // }

            while (!models_queue.empty()) {
                auto model = models_queue.top(); // ext_ulm_t object
                // Store weight, prediction and error for the given vertex from the given 'model' in wpe[ model.vertices[i] ]
                for (size_t i = 0; i < model.vertices.size(); ++i) {
                    wpe[ model.vertices[i] ].emplace_back(model.w_path[i], model.predictions[i], model.errors[i], model.mean_error, model.bw);
                }

                if (grid_predictions_map) {
                    for (size_t i = 0; i < model.grid_vertices.size(); ++i) {
                        size_t grid_idx = model.grid_vertices[i];
                        grid_wpe_map[grid_idx].emplace_back(model.grid_w_path[i], model.grid_predictions[i], model.mean_error, model.bw);
                    }
                }

                models_queue.pop();
            }
        }

        // Model averaging over the vertices of the original graph
        for (size_t i = 0; i < n_vertices; i++ ) {

            double prediction_sum = 0.0;
            double weight_sum     = 0.0;
            double error_sum      = 0.0;
            double local_scale    = 0.0; //INFINITY;
            double effective_weight;

            for (const auto& x : wpe[i]) {
                if (model_blending_coef == 0.0) {
                    // Pure position weight
                    effective_weight = x.weight;
                }
                else if (model_blending_coef == 1.0) {
                    // Full mean error influence
                    effective_weight = x.mean_error * x.weight;
                }
                else {
                    // Smooth interpolation between the two approaches
                    // Method 1: Linear interpolation of weights
                    effective_weight = (1.0 - model_blending_coef) * x.weight + model_blending_coef * (x.mean_error * x.weight);
                    // Method 2: Power-based scaling
                    //effective_weight = x.weight * pow(x.mean_error, model_blending_coef);
                }

                prediction_sum += effective_weight * x.prediction;
                weight_sum     += effective_weight;

                if (process_errors) error_sum += effective_weight * x.error;
                // if (process_scale && x.bw < local_scale) local_scale = x.bw; // min bw method
                if (process_scale) local_scale += effective_weight * x.bw; // weighted mean bw method
            }

            if (weight_sum > 0) {
                if (weight_sum > 1e-10) {
                    predictions[i] = prediction_sum / weight_sum;
                    if (process_errors)
                        errors[i] = error_sum / weight_sum;
                    if (process_scale)
                        scale[i] = local_scale / weight_sum;
                } else {
                    if (verbose) {
                        REPORT_WARNING("Very small weight sum encountered for vertex %zu. Setting predictions[i] to y[i]\n", i);
                    }
                    predictions[i] = y[i];  // Fall back to original value
                }
            } else {
                REPORT_ERROR("Weight sum = 0 for vertex %zu in predictions\n", i);
                //predictions[i] = y[i];
            }
        }

        if (grid_predictions_map) {

            //Rprintf("grid_wpe_map.size(): %zu\n", grid_wpe_map.size());

            // Model averaging over the grid vertices
            for (auto& [i, wpe_entries] : grid_wpe_map) {

                //Rprintf("i: %zu\twpe_entries.size(): %zu\n", i, wpe_entries.size());

                auto x = wpe_entries[0];
                (*grid_predictions_map)[i] = x.prediction;

                if (wpe_entries.size() < 2) {
                    continue;
                }

                double prediction_sum = 0.0;
                double weight_sum = 0.0;
                double effective_weight;

                for (const auto& x : wpe_entries) {

                    if (std::isnan(x.weight) || (model_blending_coef > 0 && std::isnan(x.mean_error))) {
                        continue;
                    }

                    if (model_blending_coef == 0.0) {
                        // Pure position weight
                        effective_weight = x.weight;
                    }
                    else if (model_blending_coef == 1.0) {
                        // Full mean error influence
                        effective_weight = x.mean_error * x.weight;
                    }
                    else {
                        // Smooth interpolation between the two approaches
                        // Method 1: Linear interpolation of weights
                        effective_weight = (1.0 - model_blending_coef) * x.weight + model_blending_coef * (x.mean_error * x.weight);
                        // Method 2: Power-based scaling
                        //effective_weight = x.weight * pow(x.mean_error, model_blending_coef);
                    }

                    prediction_sum += effective_weight * x.prediction;
                    weight_sum += effective_weight;
                }

                if (weight_sum > 0) {
                    if (weight_sum > 1e-10) {
                        (*grid_predictions_map)[i] = prediction_sum / weight_sum;
                    } else {
                        REPORT_ERROR("Weight sum < 1e-10 for vertex %zu in grid predictions\n", i);
                    }
                } else {
                    Rprintf("\n\nERROR: weight_sum: %f for vertex %zu in grid predictions\n", weight_sum, i);
                    Rprintf("grid_wpe_map[%zu].size(): %zu\n", i, wpe_entries.size());

                    size_t counter = 0;
                    for (const auto& x : wpe_entries) {
                        Rprintf("%zu\tx.weight: %f\tx.prediction: %f\n", counter++, x.weight, x.prediction);
                    }

                    REPORT_ERROR("Weight sum = 0 for vertex %zu in grid predictions\n", i);
                }

                // Rprintf("grid_predictions_map.size(): %zu\n", (*grid_predictions_map).size());
            }
        } // END OF if (grid_predictions_map)
    }; // END OF auto average_models


    // Calculates predictions using grid vertex models. Returns void. Accepts memory allocated vector of predictions as a 'return value'
	auto grid_models_predictions = [&](
		const std::vector<double>& y,
		const std::optional<std::vector<double>>& weights,
		std::vector<double>& predictions,
		std::vector<double>& errors,
		std::vector<double>& scale,
		std::optional<std::unordered_map<size_t, double>>& grid_predictions_map
		) {

		// Extract grid vertices into a vector for use inside create_grid_vertex_models
		std::vector<size_t> grid_verts(grid_graph.grid_vertices.begin(), grid_graph.grid_vertices.end());

		// Create a temporary vector to store the model queues for each grid vertex.
		std::vector<ext_ulm_priority_queue> grid_vertex_models_vec(grid_verts.size());

		// Build an index vector for parallelization.
		std::vector<size_t> indices(grid_verts.size());
		std::iota(indices.begin(), indices.end(), 0);

		// Parallel loop: each index corresponds to one grid vertex.
		std::for_each(std::execution::seq, indices.begin(), indices.end(),
					  [&](size_t i) {
						  size_t grid_vertex = grid_verts[i];
						  // Compute models for grid vertex 'gv'
						  grid_vertex_models_vec[i] = create_grid_vertex_models(grid_vertex, y, weights);
					  }
			);

		// Combine the per-vertex results into an unordered_map keyed by grid vertex.
		std::unordered_map<size_t, ext_ulm_priority_queue> grid_vertex_models_map;
		for (size_t i = 0; i < grid_verts.size(); ++i) {
			grid_vertex_models_map[ grid_verts[i] ] = std::move(grid_vertex_models_vec[i]);
		}


		// debugging
		#if 0

		// Use the aggregated models to compute predictions, errors, and scale.
		average_models(grid_vertex_models_map, y, predictions, errors, scale, grid_predictions_map);

		#endif
	};

    std::vector<double> empty_errors;
    std::vector<double> empty_scale;

    if (n_bb > 0) {
        // Perform Bayesian bootstrap estimation

        // Initialize storage for predictions and errors
		result.bb_predictions.resize(n_bb, std::vector<double>(n_vertices, 0.0));

        // Create indices for parallel iteration
        std::vector<int> bb_indices(n_bb);
        std::iota(bb_indices.begin(), bb_indices.end(), 0);

        // Progress tracking
        std::atomic<int> bootstrap_counter{0};

        std::for_each(std::execution::seq, bb_indices.begin(), bb_indices.end(),
                      [&](int iboot) {

						  // Generate bootstrap weights
						  std::mt19937 local_rng((unsigned int)(std::hash<int>{}(iboot))); // Get thread-local RNG
                          std::vector<double> weights(n_vertices);
                          {
                              // Use thread-local RNG to generate weights
                              std::gamma_distribution<double> gamma(1.0, 1.0);
                              double sum = 0.0;
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] = gamma(local_rng);
                                  sum += weights[i];
                              }
                              // Normalize weights
                              for (int i = 0; i < n_vertices; ++i) {
                                  weights[i] /= sum;
                              }
                          }

						  // Generate bootstrapped predictions
                          std::optional<std::vector<double>> opt_weights(weights);
                          std::optional<std::unordered_map<size_t, double>> opt_nullptr_grid_pred;
                          grid_models_predictions(
                              y,
                              opt_weights,
                              result.bb_predictions[iboot],
                              empty_errors,
                              empty_scale,
                              opt_nullptr_grid_pred
                              );

                          // Thread-safe progress update without mutex
                          if (verbose) {
                              int current_count = ++bootstrap_counter;
                              //if (current_count % progress_chunk == 0) {
                              REprintf("\rBootstrap progress: %d%%",
                                       static_cast<int>((100.0 * current_count) / n_bb));
                              //}
                          }
                      });

        bool use_median = true;
        bb_cri_t bb_cri_res = bb_cri(result.bb_predictions, use_median, cri_probability);

        result.predictions = std::move(bb_cri_res.bb_Ey);
        result.cri_lower   = std::move(bb_cri_res.cri_L);
        result.cri_upper   = std::move(bb_cri_res.cri_U);

    } else {
        // Perform standard prediction without bootstrap
        std::optional<std::vector<double>> opt_nullptr;
        std::optional<std::unordered_map<size_t, double>> opt_grid_predictions_map(result.grid_predictions_map);

        grid_models_predictions(
            y,
            opt_nullptr,
            result.predictions,
            result.errors,
            result.scale,
            opt_grid_predictions_map
            );

		// debugging
		#if 0
        if (opt_grid_predictions_map) {
            result.grid_predictions_map = *opt_grid_predictions_map;

			#if 0
            print_umap(result.grid_predictions_map, "result.grid_predictions_map");
            Rprintf("result.grid_predictions_map.size(): %zu\n", result.grid_predictions_map.size());
            Rprintf("n_grid_vertices: %zu\ngrid_graph.grid_vertices.size(): %zu\n",
                    n_grid_vertices, grid_graph.grid_vertices.size());
			#endif
        }
		#endif
    }

    if (n_perms > 0) {
        // Perform permutation testing if requested

        // Initialize permutation test results
		result.null_predictions.resize(n_perms, std::vector<double>(n_vertices, 0.0));

        // Create indices for parallel iteration
        std::vector<int> perm_indices(n_perms);
        std::iota(perm_indices.begin(), perm_indices.end(), 0);

        // Atomic counter for tracking progress
        std::atomic<int> permutation_counter{0};

        // Parallel execution of permutation iterations
		std::for_each(std::execution::seq, perm_indices.begin(), perm_indices.end(),
					  [&](int iperm) {

						  // Create permuted y vector
						  // Create a local RNG for shuffling, seeded with iperm for reproducibility.
						  std::mt19937 local_rng((unsigned int)(std::hash<int>{}(iperm)));
						  std::vector<size_t> local_indices(n_vertices);
						  std::iota(local_indices.begin(), local_indices.end(), 0);
						  std::shuffle(local_indices.begin(), local_indices.end(), local_rng);

						  // Apply permutation using local_indices
						  std::vector<double> y_shuffled(n_vertices);
						  for (size_t i = 0; i < n_vertices; ++i) {
							  y_shuffled[i] = y[local_indices[i]];
						  }

						  // Generate weights using the same local RNG.
						  std::vector<double> weights(n_vertices);
						  {
							  std::gamma_distribution<double> gamma(1.0, 1.0);
							  double sum = 0.0;
							  for (size_t i = 0; i < n_vertices; ++i) {
								  weights[i] = gamma(local_rng);
								  sum += weights[i];
							  }
							  for (size_t i = 0; i < n_vertices; ++i) {
								  weights[i] /= sum;
							  }
						  }

						  // Generate predictions for shuffled y
						  std::optional<std::vector<double>> opt_weights(weights);
						  std::optional<std::unordered_map<size_t, double>> opt_nullptr_grid_pred;
						  grid_models_predictions(
							  y_shuffled,
							  opt_weights,
							  result.null_predictions[iperm],
							  empty_errors,
							  empty_scale,
							  opt_nullptr_grid_pred
							  );

						  if (verbose) {
							  int current_count = ++permutation_counter;
							  REprintf("\rPermutation Test Progress: %d%%",
									   static_cast<int>((100.0 * current_count) / n_perms));
						  }
					  }
			);

        bool use_median = true;
        bb_cri_t null_cri_res = bb_cri(result.null_predictions, use_median, cri_probability);

        result.null_predictions_cri_lower = std::move(null_cri_res.cri_L);
        result.null_predictions_cri_upper = std::move(null_cri_res.cri_U);

        if (n_bb > 0) {
            size_t n_bootstraps = 1000;
            double alpha = 0.05;
            result.permutation_tests = vertex_wasserstein_perm_test(
                result.bb_predictions,
                result.null_predictions,
                n_bootstraps,
                alpha
                );
        }
        if (verbose) {
            Rprintf("\nPermutation Test Completed.\n");
        }
    }

	return result;
}
