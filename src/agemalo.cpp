#include <thread>      // For std::thread
#include <random>      // For std::mt19937, std::random_device, std::gamma_distribution

// #include <filesystem>  // for debugging
#include <fstream>

#include "agemalo.hpp"
#include "kernels.h"
#include "error_utils.h"         // For REPORT_ERROR()
#include "predictive_errors.hpp" // For struct bb_cri_t
#include "opt_bw.hpp"
#include "cpp_utils.hpp"         // For debugging and elapsed.time
#include "progress_utils.hpp"    // For elapsed.time
#include "bandwidth_utils.hpp"   // For get_candidate_bws()

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
 * @param log_grid If true, use logarithmic spacing; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param snap_tolerance Maximum distance for snapping original vertices to grid vertices
 * @param dist_normalization_factor Factor for normalizing distances in the graph
 * @param kernel_type Type of kernel function used for weight calculation (e.g., Gaussian, triangular)
 * @param tolerance Convergence tolerance for linear model fitting
 * @param n_bb Number of Bayesian bootstrap iterations (0 for no bootstrap)
 * @param cri_probability Probability level for confidence intervals in bootstrap
 * @param n_perms Number of permutation test iterations (0 for no permutation testing)
 * @param blending_coef Number between 0 and 1 allowing for smooth interpolation between position-based weights and mean-Rf_error times position-based weights using the formula effective_weight = x.weight * pow(x.mean_error, blending_coef);
 * @param verbose Whether to print progress information
 *
 * @return adaptive_uggmalo_result_t structure containing:
 * - graph_diameter: Computed diameter of the input graph
 * - grid_opt_bw: Optimal bandwidth for each grid vertex
 * - predictions: Model-averaged predictions for original graph vertices
 * - grid_predictions: Model-averaged predictions for grid vertices
 * - Bootstrap results (if n_bb > 0):
 *   - bb_predictions: Bootstrap replication predictions
 *   - cri_lower: Lower confidence interval bounds
 *   - cri_upper: Upper confidence interval bounds
 * - Permutation test results (if n_perms > 0):
 *   - null_predictions: Predictions under null hypothesis
 *   - permutation_tests: Statistical test results including p-values and effect sizes
 *
 */
agemalo_result_t agemalo(
	const uniform_grid_graph_t& grid_graph,
	const std::vector<double>& y,
	// geodesic parameter
	size_t min_path_size,
	// bw parameters
	size_t n_bws,
	bool log_grid,
	double min_bw_factor,
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

	// Start timer for the entire function
	auto total_start_time = std::chrono::steady_clock::now();

	// Initialize result structure
	agemalo_result_t result;
	result.graph_diameter = grid_graph.graph_diameter;
	result.max_packing_radius = grid_graph.max_packing_radius;

	// Initialize vectors to store predictions for all vertices
	const size_t n_vertices = grid_graph.adjacency_list.size();
	result.predictions.resize(n_vertices, INFINITY);
	result.errors.resize(n_vertices, INFINITY);
	result.scale.resize(n_vertices, INFINITY);

	auto is_binary01 = [](const std::vector<double>& yy, double tol = 1e-12) -> bool {
		for (double v : yy) {
            if (!(std::fabs(v) <= tol || std::fabs(v - 1.0) <= tol)) {
                return false;
            }
        }
        return true;
    };

    const bool y_binary = is_binary01(y);

	// Define minimum and maximum bandwidth as a fraction of graph diameter
	double max_bw = max_bw_factor * grid_graph.graph_diameter;
	double min_bw = min_bw_factor * grid_graph.graph_diameter;
	// Ensure the min_bw is at least the maximum packing radius
	// This guarantees all vertices are covered by at least one ball
	min_bw = std::max(min_bw, grid_graph.max_packing_radius);

	initialize_kernel(kernel_type, 1.0);

	if (verbose) {
		Rprintf("Starting AGEMALO algorithm\n");
		Rprintf("n_vertices(grid_graph): %zu\n", grid_graph.adjacency_list.size());
		// packing parameters
		Rprintf("graph_diameter: %f\n", result.graph_diameter);
		Rprintf("max_packing_radius: %f\n", result.max_packing_radius);
		Rprintf("grid_graph.grid_vertices.size(): %zu\n", grid_graph.grid_vertices.size());

		// bw parameters
		Rprintf("n_bws: %zu\n", n_bws);
		Rprintf("min_bw_factor: %f\n", min_bw_factor);
		Rprintf("min_bw: %f\n", min_bw);
		Rprintf("max_bw_factor: %f\n", max_bw_factor);
		Rprintf("max_bw: %f\n", max_bw);
	}

#define DEBUG__agemalo 0
#if DEBUG__agemalo

	std::string debug_dir = "/Users/pgajer/current_projects/msr2/debugging_data/";
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
#endif

	// weight/prediction/Rf_error/mean_error/be struct needed for mode averaging and local scale estimation; we use it to record weight/prediction/Rf_error/bw of the given vertex in each model where the vertext is in the support of the model
	struct wpe_t {
		double weight;
		double prediction;
		double Rf_error;
		double mean_error;
		double bw;

		// Constructor needed for emplace_back(x,y,z)
		wpe_t(double w, double p, double e, double me, double bw)
			: weight(w), prediction(p), Rf_error(e), mean_error(me), bw(bw) {}
	};

	//
	// Lambda for creating models at grid vertices - designed to be thread safe
	//
	auto create_grid_vertex_models = [&](
		size_t grid_vertex,
		const std::vector<double>& y,
		const std::optional<std::vector<double>>& weights
		) {

		auto vertex_start_time = std::chrono::steady_clock::now();
		bool timing_needed = false;

		// Find the minimum bandwidth that ensures at least one geodesic ray has the minimum required size
		double precision = 0.001;
		double min_min_bw = grid_graph.find_minimum_bandwidth(
			grid_vertex,
			min_bw,
			max_bw,
			min_path_size,
			precision
			);

		if (min_bw < min_min_bw) {
			min_bw = min_min_bw;
		}

		if (min_bw >= max_bw) {
			// max_bw = min_bw * 1.1; // Increase max slightly
			// Rprintf("Warning: Minimum bandwidth (%f) >= maximum bandwidth (%f) for vertex %zu. Adjusted max to %f\n",
			// 		min_bw, max_bw_factor * grid_graph.graph_diameter, grid_vertex, max_bw);
			// timing_needed = true;
			REPORT_ERROR("Minimum bandwidth (%f) >= maximum bandwidth (%f)\n",
						 min_bw, max_bw_factor * grid_graph.graph_diameter);
		}

		if (verbose && grid_vertex % 10 == 0) {
			Rprintf("Processing grid vertex %zu, min_bw: %.4f, max_bw: %.4f\n",
					grid_vertex, min_bw, max_bw);
			timing_needed = true;
		}

#if 0
		std::vector<double> candidate_bws = get_candidate_bws(
			min_bw,
			max_bw,
			n_bws,
			log_grid,
			min_spacing
			);


		// ToDo: change the signature of find_optimal_bandwidth_over_grid() to

		opt_bw_t uniform_grid_graph_t::find_optimal_bandwidth_over_grid(
			size_t grid_vertex,
			const std::vector<double>& y,
			const std::vector<double>& bws,
			double dist_normalization_factor,
			bool y_binary,
			double tolerance,
			double precision,
			bool verbose,
			const std::optional<std::vector<double>>& weights) const;

#endif


		// Find optimal bandwidth and corresponding models using binary search
		auto opt_result = grid_graph.find_optimal_bandwidth_over_grid(
			grid_vertex,
			y,
			min_path_size,
			min_bw,
			max_bw,
			n_bws,
			log_grid,
			dist_normalization_factor,
			y_binary,
			model_tolerance,
			precision,
			verbose,
			weights
			);

		double opt_bw = opt_result.opt_bw;
		std::vector<ext_ulm_t> opt_models = std::move(opt_result.models);
		// we can also extract bws and errors for debugging

		if (opt_models.size() == 0) {
			REPORT_ERROR("No models found for grid vertex %zu\n", grid_vertex);
		}


#if 0

		// Track models originating from this grid vertex for debugging
		static std::unordered_map<size_t, std::vector<std::pair<size_t, ext_ulm_t>>> grid_vertex_to_models;

		// Write out bws and errors for each grid_vertex model
		std::ofstream model_file(debug_dir + "vertex_" + std::to_string(debug_vertex) +
								 "_model_" + std::to_string(model_idx) + "_from_grid_" +
								 std::to_string(grid_vertex) + ".csv");
		if (model_file.is_open()) {
			// ... write model details ...
		}


		// for debugging
		static std::unordered_map<size_t, std::vector<size_t>> grid_vertex_to_debug_models;
		const size_t debug_vertex = 4; // The vertex we're debugging

		for (size_t model_idx = 0; model_idx < opt_models.size(); ++model_idx) {
			const auto& model = opt_models[model_idx];

			// Check if this model affects our debug vertex
			bool affects_debug = false;
			for (size_t vertex : model.vertices) {
				if (vertex == debug_vertex) {
					affects_debug = true;
					break;
				}
			}

			// If this model affects our debug vertex, add it to our mapping
			if (affects_debug) {
				grid_vertex_to_debug_models[grid_vertex].push_back(model_idx);

				// Also write out the model details
				std::ofstream model_file(debug_dir + "vertex_" + std::to_string(debug_vertex) +
										 "_model_" + std::to_string(model_idx) + "_from_grid_" +
										 std::to_string(grid_vertex) + ".csv");
				if (model_file.is_open()) {
					// ... write model details ...
				}
			}
		}
#endif
		
		// Store the optimal bandwidth for this grid vertex
		result.grid_opt_bw.push_back(opt_bw);

		std::unordered_map<size_t, std::vector<wpe_t>> wpe_map; // wpe[i] stores a vector of {weight, prediction, Rf_error, mean_error} values for each model that contains the i-th vertex in its support; these values will be used to compute the model averaged predictions
		for (const auto& model : opt_models) {
			for (size_t i = 0; i < model.vertices.size(); ++i) {
				wpe_map[ model.vertices[i] ].emplace_back(
					model.w_path[i],
					model.predictions[i],
					model.errors[i],
					model.mean_error,
					model.bw
					);
			}
		}

		if (verbose && (grid_vertex % 10 == 0 || timing_needed)) {
			Rprintf("Found opt_bw: %.4f with n(opt_models): %zu for vertex %zu\n",
					opt_bw, opt_models.size(), grid_vertex);
			elapsed_time(vertex_start_time, "Grid vertex processing completed", true);
		}

#if 0
		// debugging
		// At the end of all grid vertex processing, write a summary of which grid vertices
		// contributed models to our debug vertex
		if (verbose && grid_vertex == grid_verts.back()) {
			std::ofstream summary_file(debug_dir + "vertex_" + std::to_string(debug_vertex) + "_grid_contributors.csv");
			if (summary_file.is_open()) {
				summary_file << "grid_vertex,num_models,model_ids\n";
				for (const auto& [gv, models] : grid_vertex_to_debug_models) {
					summary_file << gv << "," << models.size() << ",\"";
					for (size_t i = 0; i < models.size(); ++i) {
						if (i > 0) summary_file << ",";
						summary_file << models[i];
					}
					summary_file << "\"\n";
				}
				summary_file.close();
			}
		}
#endif

		return wpe_map;
	};

	auto average_models = [&model_blending_coef,&n_vertices](
		std::vector<std::vector<wpe_t>>& wpe,
		std::vector<double>& predictions,
		std::vector<double>& errors,
		std::vector<double>& scale
		) {

		bool process_errors = errors.size() == predictions.size();
		bool process_scale  = scale.size() == predictions.size();

		// Model averaging over the vertices of the original graph
		for (size_t i = 0; i < n_vertices; i++ ) {

			double prediction_sum = 0.0;
			double weight_sum = 0.0;
			double error_sum  = 0.0;
			double local_scale= 0.0; //INFINITY;
			double effective_weight;

			for (const auto& x : wpe[i]) {

				if (wpe[i].empty()) {

					REPORT_ERROR("Vertex %zu was not in the domain of any model\n",i);

					// Handle vertices with no model predictions
					// predictions[i] = y[i]; // Use original value or other fallback
					// 	if (process_errors) errors[i] = INFINITY; // Mark as high Rf_error
					// 	if (process_scale) scale[i] = INFINITY;  // Mark as no reliable scale
					// 	continue; // Skip to next vertex
				}

				if (model_blending_coef == 0.0) {
					// Pure position weight
					effective_weight = x.weight;
				}
				else if (model_blending_coef == 1.0) {
					// Full mean Rf_error influence
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

				if (process_errors) error_sum += effective_weight * x.Rf_error;
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
					REPORT_ERROR("Very small weight sum encountered for vertex %zu.\n", i); // Setting predictions[i] to y[i]
					//predictions[i] = std::numeric_limits<double>::quiet_NaN();
				}
			} else {
				REPORT_ERROR("Weight sum = 0 for vertex %zu in predictions\n", i);
				//predictions[i] = y[i];
			}
		}
	}; // END OF auto average_models


	auto grid_models_start_time = std::chrono::steady_clock::now();

	// Calculates predictions using grid vertex models. Returns void. Accepts memory allocated vector of predictions as a 'return value'
	auto grid_models_predictions = [&](
		const std::vector<double>& y,
		const std::optional<std::vector<double>>& weights,
		std::vector<double>& predictions,
		std::vector<double>& errors,
		std::vector<double>& scale
		) {

		// Extract grid vertices into a vector for use
		std::vector<size_t> grid_verts(grid_graph.grid_vertices.begin(), grid_graph.grid_vertices.end());
		if (verbose) {
			Rprintf("Processing %zu grid vertices\n", grid_verts.size());
		}

		// wpe[i] stores a vector of {weight, prediction, Rf_error, mean_error}
		// values for each model that contains the i-th vertex in its support;
		// these values will be used to compute the model-averaged predictions
		std::vector<std::vector<wpe_t>> wpe(n_vertices);

		// Progress counter
		std::atomic<size_t> progress_counter{0};
		const size_t progress_step = std::max(size_t(1), grid_verts.size() / 10);

		// Serial loop: each index corresponds to one grid vertex
		for (size_t i = 0; i < grid_verts.size(); ++i) {
			size_t grid_vertex = grid_verts[i];
			auto wpe_map = create_grid_vertex_models(grid_vertex, y, weights);

			for (const auto& [vertex, vertex_wpe_vector] : wpe_map) {
				wpe[vertex].insert(wpe[vertex].end(),
								   vertex_wpe_vector.begin(),
								   vertex_wpe_vector.end());
			}

			// Report progress
			if (verbose) {
				size_t completed = ++progress_counter;
				if (completed % progress_step == 0 || completed == grid_verts.size()) {
					double percentage = (100.0 * completed) / grid_verts.size();
					Rprintf("\rProcessing grid vertices: %.1f%% (%zu/%zu)",
							percentage, completed, grid_verts.size());
					R_FlushConsole();
				}
			}
		}

		if (verbose) {
			Rprintf("\n");
			elapsed_time(grid_models_start_time, "Grid vertex models creation completed", true);

			// Collect statistics on number of models per vertex
			std::vector<size_t> model_counts(n_vertices);
			for (size_t i = 0; i < n_vertices; ++i) {
				model_counts[i] = wpe[i].size();
			}

			// Print summary statistics
			print_vector_quantiles(model_counts, "Models per vertex", {0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0});

			// Print Rf_warning for vertices with no models
			size_t zero_model_count = std::count(model_counts.begin(), model_counts.end(), 0);
			if (zero_model_count > 0) {
				REPORT_WARNING("\nWarning: %zu vertices (%.1f%%) have no models\n",
							   zero_model_count, 100.0 * zero_model_count / n_vertices);
			}
		}

		// Use the aggregated models to compute predictions, errors, and scale.
		auto averaging_start_time = std::chrono::steady_clock::now();
		average_models(wpe, predictions, errors, scale);

		#if 0
		// debugging
		{
			// Debug specific vertex (vertex 5 in this case)
			size_t debug_vertex = 4;

			// Open debug file
			std::ofstream debug_file(debug_dir + "vertex_" + std::to_string(debug_vertex) + "_debug.csv");
			if (debug_file.is_open()) {
				debug_file << "model_id,grid_vertex,path_length,n_vertices,weight,prediction,Rf_error,mean_error,bandwidth\n";

				// For each model that influences this vertex
				size_t model_id = 0;
				for (const auto& x : wpe[debug_vertex]) {
					// Calculate effective weight based on model_blending_coef
					double effective_weight;
					if (model_blending_coef == 0.0) {
						effective_weight = x.weight;
					} else if (model_blending_coef == 1.0) {
						effective_weight = x.weight * x.mean_error;
					} else {
						effective_weight = (1.0 - model_blending_coef) * x.weight +
							model_blending_coef * (x.weight * x.mean_error);
					}

					// Write model details
					debug_file << model_id << ","
							   << "," // We don't directly know which grid vertex this came from
							   << "," // Path length not directly available
							   << "," // Number of vertices not directly available
							   << effective_weight << ","
							   << x.prediction << ","
							   << x.Rf_error << ","
							   << x.mean_error << ","
							   << x.bw << "\n";

					model_id++;
				}
				debug_file.close();

				// Write summary file with overall prediction information
				std::ofstream summary_file(debug_dir + "vertex_" + std::to_string(debug_vertex) + "_summary.txt");
				if (summary_file.is_open()) {
					summary_file << "Vertex: " << debug_vertex << "\n";
					summary_file << "True y value: " << y[debug_vertex] << "\n";
					summary_file << "AGEMALO prediction: " << result.predictions[debug_vertex] << "\n";
					summary_file << "Number of contributing models: " << wpe[debug_vertex].size() << "\n";
					//summary_file << "Error: " << result.errors[debug_vertex] << "\n";
					//summary_file << "Scale: " << result.scale[debug_vertex] << "\n";
					summary_file.close();
				}

				if (verbose) {
					Rprintf("Debug information for vertex %zu written to %s\n",
							debug_vertex, (debug_dir + "vertex_" + std::to_string(debug_vertex) + "_*.csv").c_str());
				}
			}
		}
		#endif

		if (verbose) {
			elapsed_time(averaging_start_time, "Model averaging completed", true);
		}
	};

	std::vector<double> empty_errors;
	std::vector<double> empty_scale;

	if (n_bb > 0) {
		// Perform Bayesian bootstrap estimation
		auto bootstrap_start_time = std::chrono::steady_clock::now();
		if (verbose) {
			Rprintf("Starting Bayesian bootstrap with %zu iterations\n", n_bb);
		}

		// Initialize storage for predictions and errors
		result.bb_predictions.resize(n_bb, std::vector<double>(n_vertices, 0.0));

		// Progress tracking
		std::atomic<int> bootstrap_counter{0};

		for (size_t iboot = 0; iboot < n_bb; ++iboot) {
			// Generate bootstrap weights
			std::mt19937 local_rng(static_cast<unsigned int>(std::hash<int>{}(iboot)));
			std::vector<double> weights(n_vertices);
			{
				std::gamma_distribution<double> gamma(1.0, 1.0);
				double sum = 0.0;
				for (size_t i = 0; i < n_vertices; ++i) {
					weights[i] = gamma(local_rng);
					sum += weights[i];
				}

				if (sum > 0.0) {
					// Normalize weights
					for (size_t i = 0; i < n_vertices; ++i) {
						weights[i] /= sum;
					}
				} else {
					// Handle edge case of all-zero weights
					std::fill(weights.begin(), weights.end(), 1.0 / static_cast<double>(n_vertices));
					REPORT_WARNING("Bootstrap instance %d generated all-zero weights. Using uniform weights.\n", iboot);
				}
			}

			// Generate bootstrapped predictions
			std::optional<std::vector<double>> opt_weights(weights);
			grid_models_predictions(
				y,
				opt_weights,
				result.bb_predictions[iboot],
				empty_errors,
				empty_scale
				);

			// Progress update
			if (verbose) {
				int current_count = ++bootstrap_counter;
				REprintf("\rBootstrap progress: %d%%",
						 static_cast<int>((100.0 * current_count) / n_bb));
			}
		}

		bool use_median = true;
		bb_cri_t bb_cri_res = bb_cri(result.bb_predictions, use_median, cri_probability);

		result.predictions = std::move(bb_cri_res.bb_Ey);
		result.cri_lower   = std::move(bb_cri_res.cri_L);
		result.cri_upper   = std::move(bb_cri_res.cri_U);

		if (verbose) {
			Rprintf("\n");
			elapsed_time(bootstrap_start_time, "Bayesian bootstrap completed", true);
		}

	} else {
// Perform standard prediction without bootstrap
		std::optional<std::vector<double>> opt_nullptr;

		grid_models_predictions(
			y,
			opt_nullptr,
			result.predictions,
			result.errors,
			result.scale
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
		auto permutation_start_time = std::chrono::steady_clock::now();

		if (verbose) {
			Rprintf("Starting permutated outcome predictions estimates with %zu permutations\n", n_perms);
		}

		// Initialize permutation test results
		result.null_predictions.resize(n_perms, std::vector<double>(n_vertices, 0.0));

		// Progress tracking
		std::atomic<int> permutation_counter{0};

		for (size_t iperm = 0; iperm < n_perms; ++iperm) {
			// Create permuted y vector
			// Local RNG for shuffling, seeded with iperm for reproducibility
			std::mt19937 local_rng(static_cast<unsigned int>(std::hash<int>{}(iperm)));
			std::vector<size_t> local_indices(n_vertices);
			std::iota(local_indices.begin(), local_indices.end(), 0);
			std::shuffle(local_indices.begin(), local_indices.end(), local_rng);

			// Apply permutation using local_indices
			std::vector<double> y_shuffled(n_vertices);
			for (size_t i = 0; i < n_vertices; ++i) {
				y_shuffled[i] = y[local_indices[i]];
			}

			// Generate weights using the same local RNG
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
			grid_models_predictions(
				y_shuffled,
				opt_weights,
				result.null_predictions[iperm],
				empty_errors,
				empty_scale
				);

			if (verbose) {
				int current_count = ++permutation_counter;
				REprintf("\rPermutation Test Progress: %d%%",
						 static_cast<int>((100.0 * current_count) / n_perms));
			}
		}

		bool use_median = true;
		bb_cri_t null_cri_res = bb_cri(result.null_predictions, use_median, cri_probability);

		result.null_predictions_cri_lower = std::move(null_cri_res.cri_L);
		result.null_predictions_cri_upper = std::move(null_cri_res.cri_U);

		if (verbose) {
			Rprintf("\n");
			elapsed_time(permutation_start_time, "Permutated outcome predictions estimates completed", true);
		}

		if (n_bb > 0) {
			auto permutation_tests_start_time = std::chrono::steady_clock::now();

			size_t n_bootstraps = 1000;
			double alpha = 0.05;
			result.permutation_tests = vertex_wasserstein_perm_test(
				result.bb_predictions,
				result.null_predictions,
				n_bootstraps,
				alpha
				);

			if (verbose) {
				Rprintf("\n");
				elapsed_time(permutation_tests_start_time, "Permutation testing completed", true);
			}
		}
	}

	if (verbose) {
		elapsed_time(total_start_time, "AGEMALO algorithm completed", true);
	}

	return result;
}
