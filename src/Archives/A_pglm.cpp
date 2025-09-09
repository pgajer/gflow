#include "A_pglm.h"

/**
 * @brief Performs Path Graph Linear Model fitting with precomputed path models
 *
 * @details This function implements an optimized version of the Path Graph Linear Model
 * algorithm that precomputes and caches linear models for all paths. The algorithm consists of two main phases:
 *
 * Phase 1 - Model Precomputation:
 * - Identifies all unique paths in the graph
 * - For each path, computes and caches linear models for all positions
 * - Uses kernel weights combined with Bayesian bootstrap weights
 *
 * Phase 2 - Vertex Processing:
 * - For each vertex, finds all paths containing it
 * - Selects valid paths based on distance to middle point criteria
 * - Uses precomputed models to find best fit
 * - Optionally performs weighted averaging of predictions
 *
 * @param path_graph Path graph structure containing topology and weights
 * @param y Vector of response variables for each vertex
 * @param weights Vector of Bayesian bootstrap weights for each vertex
 * @param ikernel Kernel function identifier for weight computation
 * @param max_distance_deviation Maximum allowed deviation from optimal position in path
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @param epsilon Small number for numerical stability (default: 1e-8)
 *
 * @return std::pair containing:
 *         - first: Vector of fitted values for each vertex (Ey)
 *         - second: Vector of LOOCV errors for each vertex
 *
 * @throws Runtime error if:
 *         - Input vectors have inconsistent sizes
 *         - No paths found for a vertex
 *         - No valid paths within deviation limits
 *
 * @pre path_graph must be properly initialized with valid topology
 * @pre y.size() == path_graph.vertex_paths.size()
 * @pre weights.size() == path_graph.vertex_paths.size()
 * @pre All weights must be non-negative and sum to 1
 */
std::pair<std::vector<double>, std::vector<double>> pglm_with_precomputed_path_models(
	const path_graph_plm_t& path_graph,
	const std::vector<double>& y,
	const std::vector<double>& weights,
	int ikernel,
	int max_distance_deviation,
	double dist_normalization_factor = 1.01,
	double epsilon = 1e-8) {

#define DEBUG__pglm_with_precomputed_path_models 0

	int h = path_graph.h;
	int desired_path_length = h + 1;
	int n_vertices = path_graph.vertex_paths.size();
	bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

// Initialize cache and collect unique paths
	std::map<itriplet_t, lm_loocv_t> model_cache;

// Initialize working vectors for model computation
	std::vector<double> y_path;
	std::vector<double> x_path;
	std::vector<double> bb_weights_path;
	std::vector<double> w_path;
	std::vector<double> d_path;
	std::vector<std::vector<double>> w_list;

// Phase 1: Precompute models for all paths

#if DEBUG__pglm_with_precomputed_path_models
	int n_shortest_paths = path_graph.shortest_paths.size();
	Rprintf("\n\nn_shortest_paths: %d\n", n_shortest_paths);
	std::vector<int> path_length_freqs(n_vertices, 0); // frequencies of path lengths
#endif

	for (const auto& path_endpoints : path_graph.longest_paths) {
		const std::vector<int>& path = path_graph.shortest_paths.at(path_endpoints);
		int path_n_vertices = path.size();

#if DEBUG__pglm_with_precomputed_path_models
		path_length_freqs[path_n_vertices]++;
#endif

// Resize working vectors for current path
		y_path.resize(path_n_vertices);
		x_path.resize(path_n_vertices);
		bb_weights_path.resize(path_n_vertices);
		w_path.resize(path_n_vertices);
		d_path.resize(path_n_vertices);
		w_list.resize(path_n_vertices);

// Initialize path coordinates
		x_path[0] = 0;

// Extract values and compute path distances
		for (int i = 0; i < path_n_vertices; ++i) {
			y_path[i] = y[path[i]];
			bb_weights_path[i] = weights[path[i]];
			if (i > 0) {
				auto neighbors = path_graph.adj_list[path[i - 1]];
				auto neighbor_weights = path_graph.weight_list[path[i - 1]];
				auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
				int neighbor_i = it - neighbors.begin();
				x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i];
			}
		}

// Compute kernel weights for each position
		for (int pos = 0; pos < path_n_vertices; ++pos) {
// Calculate distances to reference point
			double max_dist = 0.0;
			for (int i = 0; i < path_n_vertices; ++i) {
				d_path[i] = std::abs(x_path[i] - x_path[pos]);
				max_dist = std::max(max_dist, d_path[i]);
			}
			if (max_dist == 0) max_dist = 1;
			max_dist *= dist_normalization_factor;

// Normalize distances and compute kernel weights
			for (int i = 0; i < path_n_vertices; ++i) {
				d_path[i] /= max_dist;
			}

#pragma omp critical(kernel_init)
			{
				initialize_kernel(ikernel, 1.0);
			}
			kernel_fn(d_path.data(), path_n_vertices, w_path.data());

// Normalize kernel weights
			double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
			for (int i = 0; i < path_n_vertices; ++i) {
				w_path[i] = (w_path[i] / total_w_path) * bb_weights_path[i];
			}

			w_list[pos].assign(w_path.begin(), w_path.end());
		}

// Compute models for all positions
		auto models = predict_lms_1d_loocv(y_path, x_path, path, w_list, y_binary, epsilon);

// Cache computed models
		for (int pos = 0; pos < path_n_vertices; ++pos) {
			model_cache.emplace(itriplet_t{path_endpoints.first, path_endpoints.second, pos},
								models[pos]);
		}
	}

#if DEBUG__pglm_with_precomputed_path_models
	Rprintf("path_length_freqs\n");
	Rprintf("len\tfreqs\n");
	for (int i = 0; i < n_vertices; ++i) {
		if (path_length_freqs[i] > 0) {
			Rprintf("%d:\t%d\n", i, path_length_freqs[i]);
		}
	}
#endif

// Phase 2: Process each vertex using precomputed models
	std::vector<double> Ey(n_vertices);
	std::vector<double> errors(n_vertices);
	std::vector<lm_loocv_t> best_models(n_vertices);

	for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
		auto vertex_path_graph_info = path_graph.vertex_paths[vertex_i];

// Get all paths containing the vertex
		std::vector<std::pair<std::vector<int>, int>> vertex_paths =
			vertex_path_graph_info.get_full_paths(desired_path_length, path_graph.shortest_paths);

		if (vertex_paths.empty()) {
			Rf_error("No paths found for vertex %d", vertex_i);
		}

// Debugging: check if there are any paths of length != h + 1
		for (const auto& path : vertex_paths) {
			int path_len = path.first.size();
			if (path_len != h + 1) {
				Rprintf("vertex_i: %d\n", vertex_i);
				REPORT_ERROR("Found path of length %d that is not equal to desired_path_length: %d\n", path_len, desired_path_length);
			}
		}

// Find valid paths based on distance to middle point
		std::vector<size_t> valid_path_indices;
		valid_path_indices.reserve(vertex_paths.size());

		int mid_pt = (desired_path_length - 1) / 2;
		int min_dist_to_mid_pt = h;

// Find minimum distance to middle point
		for (const auto& path : vertex_paths) {
			min_dist_to_mid_pt = std::min(min_dist_to_mid_pt, std::abs(path.second - mid_pt));
		}

// Collect paths within allowed deviation
		for (size_t path_i = 0; path_i < vertex_paths.size(); ++path_i) {
			if (std::abs(vertex_paths[path_i].second - mid_pt) <=
				min_dist_to_mid_pt + max_distance_deviation) {
				valid_path_indices.push_back(path_i);
			}
		}

		if (valid_path_indices.empty()) {
			Rf_error("No valid paths found for vertex %d within deviation limits", vertex_i);
		}

// Find best model among valid paths
		double best_loocv_error = std::numeric_limits<double>::infinity();
		lm_loocv_t best_model;

		for (const auto& path_i : valid_path_indices) {
			const auto& path = vertex_paths[path_i].first;
			int position_in_path = vertex_paths[path_i].second;

			auto model = model_cache.at(itriplet_t{path[0], path.back(), position_in_path});

			if (model.loocv_at_ref_vertex < best_loocv_error) {
				best_loocv_error = model.loocv_at_ref_vertex;
				best_model = model;
			}
		}

		best_models[vertex_i] = best_model;
		Ey[vertex_i] = best_model.predicted_value;
		errors[vertex_i] = best_model.loocv_at_ref_vertex;
	}

// Phase 3: Perform model averaging
	double i_best_model_predicted_value = 0;
	for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
		const auto& best_model = best_models[vertex_i];
		const auto& path_vertex_indices = best_model.vertex_indices;
		const auto& path_vertex_weights = best_model.w_values;
#if 0
		Rprintf("vertex_i: %d\n", vertex_i);
		print_vect(path_vertex_indices, "path_vertex_indices");
		print_vect(path_vertex_weights, "path_vertex_weights");
#endif

// Compute weighted average of predictions
		double weight = path_vertex_weights[best_model.ref_index];
		double total_weight = weight;
		double vertex_Ey = weight * Ey[vertex_i];

		for (size_t i = 0; i < path_vertex_indices.size(); ++i) {
			int index = path_vertex_indices[i];
			if (index != vertex_i) {
				const auto& i_best_model = best_models[index];

// Check if vertex_i is in the model's training set
				auto it = std::find(i_best_model.vertex_indices.begin(),
									i_best_model.vertex_indices.end(),
									vertex_i);

				if (it != i_best_model.vertex_indices.end()) {
					weight = path_vertex_weights[i];
					total_weight += weight;
					i_best_model_predicted_value = i_best_model.predict(vertex_i);
					if (y_binary) i_best_model_predicted_value = std::clamp(i_best_model_predicted_value, 0.0, 1.0);
					vertex_Ey += weight * i_best_model_predicted_value;
				}
			}
		}

		Ey[vertex_i] = vertex_Ey / total_weight;
	}

return std::make_pair(Ey, errors);
}


#if 0
// New helper struct to store model information for each vertex
struct vertex_model_info_t {
    std::vector<lm_loocv_t> models;           // All models containing this vertex
    std::vector<double> model_weights;        // Original kernel weights for each model
    std::vector<int> positions_in_models;     // Position of vertex in each model
};

std::pair<std::vector<double>, std::vector<double>> pglm_with_precomputed_path_models(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon) {

    int h = path_graph.h;
    int path_n_vertices = h + 1;
    int n_vertices = path_graph.vertex_paths.size();
    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    // Initialize storage for all vertex models
    std::vector<vertex_model_info_t> vertex_models(n_vertices);

    // Initialize cache and working vectors
    std::map<itriplet_t, lm_loocv_t> model_cache;
    std::vector<double> y_path, x_path, bb_weights_path, w_path, d_path;
    std::vector<std::vector<double>> w_list;

    // Phase 1: Precompute models and store relevant information for each vertex
    for (const auto& path_endpoints : path_graph.longest_paths) {
        const std::vector<int>& path = path_graph.shortest_paths.at(path_endpoints);
        int path_n_vertices = path.size();

        // Resize working vectors
        y_path.resize(path_n_vertices);
        x_path.resize(path_n_vertices);
        bb_weights_path.resize(path_n_vertices);
        w_path.resize(path_n_vertices);
        d_path.resize(path_n_vertices);
        w_list.resize(path_n_vertices);

        // Initialize path coordinates and extract values
        x_path[0] = 0;
        for (int i = 0; i < path_n_vertices; ++i) {
            y_path[i] = y[path[i]];
            bb_weights_path[i] = weights[path[i]];
            if (i > 0) {
                auto neighbors = path_graph.adj_list[path[i - 1]];
                auto neighbor_weights = path_graph.weight_list[path[i - 1]];
                auto it = std::find(neighbors.begin(), neighbors.end(), path[i]);
                int neighbor_i = it - neighbors.begin();
                x_path[i] = x_path[i - 1] + neighbor_weights[neighbor_i];
            }
        }

        // Compute kernel weights for each position
        for (int pos = 0; pos < path_n_vertices; ++pos) {
            // Calculate distances to reference point
            double max_dist = 0.0;
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] = std::abs(x_path[i] - x_path[pos]);
                max_dist = std::max(max_dist, d_path[i]);
            }
            if (max_dist == 0) max_dist = 1;
            max_dist *= dist_normalization_factor;

            // Normalize distances
            for (int i = 0; i < path_n_vertices; ++i) {
                d_path[i] /= max_dist;
            }

            #pragma omp critical(kernel_init)
            {
                initialize_kernel(ikernel, 1.0);
            }
            kernel_fn(d_path.data(), path_n_vertices, w_path.data());

            // Normalize kernel weights
            double total_w_path = std::accumulate(w_path.begin(), w_path.end(), 0.0);
            for (int i = 0; i < path_n_vertices; ++i) {
                w_path[i] = (w_path[i] / total_w_path) * bb_weights_path[i];
            }

            w_list[pos].assign(w_path.begin(), w_path.end());
        }

        // Compute models for all positions
        auto models = predict_lms_1d_loocv(y_path, x_path, path, w_list, y_binary, epsilon);

        // Store models and associate them with vertices
        int mid_vertex = (path_n_vertices - 1) / 2;

        for (int pos = 0; pos < path_n_vertices; ++pos) {
            auto& model = models[pos];
            model_cache.emplace(itriplet_t{path_endpoints.first, path_endpoints.second, pos}, model);

            // Check if this position is within allowed deviation from middle point
            int dist_to_mid = std::abs(pos - mid_vertex);
            if (dist_to_mid <= max_distance_deviation) {
                // Store model information for each vertex in the path
                for (size_t i = 0; i < path.size(); ++i) {
                    int vertex = path[i];
                    vertex_models[vertex].models.push_back(model);
                    vertex_models[vertex].model_weights.push_back(w_list[pos][i]);
                    vertex_models[vertex].positions_in_models.push_back(i);
                }
            }
        }
    }

    // Phase 2: Compute model-averaged predictions for each vertex
    std::vector<double> Ey(n_vertices);
    std::vector<double> errors(n_vertices);

    for (int vertex_i = 0; vertex_i < n_vertices; ++vertex_i) {
        const auto& vertex_info = vertex_models[vertex_i];

        if (vertex_info.models.empty()) {
            Rf_error("No valid models found for vertex %d", vertex_i);
        }

        // Compute weighted average prediction and error
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        double weighted_error = 0.0;

        for (size_t i = 0; i < vertex_info.models.size(); ++i) {
            const auto& model = vertex_info.models[i];
            double weight = vertex_info.model_weights[i];
            int pos = vertex_info.positions_in_models[i];

            double prediction = model.predictions[pos];
            if (y_binary) {
                prediction = std::clamp(prediction, 0.0, 1.0);
            }

            weighted_sum += weight * prediction;
            weighted_error += weight * model.errors[pos];
            weight_sum += weight;
        }

        Ey[vertex_i] = weighted_sum / weight_sum;
        errors[vertex_i] = weighted_error / weight_sum;
    }

    return std::make_pair(Ey, errors);
}
#endif


#if 0
std::pair<std::vector<double>, std::vector<double>> pglm_with_precomputed_path_models(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int ikernel,
    int max_distance_deviation,
    double dist_normalization_factor = 1.01,
    double epsilon = 1e-8);
#endif
