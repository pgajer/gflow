#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>        // For std::vector used throughout the code
#include <unordered_map> // For std::unordered_map used for boundary_vertices
#include <map>
#include <unordered_set>
#include <set>
#include <optional>
#include <algorithm>
#include <span>
#include <ranges>
#include <limits>
#include <utility>            // For std::pair

#include <Eigen/Core>
#include <Eigen/Dense>  // For Eigen::MatrixXd
#include <Eigen/Sparse> // For Eigen::SparseMatrix, Triplet

#include "iknn_graphs.hpp"
#include "edge_info.hpp"
#include "edge_weights.hpp"
#include "vertex_info.hpp"
#include "vertex_path.hpp"
#include "reachability_map.hpp"
#include "explored_tracker.hpp"
#include "circular_param_result.hpp"
#include "weighted_correlation.hpp"
#include "gradient_flow.hpp"
#include "ulm.hpp"
#include "opt_bw.hpp"
#include "graph_spectral_lowess.hpp"     // For graph_spectral_lowess_t
#include "graph_spectral_lowess_mat.hpp" // For graph_spectral_lowess_mat_t
#include "edge_pruning_stats.hpp"        // For edge_pruning_stats_t
#include "nada_graph_spectral_lowess.hpp"// For nada_graph_spectral_lowess_t
#include "graph_deg0_lowess_cv.hpp"
#include "graph_deg0_lowess_cv_mat.hpp"
#include "graph_deg0_lowess_buffer_cv.hpp"
#include "graph_kernel_smoother.hpp"
#include "graph_bw_adaptive_spectral_smoother.hpp"
#include "klaps_low_pass_smoother.hpp"
#include "klaps_spectral_filter.hpp"
#include "basin.hpp"
#include "invalid_vertex.hpp"
#include "harmonic_smoother.hpp"

// Undefine conflicting macros after including R headers
#undef length


struct path_t {
	std::vector<size_t> vertices;  ///< vertices of the path
	std::vector<double> distances; ///< distance of each vertex to the reference vertex
	size_t ref_vertex_index = 0;   ///< index of the reference vertex within vertices (and hence distances)
	double total_weight = 0.0;

	// For priority queue ordering; sort paths in the descending order of their total_length's
	bool operator<(const path_t& other) const {
		return total_weight > other.total_weight;
	}
};

// shortest_paths_t allows reconstruction of shortest path between the given ref
// vertex and any vertex reachable from the ref vertex
struct subpath_t {
   size_t path_idx;   ///< path index
   size_t vertex_idx; ///< index of the given vertex of reachable_vertices set within the path with index path_idx; for example if path = {4,20,31,22} and v = 31, then vertex_idx for v is 2 as path[2] = 31
};
struct shortest_paths_t {
	std::vector<path_t> paths;
	std::unordered_set<int> reachable_vertices;
	std::unordered_map<size_t, subpath_t> vertex_to_path_map; // Maps vertex to path index
};

struct gray_xyw_t {
	std::vector<size_t> vertices;  ///< Indices of the original graph vertices forming the path
	std::vector<double> x_path;///< Cumulative distance along path from initial vertex
	std::vector<double> y_path;///< Y-values of vertices restricted to the path
	std::vector<double> w_path;///< Kernel-weighted distances from reference vertex

// Method to evaluate linearity
	double evaluate_linearity() const {
		double corr = calculate_weighted_correlation(x_path, y_path, w_path);
		return corr * corr;
	}
};

struct composite_shortest_paths_t : public shortest_paths_t {
	std::vector<std::pair<size_t, size_t>> composite_paths;

	explicit composite_shortest_paths_t(const shortest_paths_t& shortest_paths)
		: shortest_paths_t(shortest_paths)
		{}

	void add_composite_shortest_path(size_t i, size_t j) {
		composite_paths.emplace_back(i, j);
	}

	bool is_single_path(size_t index) const {
		return composite_paths[index].second == INVALID_VERTEX;
	}
};

// Structure containing vectors of absolute and relative deviations for all
// edges (i,k) of the given graph that have triangle decomposition
// (i,k) = (i,j) \circ (j,k).
// That is, (i,k) is triangle decomposable if there is a vertex j, such that (i,j) and (j,k) are edges of the graph
struct edge_weight_deviations_t {
	std::vector<double> absolute_deviations;
	std::vector<double> relative_deviations;
};

struct edge_weight_rel_deviation_t {
	double rel_deviation;
	size_t source;  // source node i
	size_t target;  // target node k
	size_t best_intermediate;   // best intermediate node j

	explicit edge_weight_rel_deviation_t(double rel_deviation,
										 size_t source,
										 size_t target,
										 size_t best_intermediate)
		: rel_deviation(rel_deviation),
		  source(source),
		  target(target),
		  best_intermediate(best_intermediate)
		{}
};

struct gradient_trajectory_t {
	std::vector<size_t> path;///< The vertices forming the trajectory
	bool ends_at_critical;   ///< Whether the trajectory ends at local maximum or minimum
	bool ends_at_lmax;   ///< Whether the trajectory ends at local maximum
	double quality_metric;   ///< Quality metric of the trajectory (monotonicity * adjusted rate)
	double total_change; ///< Total function value change along the trajectory
};

/**
 * @brief Structure to hold a prediction result with a flag indicating exclusion
 */
struct prediction_result_t {
	double value;   ///< Predicted value
	bool is_excluded;   ///< Flag indicating if prediction could not be made
};

struct set_wgraph_t {

	std::vector<std::set<edge_info_t>> adjacency_list; 	// Core graph structure
	double graph_diameter;
	double max_packing_radius;

	// Default constructor
	set_wgraph_t() : graph_diameter(-1.0), max_packing_radius(-1.0) {}

	explicit set_wgraph_t(
		const std::vector<std::vector<int>>& adj_list,
		const std::vector<std::vector<double>>& weight_list
		);

	explicit set_wgraph_t(
		const iknn_graph_t& iknn_graph
		);

	// Constructor that initializes an empty graph with n_vertices
	explicit set_wgraph_t(size_t n_vertices)
		: adjacency_list(n_vertices),
		  graph_diameter(-1.0),
		  max_packing_radius(-1.0)
		{}

	// Core graph operations
	void print(const std::string& name,
			   bool split,
			   size_t shift
		) const;

	size_t num_vertices() const {
		return adjacency_list.size();
	}

	/**
	 * @brief Compute the median edge length in the graph.
	 */
	double compute_median_edge_length() const;

	void compute_graph_diameter();

	// ----------------------------------------------------------------
	//
	// precompute edge weights members
	//
	// ----------------------------------------------------------------
	/**
	 * @brief Ensures edge weights are computed and available for use
	 */
	void ensure_edge_weights_computed() const {
		if (!edge_weights_computed) {
			precompute_edge_weights();
		}
	}

	/**
	 * @brief Precomputes all edge weights in the graph for efficient lookup
	 */
 	void precompute_edge_weights() const;

	/**
	 * @brief Invalidates cached edge weights
	 */
	void invalidate_edge_weights() {
		edge_weights_computed = false;
		edge_weights.clear();
	}

	// ------- other ----------

	/**
	 * @brief Creates a subgraph containing only the specified vertices
	 */
	set_wgraph_t create_subgraph(
		const std::vector<size_t>& vertices
		) const;

	/**
	 * @brief Counts the number of connected components in the graph
	 */
	size_t count_connected_components() const;

	/**
	 * @brief Get all connected components in the graph
	 */
	std::vector<std::vector<size_t>> get_connected_components() const;

	void add_edge(
		size_t v1,
		size_t v2,
		double weight
		);

	// ----------------------------------------------------------------
	//
	// geometric edge pruning functions
	//
	// ----------------------------------------------------------------
	/**
	 * @brief Compute the length of the shortest path between two vertices while excluding a specific edge
	 */
	double bidirectional_dijkstra_excluding_edge(
		size_t source,
		size_t target
		) const;

	/**
	 * pure bidirectional Dijkstra algorithm
	 */
	double bidirectional_dijkstra(
		size_t source,
		size_t target
		) const;


	/**
	 * @brief Compute statistics for potential edge pruning based on geometric criteria
	 */
	edge_pruning_stats_t compute_edge_pruning_stats(
		double threshold_percentile = 0.5
		) const;

	/**
	 * @brief Prune edges geometrically based on alternative path ratio
	 */
	set_wgraph_t prune_edges_geometrically(
		double max_ratio_threshold = 1.2,
		double threshold_percentile = 0.5
		) const;

	set_wgraph_t prune_long_edges(double threshold_percentile = 0.5) const;

	edge_weight_deviations_t compute_edge_weight_deviations() const;
	std::vector<edge_weight_rel_deviation_t> compute_edge_weight_rel_deviations() const;

	bool is_composite_path_geodesic(
		size_t i,
		size_t j,
		const shortest_paths_t& shortest_paths
		) const;

	// ----------------------------------------------------------------
	//
	// graph gradient flow related functions
	//
	// ----------------------------------------------------------------
	monotonic_reachability_map_t compute_monotonic_reachability_map(
		size_t ref_vertex,
		const std::vector<double>& y,
		double radius,
		bool ascending
		) const;

	path_t reconstruct_monotonic_path(
		const monotonic_reachability_map_t& map,
		size_t target_vertex
		) const;

	std::pair<size_t, double> find_best_gradient_vertex(
		const monotonic_reachability_map_t& map,
		double min_distance,
		size_t min_path_length,
		bool ascending
		) const;

	double calculate_path_evenness(
		const std::vector<double>& edge_lengths
		) const;

	double calculate_monotonicity_index(
		double total_change,
		double cumulative_absolute_changes
		) const;

	std::vector<double> compute_weight_percentiles(
		const std::vector<double>& probs
		) const;

	std::vector<double> extract_edge_lengths(
		const path_t& path
		) const;

	std::vector<std::tuple<vertex_shortest_path_info_t, double, double>> evaluate_paths(
		const std::vector<vertex_shortest_path_info_t>& paths,
		size_t current_vertex,
		const std::vector<double>& y
		) const;


	std::vector<local_extremum_t> detect_local_minima(
		const std::vector<double>& y,
		double max_radius,
		size_t min_neighborhood_size
		) const;

	std::vector<local_extremum_t> detect_local_maxima(
		const std::vector<double>& y,
		double max_radius,
		size_t min_neighborhood_size
		) const;

	std::vector<local_extremum_t> detect_local_extrema(
		const std::vector<double>& y,
		double max_radius,
		size_t min_neighborhood_size,
		bool detect_maxima
		) const;

	basin_t find_local_extremum_geodesic_basin(
		size_t vertex,
		const std::vector<double>& y,
		bool detect_maxima
		) const;

	basin_t find_local_extremum_bfs_basin(
		size_t vertex,
		const std::vector<double>& y,
		bool detect_maxima
		) const;

	basin_cx_t create_basin_cx(
		const std::vector<double>& y
		) const;

	std::unordered_map<size_t, basin_t> compute_basins(
		const std::vector<double>& y
		) const;

	harmonic_smoother_t harmonic_smoother(
		std::vector<double>& harmonic_predictions,
		std::unordered_set<size_t>& region_vertices,
		int max_iterations = 100,
		double tolerance = 1e-6,
		int record_frequency = 1,
		size_t stability_window = 3,
		double stability_threshold = 0.05
		) const;

	void perform_harmonic_repair(
		std::vector<double>& harmonic_predictions,
		const basin_t& absorbing_basin,
		const basin_t& absorbed_basin
		) const;

	void perform_harmonic_smoothing(
		std::vector<double>& harmonic_predictions,
		std::unordered_set<size_t>& region_vertices,
		int max_iterations,
		double tolerance
		) const;

	double basin_cx_difference(
		const std::unordered_map<size_t, basin_t>& basins1,
		const std::unordered_map<size_t, basin_t>& basins2
		) const;

	basin_t find_gflow_basin(
		size_t vertex,
		const std::vector<double>& y,
		size_t min_basin_size,
		size_t min_path_size,
		double q_edge_thld,
		bool detect_maxima
		) const;

	std::pair<std::vector<basin_t>, std::vector<basin_t>>
	find_gflow_basins(
		const std::vector<double>& y,
		size_t min_basin_size,
		size_t min_path_size,
		double q_edge_thld
		) const;

	basin_t find_local_extremum(
		size_t vertex,
		const std::vector<double>& y,
		size_t min_basin_size,
		bool detect_maxima
		) const;

	std::pair<std::vector<size_t>, std::vector<size_t>> find_nbr_extrema(
		const std::vector<double>& y
		) const;

	std::pair<std::vector<basin_t>, std::vector<basin_t>>
	find_local_extrema(
		const std::vector<double>& y,
		size_t min_basin_size
		) const;

	std::vector<int> watershed_edge_weighted(
		const std::vector<double>& y,
		size_t min_basin_size
		) const;

	// -------- other --------

	gradient_flow_t compute_gradient_flow(
		std::vector<double>& y,
		std::vector<double>& scale,
		double quantile_scale_thld
		) const;

	void remove_edge(size_t v1, size_t v2);

	std::vector<vertex_shortest_path_info_t> get_vertex_shortest_paths(
		const reachability_map_t& reachability_map
		) const;

	std::vector<vertex_path_t> reconstruct_graph_paths(
		const reachability_map_t& reachability_map
		) const;

	reachability_map_t compute_graph_reachability_map(
		size_t ref_vertex,
		double radius
		) const;

	circular_param_result_t parameterize_circular_graph(
		bool use_edge_lengths = true
		) const;

	circular_param_result_t parameterize_circular_graph_with_reference(
		size_t reference_vertex,
		bool use_edge_lengths = true
		) const;

	// ----------------------------------------------------------------
	//
	// graph_spectral_lowess related functions
	//
	// ----------------------------------------------------------------
	/**
	 * @brief Original spectral LOWESS implementation (no model averaging)
	 */
	graph_spectral_lowess_t graph_spectral_lowess(
		const std::vector<double>& y,
		size_t n_evectors,
		// bw parameters
		size_t n_bws,
		bool log_grid,
		double min_bw_factor,
		double max_bw_factor,
		// kernel parameters
		double dist_normalization_factor,
		size_t kernel_type,
		// other
		double precision,
		size_t n_cleveland_iterations,
		bool verbose
		) const;

	/**
	 * @brief Non-adaptive (no model averaging) spectral LOWESS implementation
	 */
	nada_graph_spectral_lowess_t nada_graph_spectral_lowess(
		const std::vector<double>& y,
		size_t n_evectors,
		// bw parameters
		size_t n_bws,
		bool log_grid,
		double min_bw_factor,
		double max_bw_factor,
		// kernel parameters
		double dist_normalization_factor,
		size_t kernel_type,
		// other
		double precision,
		size_t n_cleveland_iterations,
		bool verbose
		) const;

	std::vector<double> graph_deg0_lowess(
		const std::vector<double>& y,
		double bandwidth,
		size_t kernel_type,
		double dist_normalization_factor,
		bool verbose
		) const;

	graph_deg0_lowess_cv_t graph_deg0_lowess_cv(
		const std::vector<double>& y,
		double min_bw_factor,
		double max_bw_factor,
		size_t n_bws,
		bool log_grid,
		size_t kernel_type,
		double dist_normalization_factor,
		bool use_uniform_weights,
		size_t n_folds,
		bool with_bw_predictions,
		double precision,
		bool verbose
		);

	graph_deg0_lowess_cv_mat_t graph_deg0_lowess_cv_mat(
		const std::vector<std::vector<double>>& Y,
		double min_bw_factor,
		double max_bw_factor,
		size_t n_bws,
		bool log_grid,
		size_t kernel_type,
		double dist_normalization_factor,
		bool use_uniform_weights,
		size_t n_folds,
		bool with_bw_predictions,
		double precision,
		bool verbose
		);

	// Matrix version of graph_spectral_lowess for multiple response variables
	graph_spectral_lowess_mat_t graph_spectral_lowess_mat(
		const std::vector<std::vector<double>>& Y,
		size_t n_evectors,
		// bw parameters
		size_t n_bws,
		bool log_grid,
		double min_bw_factor,
		double max_bw_factor,
		// kernel parameters
		double dist_normalization_factor,
		size_t kernel_type,
		// other
		double precision,
		size_t n_cleveland_iterations,
		bool with_errors,
		bool with_scale,
		bool verbose
		) const;

	/**
	 * @brief Perform degree-0 LOWESS with buffer zone cross-validation for bandwidth selection
	 */
	graph_deg0_lowess_buffer_cv_t graph_deg0_lowess_buffer_cv(
		const std::vector<double>& y,
		double min_bw_factor,
		double max_bw_factor,
		size_t n_bws,
		bool log_grid,
		size_t kernel_type,
		double dist_normalization_factor,
		bool use_uniform_weights,
		size_t buffer_hops,
		bool auto_buffer_hops,
		size_t n_folds,
		bool with_bw_predictions,
		double precision,
		bool verbose
		);

	/**
	 * @brief Perform graph kernel smoothing with buffer zone cross-validation for bandwidth selection
	 */
	graph_kernel_smoother_t graph_kernel_smoother(
		const std::vector<double>& y,
		double min_bw_factor,
		double max_bw_factor,
		size_t n_bws,
		bool log_grid,
		size_t vertex_hbhd_min_size,
		size_t kernel_type,
		double dist_normalization_factor,
		bool use_uniform_weights,
		size_t buffer_hops,
		bool auto_buffer_hops,
		size_t n_folds,
		bool with_bw_predictions,
		double precision,
		bool verbose
		);

	/**
	 * @brief Spectral-based graph smoother with global bandwidth selection and optional diagnostics
	 */
	graph_bw_adaptive_spectral_smoother_t graph_bw_adaptive_spectral_smoother(
		const std::vector<double>& y,
		size_t n_evectors,
		double min_bw_factor,
		double max_bw_factor,
		size_t n_bws,
		bool log_grid,
		size_t kernel_type,
		double dist_normalization_factor,
		double precision,
		bool use_global_bw_grid,
		bool with_bw_predictions,
		bool with_vertex_bw_errors,
		bool verbose
	);

	klaps_low_pass_smoother_t klaps_low_pass_smoother(
		const std::vector<double>& y,
		size_t n_evectors_to_compute,
		size_t min_num_eigenvectors,
		size_t max_num_eigenvectors,
		double tau_factor,
		double radius_factor,
		size_t laplacian_power,
		size_t n_candidates,
		bool   log_grid,
		double energy_threshold,
		bool   with_k_predictions,
		bool   verbose
		) const;

	klaps_spectral_filter_t
	klaps_spectral_filter(
		const std::vector<double>& y,
		filter_type_t filter_type,
		size_t n_evectors_to_compute,
		double tau_factor,
		double radius_factor,
		size_t laplacian_power,
		size_t n_candidates,
		bool   log_grid,
		bool   with_t_predictions,
		bool   verbose
		) const;

	/**
	 * @brief Computes the eigenvectors of the graph Laplacian.
	 */
	Eigen::MatrixXd compute_graph_laplacian_eigenvectors(
		size_t n_evectors,
		bool verbose
		) const;

	Eigen::MatrixXd compute_graph_shifted_kernel_laplacian_eigenvectors(
		size_t n_evectors,
		double tau,
		size_t k,
		bool verbose
		) const;

	std::pair<Eigen::VectorXd, Eigen::MatrixXd>
	compute_graph_shifted_kernel_laplacian_spectrum(
		size_t n_evectors,
		double tau,
		double radius_factor,
		size_t laplacian_power,
		bool verbose
		) const;

	/**
	 * @brief Returns both eigenvalues and eigenvectors
	 */
	std::pair<Eigen::VectorXd, Eigen::MatrixXd>
	compute_graph_laplacian_spectrum(
		size_t n_evectors,
		bool   verbose
	) const;

	/**
	 * @brief Compute total squared curvature using the random‑walk normalized Laplacian on an unweighted graph.
	 */
	double get_total_sq_curvature(
		const std::vector<double>& predictions
		) const;

	/**
	 * @brief Compute sum of squared normalized curvature (Laplacian) with additive regularization.
	 */
 	double get_total_sq_normalized_curvature(
		const std::vector<double>& predictions,
		double delta
		) const;



	/**
	 * @brief Create a buffer zone around test vertices up to a specified hop distance
	 */
	std::unordered_set<size_t> create_buffer_zone(
		const std::vector<size_t>& test_vertices,
		size_t buffer_hops
		);

	/**
	 * @brief Creates spatially stratified cross-validation folds for graph-based data
	 */
	std::vector<std::vector<size_t>> create_spatially_stratified_folds(
		size_t n_folds
		);

	/**
	 * @brief Predicts the value at a test vertex using training vertices outside a buffer zone
	 */
	double predict_test_vertex_with_buffer(
		size_t test_vertex,
		size_t buffer_hops,
		const std::vector<double>& y,
		const std::vector<size_t>& current_fold,
		double dist_normalization_factor,
		bool use_uniform_weights
		);

	/**
	 * @brief Predict values for test vertices while respecting buffer zones
	 */
	std::vector<prediction_result_t> predict_with_buffer_zone(
		double bandwidth,
		const std::vector<double>& weights,
		const std::vector<double>& y,
		const std::vector<size_t>& test_vertices,
		const std::unordered_set<size_t>& buffer_zone,
		double dist_normalization_factor,
		bool use_uniform_weights
		);

	/**
	 * @brief Generate predictions for all vertices using the optimal bandwidth
	 */
	std::vector<double> predict_all_with_optimal_bandwidth(
		double bandwidth,
		const std::vector<double>& y,
		double dist_normalization_factor,
		bool use_uniform_weights,
		double lower_bound,
		double upper_bound,
		size_t domain_min_size
		);

	/**
	 * @brief Determine the optimal buffer hop distance based on spatial autocorrelation
	 */
	size_t determine_optimal_buffer_hops(
		const std::vector<double>& y,
		bool verbose
		);

	/**
	 * @brief Calculate Moran's I spatial autocorrelation statistic
	 */
	double calculate_morans_i(
		const std::vector<double>& y,
		size_t hop_distance
		);

	/**
	 * @brief Model-averaged spectral LOWESS implementation
	 */
	graph_spectral_lowess_t graph_spectral_ma_lowess(
		const std::vector<double>& y,
		size_t n_evectors,
		// bw parameters
		size_t n_bws,
		bool log_grid,
		double min_bw_factor,
		double max_bw_factor,
		// kernel parameters
		double dist_normalization_factor,
		size_t kernel_type,
		// model parameters
		double model_blending_coef,
		// other
		double precision,
		bool verbose
		) const;

	/**
	 * @brief Find vertices within a specified radius of a reference vertex
	 */
	std::unordered_map<size_t, double> find_vertices_within_radius(
		size_t vertex,
		double radius
		) const;

	/**
	 * @brief Find minimum radius that includes a specified number of vertices
	 */
	double find_minimum_radius_for_domain_min_size(
		size_t vertex,
		double lower_bound,
		double upper_bound,
		size_t domain_min_size,
		double precision
		) const;

	double find_min_radius_for_neighbor_count(
		std::unordered_map<size_t, double>& ngbr_dist_map,
		size_t vertex_hbhd_min_size
		) const;

	std::pair<std::vector<std::pair<size_t, double>>, double>
	get_sorted_vertices_and_min_radius(
		const std::unordered_map<size_t, double>& ngbr_dist_map,
		size_t vertex_hbhd_min_size
		) const;

	std::pair<std::vector<std::pair<size_t, double>>, double>
	get_sorted_vertices_and_min_radius(
		const std::unordered_map<size_t, double>& ngbr_dist_map,
		size_t vertex_hbhd_min_size,
		std::unordered_set<size_t>& training_set
		) const;

	size_t get_lextr_count(
		std::vector<double>& y
		) const;

	// ----------------------------------------------------------------


	std::pair<size_t, double> get_vertex_eccentricity(
		size_t start_vertex
		) const;

	std::vector<size_t> create_maximal_packing(
		double radius,
		size_t start_vertex
		) const;

	std::vector<size_t> create_maximal_packing(
		size_t grid_size,
		size_t max_iterations,
		double precision
		);

	std::vector<size_t> find_boundary_vertices_outside_radius(
		size_t start,
		double radius,
		explored_tracker_t& explored_tracker
		) const;

	std::pair<size_t, double> find_first_vertex_outside_radius(
		size_t start,
		double radius
		) const;

	std::pair<size_t, double> find_first_vertex_outside_radius(
		size_t start,
		double radius,
		explored_tracker_t& explored_tracker
		) const;

	double compute_shortest_path_distance(
		size_t from,
		size_t to
		) const;

	std::unordered_map<size_t,double>
	compute_shortest_path_distances(
		size_t from,
		const std::unordered_set<size_t>& to_set
		) const;

	void trace_exploration_path(
		size_t from_vertex,
		size_t to_vertex,
		double radius
		) const;

	//
	// AGEMALO helper functions
	//
	shortest_paths_t find_graph_paths_within_radius(
		size_t start,
		double radius
		) const;

	gray_xyw_t get_xyw_along_path(
		const std::vector<double>& y,
		path_t& path,
		double dist_normalization_factor
		) const;

	/**
	 * @brief Find the minimum bandwidth that ensures at least one geodesic ray has the minimum required size
	 */
	double find_minimum_bandwidth(
		size_t grid_vertex,
		double lower_bound,
		double upper_bound,
		size_t min_path_size,
		double precision
		) const;

	/**
	 * @brief Check if any path in the shortest_paths structure has the minimum required size
	 */
	bool has_sufficient_path_size(
		const shortest_paths_t& paths,
		size_t min_path_size
		) const;


	/**
	 * @brief Cluster similar geodesic rays and select representatives
	 */
	std::vector<path_t> cluster_and_select_rays(
		size_t grid_vertex,
		const shortest_paths_t& shortest_paths,
		const std::vector<double>& y,
		double directional_threshold,
		double jaccard_thld,
		size_t min_path_size,
		double dist_normalization_factor
		) const;

	/**
	 * @brief Determine if two paths explore similar directions from a reference vertex
	 */
	bool paths_share_direction(
		const path_t& path_i,
		const path_t& path_j,
		size_t ref_vertex,
		double directional_threshold,
		double jaccard_thld
		) const;

	/**
	 * @brief Select the best path from a cluster based on weighted correlation R-squared
	 */
	size_t select_best_path(
		const std::vector<size_t>& cluster,
		const shortest_paths_t& shortest_paths,
		const std::vector<double>& y,
		size_t min_path_size,
		double dist_normalization_factor
		) const;

	/**
	 * @brief Find the optimal bandwidth that minimizes mean prediction error using binary search
	 */
	std::pair<double, std::vector<ext_ulm_t>> find_optimal_bandwidth(
		size_t grid_vertex,
		const std::vector<double>& y,
		size_t min_path_size,
		double min_bw,
		double max_bw,
		double dist_normalization_factor,
		bool y_binary,
		double tolerance,
		double precision = 0.001,
		bool verbose = false,
		const std::optional<std::vector<double>>& weights = std::nullopt
		) const;

	/**
	 * @brief Find the optimal bandwidth that minimizes mean prediction error using candidate bw grid
	 */
	opt_bw_t find_optimal_bandwidth_over_grid(
		size_t grid_vertex,
		const std::vector<double>& y,
		size_t min_path_size,
		double min_bw,
		double max_bw,
		size_t n_bws,
		bool log_grid,
		double dist_normalization_factor,
		bool y_binary,
		double tolerance,
		double precision = 1e-6,
		bool verbose = false,
		const std::optional<std::vector<double>>& weights = std::nullopt
		) const;

	/**
	 * @brief Evaluate the prediction error for models at a specific bandwidth
	 */
	std::pair<double, std::vector<ext_ulm_t>> evaluate_bandwidth(
		size_t grid_vertex,
		const std::vector<double>& y,
		size_t min_path_size,
		double bw,
		double min_bw,
		double dist_normalization_factor,
		bool y_binary,
		double tolerance,
		const std::optional<std::vector<double>>& weights = std::nullopt
		) const;


private:
	// Cache for computation - updated to use shortest_paths_t
	mutable std::unordered_map<size_t, shortest_paths_t> paths_cache;
	mutable std::unordered_set<size_t> unprocessed_vertices;
	mutable edge_weights_t edge_weights; // using edge_weights_t = std::unordered_map<std::pair<size_t,size_t>, double, size_t_pair_hash_t>;
	mutable bool edge_weights_computed = false;

	// Helper methods for gradient flow computation
	std::optional<std::pair<size_t, bool>> check_local_extremum(
		size_t vertex,
		const shortest_paths_t& paths_result,
		const std::vector<double>& y
		) const;

	gradient_trajectory_t construct_trajectory(
		size_t start,
		bool ascending_init,
		const std::vector<double>& scale,
		const std::vector<double>& y,
		const std::unordered_map<size_t, bool>& extrema_map,
		double long_edge_lower_thld,
		double long_edge_upper_thld
		) const;

	// Helper method to find connected components in a pro-cell
	std::vector<std::vector<size_t>> find_connected_components(
		const std::vector<size_t>& procell_vertices
		) const;
};

#endif   // GRAPH_HPP__
