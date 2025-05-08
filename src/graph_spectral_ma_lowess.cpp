#include <SymEigsSolver.h>             // For eigenvalue computation
#include <MatOp/DenseSymMatProd.h>
#include <MatOp/SparseSymMatProd.h>
#include <Eigen/Dense>

#include <execution>                   // For std::execution::seq/par
#include <thread>                      // For std::thread::hardware_concurrency
#include <algorithm>                   // For std::min_element
#include <numeric>                     // For std::accumulate
#include <map>                         // For std::map

// for debugging
#include <filesystem>
#include <fstream>
#include "cpp_utils.hpp"               // For debugging and elapsed.time

#include "graph_spectral_lowess.hpp"   // For graph_spectral_lowess_t
#include "bandwidth_utils.hpp"         // For get_candidate_bws
#include "kernels.h"                   // For kernel functions
#include "error_utils.h"               // For REPORT_ERROR
#include "mlm.hpp"                     // For lm_t structure
#include "set_wgraph.hpp"              // For set_wgraph_t
#include "mabilo.hpp"                  // For uwmabilo()
#include "cpp_stats_utils.hpp"         // For running_window_average()


/**
 * @brief Creates a spectral embedding for a set of vertices using Laplacian eigenvectors
 */
Eigen::MatrixXd create_spectral_embedding(
    const std::map<size_t, double>& vertex_map,
    const Eigen::MatrixXd& eigenvectors,
    size_t n_evectors);

/**
 * @brief Model-averaged spectral-based locally weighted regression for graph data
 *
 * @details This function computes local linear approximations of the response variable
 * at each vertex using spectral embedding of local neighborhoods, then averages multiple
 * models to produce robust predictions. The algorithm:
 *
 * 1. Computes Laplacian eigenvectors for dimension reduction
 * 2. For each grid vertex:
 *    - Identifies neighboring vertices within varying bandwidths
 *    - Creates spectral embeddings of these local neighborhoods
 *    - Fits weighted linear models at multiple bandwidths
 *    - Records model information for all vertices in each model's domain
 * 3. For each vertex in the graph:
 *    - Combines predictions from all models containing the vertex
 *    - Weights models based on proximity and prediction error
 *    - Produces model-averaged predictions, errors, and scale estimates
 *
 * This approach improves robustness to bandwidth selection and provides better
 * predictions in regions with complex geometry or varying density.
 *
 * @param y Response values at each vertex in the graph
 * @param n_evectors Number of eigenvectors to use for the spectral embedding
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing for bandwidth grid; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances in kernel weight calculations
 * @param kernel_type Type of kernel function to use for weighting (e.g., Gaussian, triangular)
 * @param precision Precision threshold for binary search and numerical comparisons
 * @param model_blending_coef Coefficient between 0 and 1 for blending position-based and error-based weights
 * @param verbose Whether to print progress information
 *
 * @return graph_spectral_lowess_t Structure containing:
 *         - predictions: Model-averaged predictions at each vertex
 *         - errors: Estimated prediction errors
 *         - scale: Local bandwidth/scale parameter for each vertex
 */
/**
 * @brief Model-averaged spectral-based locally weighted regression for graph data
 *
 * @details This function computes local linear approximations of the response variable
 * at each vertex using spectral embedding of local neighborhoods, then averages multiple
 * models to produce robust predictions. The algorithm:
 *
 * 1. Computes Laplacian eigenvectors for dimension reduction
 * 2. For each grid vertex:
 *    - Identifies neighboring vertices within varying bandwidths
 *    - Creates spectral embeddings of these local neighborhoods
 *    - Fits weighted linear models at multiple bandwidths
 *    - Records model information for all vertices in each model's domain
 * 3. For each vertex in the graph:
 *    - Combines predictions from all models containing the vertex
 *    - Weights models based on proximity and prediction error
 *    - Produces model-averaged predictions, errors, and scale estimates
 *
 * This approach improves robustness to bandwidth selection and provides better
 * predictions in regions with complex geometry or varying density.
 *
 * @param y Response values at each vertex in the graph
 * @param n_evectors Number of eigenvectors to use for the spectral embedding
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing for bandwidth grid; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances in kernel weight calculations
 * @param kernel_type Type of kernel function to use for weighting (e.g., Gaussian, triangular)
 * @param precision Precision threshold for binary search and numerical comparisons
 * @param model_blending_coef Coefficient between 0 and 1 for blending position-based and error-based weights
 * @param verbose Whether to print progress information
 *
 * @return graph_spectral_lowess_result_t Structure containing:
 *         - predictions: Model-averaged predictions at each vertex
 *         - errors: Estimated prediction errors
 *         - scale: Local bandwidth/scale parameter for each vertex
 */
graph_spectral_lowess_t set_wgraph_t::graph_spectral_ma_lowess(
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
) const {
    auto total_start_time = std::chrono::steady_clock::now();

    // Initialize kernel function
    initialize_kernel(kernel_type, 1.0);

    // Set Eigen to use available threads for parallel computation
    unsigned int available_threads = std::thread::hardware_concurrency();
    if (available_threads == 0) available_threads = 4; // Fallback if detection fails
    Eigen::setNbThreads(available_threads);

    // Define minimum and maximum bandwidth based on graph diameter
    if (graph_diameter <= 0) {
        REPORT_ERROR("Invalid graph diameter: %f. Must be positive.", graph_diameter);
    }

    double min_bw = min_bw_factor * graph_diameter;
    double max_bw = max_bw_factor * graph_diameter;

    // Minimum number of vertices needed for a robust n_evectors-dimensional model
    size_t domain_min_size = n_evectors + 5; // Adding extra points for stability

    if (verbose) {
        Rprintf("Starting graph_spectral_ma_lowess algorithm\n");
        Rprintf("Number of vertices: %zu\n", adjacency_list.size());
        Rprintf("Graph diameter: %f\n", graph_diameter);
        Rprintf("min_bw: %.4f\n", min_bw);
        Rprintf("max_bw: %.4f\n", max_bw);
        Rprintf("\nNumber of eigenvectors: %zu\n", n_evectors);
        Rprintf("domain_min_size: %zu\n\n", domain_min_size);
        Rprintf("Using %u threads for Eigen operations\n", available_threads);
        Rprintf("\n");
    }

    // Initialize result structure
    size_t n_vertices = adjacency_list.size();
    graph_spectral_lowess_t result;
    result.predictions.resize(n_vertices, 0.0);
    result.errors.resize(n_vertices, 0.0);
    result.scale.resize(n_vertices, 0.0);

    // Construct the Laplacian matrix
    if (verbose) {
        Rprintf("Constructing Laplacian matrix...\n");
    }

    // Create adjacency matrix as a sparse matrix
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;
    triples.reserve(adjacency_list.size() * 2); // Estimate for undirected graph

    for (size_t i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            size_t j = edge.vertex;
            if (i < j) {  // Ensure each edge is added only once
                triples.push_back(Eigen::Triplet<double>(i, j, edge.weight));
                triples.push_back(Eigen::Triplet<double>(j, i, edge.weight));
            }
        }
    }
    A.setFromTriplets(triples.begin(), triples.end());

    // Compute degree matrix
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        double sum = 0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            sum += it.value();
        }
        D.insert(k, k) = sum;
    }

    // Compute Laplacian: L = D - A
    Eigen::SparseMatrix<double> L = D - A;

    // Add small regularization for numerical stability
    for (int i = 0; i < n_vertices; i++) {
        L.coeffRef(i, i) += 1e-8;
    }

    // Compute eigenvalues and eigenvectors of the Laplacian
    if (verbose) {
        Rprintf("Computing Laplacian eigenvectors...\n");
    }

    Spectra::SparseSymMatProd<double> op(L);

    // Ensure we request enough eigenvectors
    int nev = std::min(n_evectors + 5, n_vertices); // Request a few extra
    int ncv = std::min(2 * nev, (int)n_vertices); // Control parameter for the algorithm

    // Ensure nev < ncv
    if (nev >= ncv) {
        ncv = nev + 1;
    }

    // Construct eigen solver to find eigenvalues closest to 0
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    int maxit = 1000; // Increase maximum iterations
    double tol = 1e-10; // Tolerance for convergence
    eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

    if (eigs.info() != Spectra::CompInfo::Successful) {
        REPORT_ERROR("Eigenvalue computation failed. Try increasing max iterations or adjusting tolerance.");
    }

    // Get eigenvectors
    Eigen::MatrixXd eigenvectors = eigs.eigenvectors();

    if (verbose) {
        Rprintf("Eigendecomposition completed. Processing vertices...\n");
    }

    // Define the weight-prediction-error structure similar to agemalo implementation
    struct wpe_t {
        double weight;       // Kernel weight based on distance
        double prediction;   // Model prediction at this vertex
        double error;        // LOOCV error estimate
        double mean_error;   // Mean error of the entire model
        double bw;           // Bandwidth used for the model

        wpe_t(double w, double p, double e, double me, double b)
            : weight(w), prediction(p), error(e), mean_error(me), bw(b) {}
    };

    // Create a vector to store model information for each vertex
    std::vector<std::vector<wpe_t>> vertex_models(n_vertices);

    // Identify the vertices that will serve as model centers
    // In this case, we'll use all vertices as potential centers
    std::vector<size_t> model_centers(n_vertices);
    std::iota(model_centers.begin(), model_centers.end(), 0);

    // Create a progress tracker
    std::atomic<size_t> progress_counter{0};
    const size_t progress_step = std::max(size_t(1), n_vertices / 20);

	// First, process all model centers to create models
    auto model_center_processing_start = std::chrono::steady_clock::now();
    std::vector<std::unordered_map<size_t, std::vector<wpe_t>>> thread_local_results(model_centers.size());

    std::for_each(std::execution::par, model_centers.begin(), model_centers.end(),
                  [&](size_t i) {
                      size_t center_vertex = model_centers[i];
                      try {
                          // Find minimum bandwidth that ensures enough vertices for modeling
                          double min_vertex_bw = find_minimum_radius_for_domain_min_size(
                              center_vertex,
                              min_bw,
                              max_bw,
                              domain_min_size,
                              precision
                          );

                          // Check bandwidth constraints
                          if (min_vertex_bw >= max_bw) {
                              if (verbose) {
                                  Rprintf("Warning: Required minimum bandwidth (%.4f) exceeds maximum bandwidth (%.4f) for vertex %zu. Skipping.\n",
                                         min_vertex_bw, max_bw, center_vertex);
                              }
                              return;
                          }

                          // Generate candidate bandwidths
                          std::vector<double> candidate_bws = get_candidate_bws(
                              min_vertex_bw,
                              max_bw,
                              n_bws,
                              log_grid,
                              precision
                          );

                          // Find all vertices within maximum radius
                          auto vertices_within_radius = find_vertices_within_radius(center_vertex, max_bw);

                          // Process each bandwidth
                          for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {
                              double current_bw = candidate_bws[bw_idx];

                              // Filter vertices by current bandwidth
                              std::map<size_t, double> local_vertex_map;
                              for (const auto& [vertex, distance] : vertices_within_radius) {
                                  if (distance <= current_bw) {
                                      local_vertex_map[vertex] = distance;
                                  }
                              }

                              // Skip if we don't have enough vertices
                              if (local_vertex_map.size() < domain_min_size) {
                                  continue;
                              }

                              // Create embedding using eigenvectors
                              Eigen::MatrixXd embedding = create_spectral_embedding(
                                  local_vertex_map,
                                  eigenvectors,
                                  n_evectors
                              );

                              // Fit weighted linear model
                              lm_t model = fit_linear_model(
                                  embedding,
                                  y,
                                  local_vertex_map,
                                  dist_normalization_factor
                              );


							  for (size_t j = 0; j < model.vertices.size(); ++j) {
								  size_t vertex = model.vertices[j];

								  // Make sure the vertex exists in the map, create it if not
								  if (thread_local_results[i].find(vertex) == thread_local_results[i].end()) {
									  thread_local_results[i][vertex] = std::vector<wpe_t>();
								  }

								  // Validate indices
								  if (j < model.weights.size() && j < model.predictions.size() && j < model.errors.size()) {
									  // Create a temporary object and use push_back instead of emplace_back
									  wpe_t temp_wpe(
										  model.weights[j],       // Weight from kernel function
										  model.predictions[j],   // Prediction at this vertex
										  model.errors[j],        // LOOCV error at this vertex
										  model.mean_error,       // Mean error across the model
										  current_bw              // Bandwidth used
										  );

									  try {
										  thread_local_results[i][vertex].push_back(temp_wpe);
									  } catch (const std::exception& e) {
										  REPORT_ERROR("Exception during push_back for vertex %zu: %s\n", vertex, e.what());
									  } catch (...) {
										  REPORT_ERROR("Unknown exception during push_back for vertex %zu\n", vertex);
									  }
								  } else {
									  REPORT_ERROR("ERROR: index %zu out of bounds for model (vertices: %zu, weights: %zu, predictions: %zu, errors: %zu)\n",
											  j, model.vertices.size(), model.weights.size(), model.predictions.size(), model.errors.size());
								  }
							  }
						  }

                          // Update progress counter
                          if (verbose) {
                              size_t current = ++progress_counter;
                              if (current % progress_step == 0 || current == model_centers.size()) {
                                  double percentage = 100.0 * current / model_centers.size();
                                  Rprintf("\rProcessing model centers: %.1f%% (%zu/%zu)",
                                         percentage, current, model_centers.size());
                                  R_FlushConsole();
                              }
                          }
                      } catch (const std::exception& e) {
                          REPORT_ERROR("Error processing center vertex %zu: %s", center_vertex, e.what());
                      }
                  });

    // Now merge the thread-local results into the global vertex_models
    if (verbose) {
        Rprintf("\nMerging results from multiple threads...\n");
    }

    vertex_models.resize(n_vertices);
    for (size_t i = 0; i < thread_local_results.size(); ++i) {
        for (const auto& [vertex, models] : thread_local_results[i]) {
            // Append the models from this thread to the global vertex_models
            vertex_models[vertex].insert(
                vertex_models[vertex].end(),
                models.begin(),
                models.end()
            );
        }
    }

    if (verbose) {
        Rprintf("\nModel center processing completed in %.2f seconds\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - model_center_processing_start).count());

        // Collect statistics on number of models per vertex
        std::vector<size_t> model_counts(n_vertices);
        for (size_t i = 0; i < n_vertices; ++i) {
            model_counts[i] = vertex_models[i].size();
        }

        // Print summary statistics
        print_vector_quantiles(model_counts, "Models per vertex", {0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0});

        // Print warning for vertices with no models
        size_t zero_model_count = std::count(model_counts.begin(), model_counts.end(), 0);
        if (zero_model_count > 0) {
            REPORT_WARNING("Warning: %zu vertices (%.1f%%) have no models\n",
                           zero_model_count, 100.0 * zero_model_count / n_vertices);
        }
    }

    // Average the models for each vertex
    auto model_averaging_start = std::chrono::steady_clock::now();

    if (verbose) {
        Rprintf("Performing model averaging...\n");
    }

    // Create a vector of vertices to process in parallel
    std::vector<size_t> all_vertices(n_vertices);
    std::iota(all_vertices.begin(), all_vertices.end(), 0);

    // Process each vertex in parallel
    std::for_each(std::execution::par, all_vertices.begin(), all_vertices.end(),
        [&](size_t vertex) {
            const auto& models = vertex_models[vertex];

            // Skip if no models for this vertex
            if (models.empty()) {
                result.predictions[vertex] = std::numeric_limits<double>::quiet_NaN();
                result.errors[vertex] = std::numeric_limits<double>::infinity();
                result.scale[vertex] = max_bw;
                return;
            }

            double prediction_sum = 0.0;
            double weight_sum = 0.0;
            double error_sum = 0.0;
            double scale_sum = 0.0;

            // Average all models that include this vertex
            for (const auto& model : models) {
                // Calculate effective weight based on blending coefficient
                double effective_weight;

                if (model_blending_coef == 0.0) {
                    // Pure position weight
                    effective_weight = model.weight;
                }
                else if (model_blending_coef == 1.0) {
                    // Full mean error influence
                    effective_weight = model.weight * model.mean_error;
                }
                else {
                    // Smooth interpolation between the two approaches
                    effective_weight = (1.0 - model_blending_coef) * model.weight +
                                       model_blending_coef * (model.weight * model.mean_error);
                }

                prediction_sum += effective_weight * model.prediction;
                error_sum += effective_weight * model.error;
                scale_sum += effective_weight * model.bw;
                weight_sum += effective_weight;
            }

            // Compute weighted averages
            if (weight_sum > 1e-10) {
                result.predictions[vertex] = prediction_sum / weight_sum;
                result.errors[vertex] = error_sum / weight_sum;
                result.scale[vertex] = scale_sum / weight_sum;
            } else {
                // Fallback for numerical stability issues
                result.predictions[vertex] = y[vertex];
                result.errors[vertex] = std::numeric_limits<double>::infinity();
                result.scale[vertex] = max_bw;
            }
        }
    );

    if (verbose) {
        elapsed_time(model_averaging_start, "Model averaging completed", true);
        elapsed_time(total_start_time, "graph_spectral_ma_lowess algorithm completed", true);
    }

    return result;
}
