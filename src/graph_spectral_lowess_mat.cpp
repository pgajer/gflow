#include <Eigen/Dense>
#include <Spectra/SymEigsSolver.h>             // For eigenvalue computation
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <execution>                   // For std::execution::seq/par
#include <thread>                      // For std::thread::hardware_concurrency
#include <algorithm>                   // For std::min_element
#include <numeric>                     // For std::accumulate
#include <map>                         // For std::map
#include <mutex>

#include "graph_spectral_lowess_mat.hpp"  // For graph_spectral_lowess_mat_t
#include "bandwidth_utils.hpp"         // For get_candidate_bws
#include "kernels.h"                   // For kernel functions
#include "error_utils.h"               // For REPORT_ERROR
#include "mlm.hpp"                     // For lm_t structure
#include "set_wgraph.hpp"              // For set_wgraph_t
#include "mabilo.hpp"                  // For uwmabilo()
#include "cpp_stats_utils.hpp"         // For running_window_average()
#include "progress_utils.hpp"          // For progress_tracker_t

// Forward declaration of the spectral embedding function already defined elsewhere
Eigen::MatrixXd create_spectral_embedding(
    const std::map<size_t, double>& vertex_map,
    const Eigen::MatrixXd& eigenvectors,
    size_t n_evectors);

/**
 * @brief Matrix version of spectral-based locally weighted regression for graph data
 *
 * @details This function computes local linear approximations of multiple response variables
 * at each vertex using spectral embedding of local neighborhoods. The algorithm:
 * 1. Computes the Laplacian eigenvectors once for dimension reduction
 * 2. For each vertex, identifies neighboring vertices within varying bandwidths
 * 3. Creates spectral embeddings of these local neighborhoods
 * 4. Fits weighted linear models for each response variable column and selects optimal bandwidth
 * 5. Produces smoothed predictions, error estimates, and local scale information for each column
 *
 * This matrix version is more efficient than calling the vector version on each column separately,
 * as it computes the graph Laplacian and eigenvectors only once.
 *
 * @param Y Matrix of response values where Y[j][i] is the j-th response variable at the i-th vertex
 * @param n_evectors Number of eigenvectors to use for the spectral embedding
 * @param n_bws Number of candidate bandwidths to evaluate
 * @param log_grid If true, use logarithmic spacing for bandwidth grid; if false, use linear spacing
 * @param min_bw_factor Factor multiplied by graph diameter to determine minimum bandwidth
 * @param max_bw_factor Factor multiplied by graph diameter to determine maximum bandwidth
 * @param dist_normalization_factor Factor for normalizing distances in kernel weight calculations
 * @param kernel_type Type of kernel function to use for weighting (e.g., Gaussian, triangular)
 * @param precision Precision threshold for binary search and numerical comparisons
 * @param n_cleveland_iterations Number of robustness iterations for Cleveland's algorithm
 * @param with_errors Whether to compute and return prediction errors
 * @param with_scale Whether to compute and return bandwidth/scale information
 * @param verbose Whether to print progress information
 *
 * @return graph_spectral_lowess_mat_t Structure containing:
 *         - predictions: Matrix of smoothed values for each response variable at each vertex
 *         - errors: Matrix of estimated prediction errors (if with_errors is true)
 *         - scale: Matrix of local bandwidth/scale parameters (if with_scale is true)
 */
graph_spectral_lowess_mat_t set_wgraph_t::graph_spectral_lowess_mat(
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
) const {
    auto start_time = std::chrono::steady_clock::now();

    // Validate input
    if (Y.empty()) {
        REPORT_ERROR("Input matrix Y is empty");
    }
    
    size_t n_response_vars = Y.size();
    if (n_response_vars == 0) {
        REPORT_ERROR("No response variables provided");
    }
    
    size_t n_vertices = Y[0].size();
    for (size_t j = 1; j < n_response_vars; ++j) {
        if (Y[j].size() != n_vertices) {
            REPORT_ERROR("Inconsistent number of vertices across response variables");
        }
    }
    
    if (n_vertices != adjacency_list.size()) {
        REPORT_ERROR("Number of vertices in Y (%zu) does not match graph size (%zu)",
                   n_vertices, adjacency_list.size());
    }

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

    size_t ncc = count_connected_components();
    if (ncc > 1) {
        REPORT_ERROR("The graph has to have only one connected component!!!\n");
    }

    if (verbose) {
        Rprintf("Starting graph_spectral_lowess_mat() algorithm\n");
        Rprintf("Number of vertices: %zu\n", n_vertices);
        Rprintf("Number of response variables: %zu\n", n_response_vars);
        Rprintf("Graph diameter: %f\n", graph_diameter);
        Rprintf("min_bw: %.4f\n", min_bw);
        Rprintf("max_bw: %.4f\n", max_bw);

        Rprintf("\nNumber of eigenvectors: %zu\n", n_evectors);
        Rprintf("domain_min_size: %zu\n\n", domain_min_size);
        Rprintf("Using %u threads for Eigen operations\n", available_threads);

        Rprintf("\nNumber of connected components: %zu\n", ncc);
        Rprintf("\n");
    }

    // Initialize result structure
    graph_spectral_lowess_mat_t result;
    
    // Initialize predictions for all response variables
    result.predictions.resize(n_response_vars);
    for (size_t j = 0; j < n_response_vars; ++j) {
        result.predictions[j].resize(n_vertices);
    }
    
    // Initialize errors and scale arrays only if needed
    if (with_errors) {
        result.errors.resize(n_response_vars);
        for (size_t j = 0; j < n_response_vars; ++j) {
            result.errors[j].resize(n_vertices);
        }
    }
    
    if (with_scale) {
        result.scale.resize(n_response_vars);
        for (size_t j = 0; j < n_response_vars; ++j) {
            result.scale[j].resize(n_vertices);
        }
    }

    // Step 1: Construct the Laplacian matrix (only done once)
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

    // Step 2: Compute eigenvalues and eigenvectors of the Laplacian (only done once)
    if (verbose) {
        Rprintf("Computing Laplacian eigenvectors...\n");
    }

    Spectra::SparseSymMatProd<double> op(L);

    // Ensure we request enough eigenvectors
    int nev = std::min(static_cast<int>(n_evectors + 5), static_cast<int>(n_vertices)); // Request a few extra
    int ncv = std::min(2 * nev, static_cast<int>(n_vertices)); // Control parameter for the algorithm

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
    Eigen::MatrixXd eigenvectors;

    if (eigs.info() != Spectra::CompInfo::Successful) {
        // if (verbose) {
        // Rprintf("Initial eigenvalue computation failed. Attempting with adjusted parameters...\n");
        // }
        // Define fallback parameters to try
        std::vector<std::pair<int, double>> attempts = {
            {2000, 1e-8},    // More iterations, slightly relaxed tolerance
            {3000, 1e-6},    // Even more iterations, more relaxed tolerance
            {5000, 1e-4}     // Final attempt with very relaxed parameters
        };
        // Try with original ncv first
        bool success = false;
        for (const auto& params : attempts) {
            int adjusted_maxit = params.first;
            double adjusted_tol = params.second;
            // if (verbose) {
            // Rprintf("Trying with original ncv=%d, maxit=%d, tol=%g\n",
            // ncv, adjusted_maxit, adjusted_tol);
            // }
            eigs.init();
            eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
            if (eigs.info() == Spectra::CompInfo::Successful) {
                if (verbose) {
                    Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                            ncv, adjusted_maxit, adjusted_tol);
                }
                eigenvectors = eigs.eigenvectors();  // Add this line to extract eigenvectors
                success = true;
                break;
            }
        }
        // If still not successful, try with increased ncv values
        if (!success) {

            // Define multipliers for ncv
            // std::vector<int> ncv_multipliers = {2, 4, 8};
            // int max_ncv = std::min(16 * nev, (int)n_vertices);

            int max_ncv = std::min((1 << (int)std::log2(n_vertices)), (int)n_vertices);
            std::vector<int> ncv_multipliers = {2, 4, 8, 16, 32, 64};

            for (const auto& multiplier : ncv_multipliers) {
                int adjusted_ncv = std::min(multiplier * ncv, max_ncv);
                // if (verbose) {
                // Rprintf("Trying with increased ncv=%d\n", adjusted_ncv);
                // }
                // Create a new solver with adjusted ncv
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> adjusted_eigs(op, nev, adjusted_ncv);
                for (const auto& params : attempts) {
                    int adjusted_maxit = params.first;
                    double adjusted_tol = params.second;
                    // if (verbose) {
                    // Rprintf("Trying with ncv=%d, maxit=%d, tol=%g\n",
                    // adjusted_ncv, adjusted_maxit, adjusted_tol);
                    // }
                    adjusted_eigs.init();
                    adjusted_eigs.compute(Spectra::SortRule::SmallestAlge, adjusted_maxit, adjusted_tol);
                    if (adjusted_eigs.info() == Spectra::CompInfo::Successful) {
                        if (verbose) {
                            Rprintf("Eigenvalue computation succeeded with adjusted parameters: ncv=%d, maxit=%d, tol=%g\n",
                                    adjusted_ncv, adjusted_maxit, adjusted_tol);
                        }
                        eigenvectors = adjusted_eigs.eigenvectors();
                        success = true;
                        break;
                    }
                }
                if (success) {
                    break;
                }
            }
        }
        // If all attempts failed, report an error
        if (!success) {
            REPORT_ERROR("Eigenvalue computation failed after multiple attempts with adjusted parameters.");
        }
    } else {
        // Get eigenvectors
        eigenvectors = eigs.eigenvectors();
    }

    // if (verbose) {
    //     Eigen::VectorXd eigenvalues = eigs.eigenvalues();
    //     Rprintf("First few eigenvalues:\n");
    //     for (int i = 0; i < std::min(5, static_cast<int>(eigenvalues.size())); i++) {
    //         Rprintf("Eigenvalue %d: %.8f\n", i+1, eigenvalues(i));
    //     }
    // }

    // Set up progress tracking
    progress_tracker_t progress(n_vertices, "Processing vertices", 1);

    // Process each vertex in parallel
    std::vector<size_t> vertices(n_vertices);
    std::iota(vertices.begin(), vertices.end(), 0);

    std::mutex output_mutex; // For synchronizing output operations only
    std::atomic<size_t> completed_vertices(0);
    std::atomic<size_t> last_reported(0);
    const size_t report_interval = 10; // Update every 10 vertices

    // Use std::execution::par_unseq for parallel processing
    std::for_each(std::execution::seq, vertices.begin(), vertices.end(),
                  [&](size_t vertex) {
                      try {
                          // Find minimum bandwidth that ensures enough vertices for modeling
                          double vertex_min_bw = min_bw;

                          // Ensure we have enough vertices for the model
                          vertex_min_bw = find_minimum_radius_for_domain_min_size(
                              vertex,
                              min_bw,
                              max_bw,
                              domain_min_size,
                              precision
                              );

                          // Check bandwidth constraints
                          if (vertex_min_bw >= max_bw) {
                              REPORT_ERROR("Required minimum bandwidth (%.4f) exceeds maximum bandwidth (%.4f) for vertex %zu",
                                           vertex_min_bw, max_bw, vertex);
                          }

                          // Generate candidate bandwidths
                          std::vector<double> candidate_bws = get_candidate_bws(
                              vertex_min_bw,
                              max_bw,
                              n_bws,
                              log_grid,
                              precision
                              );

                          // Find all vertices within maximum radius and sort by distance
                          auto reachable_vertices = find_vertices_within_radius(vertex, max_bw);
                          std::vector<std::pair<size_t, double>> sorted_vertices;
                          sorted_vertices.reserve(reachable_vertices.size());

                          // Convert map to vector for sorting
                          for (const auto& [v, dist] : reachable_vertices) {
                              sorted_vertices.emplace_back(v, dist);
                          }

                          // Sort by distance
                          std::sort(sorted_vertices.begin(), sorted_vertices.end(),
                                    [](const auto& a, const auto& b) { return a.second < b.second; });

                          // Create storage structures for all response variables
                          std::vector<std::vector<double>> bandwidth_errors(n_response_vars,
                                                                            std::vector<double>(candidate_bws.size(), std::numeric_limits<double>::infinity()));
                          std::vector<std::vector<lm_t>> bandwidth_models(n_response_vars,
                                                                          std::vector<lm_t>(candidate_bws.size()));
                          #if 0
                          static std::mutex init_mutex;
                          std::vector<std::vector<double>> bandwidth_errors;
                          std::vector<std::vector<lm_t>> bandwidth_models;
                          {
                              std::lock_guard<std::mutex> lock(init_mutex);
                              // Initialize data structures
                              bandwidth_errors.resize(n_response_vars);
                              bandwidth_models.resize(n_response_vars);
                              for (size_t j = 0; j < n_response_vars; ++j) {
                                  bandwidth_errors[j].resize(candidate_bws.size(), std::numeric_limits<double>::infinity());
                                  bandwidth_models[j].resize(candidate_bws.size());
                              }
                          }
                          #endif

                          // Later when filtering by bandwidth:
                          for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx) {
                              double current_bw = candidate_bws[bw_idx];

                              // Filter vertices by current bandwidth with early termination
                              std::map<size_t, double> local_vertex_map;
                              for (const auto& [v, dist] : sorted_vertices) {
                                  if (dist > current_bw) {
                                      break;  // Early termination - remaining vertices are all farther
                                  }
                                  local_vertex_map[v] = dist;
                              }

                              // Skip if we don't have enough vertices
                              if (local_vertex_map.size() < domain_min_size) {
                                  continue;
                              }

                              // Create embedding using eigenvectors (same for all response variables)
                              Eigen::MatrixXd embedding = create_spectral_embedding(
                                  local_vertex_map,
                                  eigenvectors,
                                  n_evectors
                                  );

                              // Process each response variable
                              for (size_t j = 0; j < n_response_vars; ++j) {
                                  // Fit weighted linear model for this response variable
                                  lm_t model = cleveland_fit_linear_model(
                                      embedding,
                                      Y[j],
                                      local_vertex_map,
                                      dist_normalization_factor,
                                      n_cleveland_iterations
                                      );

                                  // Store model and error
                                  bandwidth_errors[j][bw_idx] = model.mean_error;
                                  bandwidth_models[j][bw_idx] = std::move(model);

                              } // END OF for (size_t j = 0; j < n_response_vars; ++j)
                          } // END OF for (size_t bw_idx = 0; bw_idx < candidate_bws.size(); ++bw_idx)


                          // Process the best bandwidth for each response variable
                          for (size_t j = 0; j < n_response_vars; ++j) {

                              auto response_bandwidth_errors = bandwidth_errors[j];

                              // Find best bandwidth (minimum error)
                              std::vector<double>::iterator min_error_it;
                              bool smooth_bandwidth_errors = true;

                              if (smooth_bandwidth_errors) {
                                  double q_thld = 1.0 / 3.0;
                                  int k = std::max(2, static_cast<int>(q_thld * n_bws));
                                  double epsilon = 1e-10;
                                  std::vector<double> null_vector;
                                  auto response_bandwidth_errors_fit = uwmabilo(candidate_bws,
                                                                       response_bandwidth_errors,
                                                                       null_vector,
                                                                       k,
                                                                       k,
                                                                       kernel_type,
                                                                       dist_normalization_factor,
                                                                       epsilon,
                                                                       false);
                                  response_bandwidth_errors = std::move(response_bandwidth_errors_fit.predictions);
                              }

                              min_error_it = std::min_element(response_bandwidth_errors.begin(), response_bandwidth_errors.end());

                              if (min_error_it != response_bandwidth_errors.end() && std::isfinite(*min_error_it)) {
                                  // Get index of best bandwidth
                                  size_t best_bw_idx = min_error_it - response_bandwidth_errors.begin();

                                  // Store results for this vertex and response variable
                                  const auto& best_model = bandwidth_models[j][best_bw_idx];

                                  if (best_model.vertices.empty()) {
                                      REPORT_ERROR("Empty vertices in best model for vertex %zu, response variable %zu",
                                                   vertex, j);
                                  }

                                  // Find this vertex's prediction in the model
                                  auto it = std::find(best_model.vertices.begin(), best_model.vertices.end(), vertex);

                                  // Store prediction, error, and scale
                                  if (it != best_model.vertices.end()) {
                                      size_t idx = it - best_model.vertices.begin();
                                      result.predictions[j][vertex] = best_model.predictions[idx];

                                      if (with_errors) {
                                          result.errors[j][vertex] = best_model.errors[idx];
                                      }
                                  } else {
                                      REPORT_ERROR("Vertex %zu not found in best model vertices for response variable %zu",
                                                   vertex, j);
                                  }

                                  if (with_scale) {
                                      result.scale[j][vertex] = candidate_bws[best_bw_idx];
                                  }
                              } else {
                                  // No valid models found
                                  result.predictions[j][vertex] = std::numeric_limits<double>::quiet_NaN();

                                  if (with_errors) {
                                      result.errors[j][vertex] = std::numeric_limits<double>::infinity();
                                  }

                                  if (with_scale) {
                                      result.scale[j][vertex] = max_bw;
                                  }

                                  REPORT_WARNING("No valid models found for vertex %zu, response variable %zu",
                                                 vertex, j);
                              }

                          } // END OF for (size_t j = 0; j < n_response_vars; ++j)

                          // Only synchronize output operations
                          if (verbose) {
                              size_t current = completed_vertices.fetch_add(1, std::memory_order_relaxed) + 1;

                              // Try to claim the reporting role for this interval
                              size_t expected = (current / report_interval) * report_interval - report_interval;
                              if (last_reported.compare_exchange_strong(expected, current,
                                                                        std::memory_order_relaxed)) {
                                  // This thread won the race to report progress for this interval
                                  progress.update(current);
                              }
                          }
                      } catch (const std::exception& e) {
                          std::lock_guard<std::mutex> lock(output_mutex);
                          REPORT_ERROR("Error processing vertex %zu: %s", vertex, e.what());
                      }
                  }
        );

    if (verbose) {
        progress.finish();
        double elapsed_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_time).count();
        Rprintf("Completed graph_spectral_lowess_mat in %.2f seconds\n", elapsed_seconds);
    }

    return result;
}
