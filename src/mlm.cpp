#include <Eigen/Dense>                 // For Eigen matrix operations

#include <numeric>                     // For std::accumulate
#include <algorithm>                   // For std::max_element
#include <limits>                      // For std::numeric_limits
#include <map>                         // For std::map
#include <optional>                    // For std::optional
//#include <mutex>                       // For std::mutex

#include "kernels.h"                   // For kernel function
#include "mlm.hpp"

// namespace {
//     std::mutex model_fitting_mutex;
// }

/**
 * @brief Safely computes matrix inversion with progressive regularization
 *
 * @param matrix The matrix to invert
 * @param identity_matrix Identity matrix of the same size
 * @return Optional matrix containing the inverse if successful
 */
std::optional<Eigen::MatrixXd> safe_matrix_inverse(
    Eigen::MatrixXd matrix,
    const Eigen::MatrixXd& identity_matrix) {

    // Initial regularization parameter
    double lambda = 1e-8;
    const double max_lambda = 1.0;
    const double lambda_factor = 10.0;

    // Try with progressively stronger regularization
    while (lambda < max_lambda) {
        // Add ridge regularization
        Eigen::MatrixXd regularized = matrix;
        regularized.diagonal().array() += lambda;

        // Try LDLT decomposition (stable for symmetric positive definite matrices)
        Eigen::LDLT<Eigen::MatrixXd> ldlt(regularized);
        if (ldlt.info() == Eigen::Success) {
            try {
                Eigen::MatrixXd inverse = ldlt.solve(identity_matrix);
                // Verify the quality of the inverse
                double relative_error = (regularized * inverse - identity_matrix).norm() / identity_matrix.norm();
                if (relative_error < 1e-5)
                    return inverse;
            } catch (...) {
                // Decomposition succeeded but solve failed, try stronger regularization
            }
        }

        // Try SVD if LDLT failed with reasonable regularization
        if (lambda >= 1e-4) {
            try {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(regularized, Eigen::ComputeThinU | Eigen::ComputeThinV);
                // Use pseudoinverse with appropriate conditioning
                double threshold = 1e-8 * svd.singularValues()(0);
                Eigen::MatrixXd inverse = svd.matrixV() *
                    (svd.singularValues().array() > threshold).select(
                        svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
                    svd.matrixU().transpose();
                return inverse;
            } catch (...) {
                // SVD failed too, increase regularization
            }
        }

        // Increase regularization for next attempt
        lambda *= lambda_factor;
    }

    // All inversion attempts failed
    return std::nullopt;
}

/**
 * @brief Efficiently computes LOOCV errors with multiple fallback strategies
 *
 * @param embedding Design matrix
 * @param XtWX_inv Inverse of X'WX if available
 * @param XtW X'W matrix
 * @param response Response vector
 * @param predictions Fitted values
 * @return Vector of LOOCV squared errors
 */
std::vector<double> compute_loocv_errors(
    const Eigen::MatrixXd& embedding,
    const std::optional<Eigen::MatrixXd>& XtWX_inv,
    const Eigen::MatrixXd& XtW,
    const Eigen::VectorXd& response,
    const Eigen::VectorXd& predictions) {

    size_t n_vertices = embedding.rows();
    std::vector<double> errors(n_vertices);

    // Try efficient hat matrix computation if we have a valid inverse
    if (XtWX_inv) {
        try {
            // Use a more conservative approach - compute hat diagonals one at a time
            // to avoid large intermediate matrix multiplications
            for (size_t i = 0; i < n_vertices; ++i) {
                try {
                    // Extract a single row from embedding
                    Eigen::RowVectorXd xi = embedding.row(i);

                    // Calculate h_ii directly without forming the entire hat matrix
                    // h_ii = xi * (X'WX)^(-1) * X'W * e_i where e_i is the i-th unit vector
                    // This simplifies to: h_ii = xi * (X'WX)^(-1) * (XtW column i)
                    double h_ii = (xi * (*XtWX_inv) * XtW.col(i))(0, 0); // <<--- line (mlm.cpp:105)

                    if (std::abs(1.0 - h_ii) > 1e-10) {
                        double residual = (response(i) - predictions(i)) / (1.0 - h_ii);
                        errors[i] = residual * residual;
                    } else {
                        errors[i] = std::numeric_limits<double>::infinity();
                    }
                } catch (...) {
                    // If calculating individual h_ii fails, use simple residuals
                    double residual = response(i) - predictions(i);
                    errors[i] = residual * residual;
                }
            }
            return errors;
        } catch (...) {
            // If entire approach fails, fall back to approximation
        }
    }

    // Fallback: Use approximate LOOCV errors based on simple residuals
    for (size_t i = 0; i < n_vertices; ++i) {
        double residual = response(i) - predictions(i);
        errors[i] = residual * residual;
    }

    return errors;
}


/**
 * @brief Fits a weighted linear model in spectral embedding space with robust numerical stability
 *
 * @details This function fits a weighted linear model to data in a spectral embedding space,
 * with weights determined by a combination of kernel-weighted distances and optional external weights.
 * It implements progressive regularization and multiple matrix decomposition strategies to ensure
 * numerical stability even with ill-conditioned matrices.
 *
 * The function performs these key steps:
 * 1. Extracts vertex indices and distances from the vertex map
 * 2. Computes kernel weights based on normalized distances from reference vertex
 * 3. Combines kernel weights with external weights if provided
 * 4. Solves the weighted least squares system using progressively stronger regularization
 * 5. Calculates predictions and LOOCV errors with stable block-based computation
 *
 * Numerical stability is ensured through:
 * - Progressive regularization that adaptively increases strength until stability is achieved
 * - Multiple decomposition methods (LDLT, SVD) with proper conditioning
 * - Block-based hat matrix computation to avoid large intermediate matrices
 * - Comprehensive error handling with appropriate fallbacks at each step
 *
 * @param embedding Matrix of spectral coordinates where each row is a vertex and first column is 1.0
 * @param y Vector of response values for all vertices in the graph
 * @param vertex_map Map of vertex indices to their distances from the reference vertex
 * @param dist_normalization_factor Factor to adjust the normalization of distances in kernel calculation
 * @param external_weights Optional vector of weights to be combined with kernel weights
 *
 * @return lm_t Structure containing:
 *   - predictions: Fitted values for each vertex
 *   - errors: Leave-one-out cross-validation squared errors
 *   - mean_error: Average of valid LOOCV errors
 *   - vertices: Indices of vertices included in the model
 *   - weights: Combined weights used for fitting
 *   - intercept: Intercept term of the fitted model
 *   - coefficients: Coefficient vector of the fitted model
 *
 * @note The function handles numerical issues by attempting progressively stronger regularization,
 * starting from 1e-8 and increasing by a factor of 10 until stability is achieved or a maximum value
 * is reached. For severely ill-conditioned matrices, SVD with appropriate conditioning threshold
 * is used as a fallback.
 *
 * @see safe_matrix_inverse For the matrix inversion with progressive regularization
 * @see compute_loocv_errors For the block-based hat matrix computation
 * @see cleveland_fit_linear_model For robust fitting using this function with iterative reweighting
 * @see kernel_fn For the kernel weighting function implementation
 */
/// new
lm_t fit_weighted_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    double dist_normalization_factor,
    const std::optional<std::vector<double>>& external_weights = std::nullopt) {

    // Eigen::setNbThreads(1);
    // // Lock the entire function to ensure thread safety
    // std::lock_guard<std::mutex> lock(model_fitting_mutex);

    lm_t result;

    // Extract vertex indices and distances in the same order as embedding rows
    std::vector<size_t> vertex_indices;
    std::vector<double> distances;

    vertex_indices.reserve(vertex_map.size());
    distances.reserve(vertex_map.size());

    for (const auto& [vertex, distance] : vertex_map) {
        vertex_indices.push_back(vertex);
        distances.push_back(distance);
    }

    size_t n_vertices = vertex_indices.size();

    // Check for empty input
    if (n_vertices == 0) {
        result.mean_error = std::numeric_limits<double>::infinity();
        return result;
    }

    // Normalize distances for kernel function
    double max_dist = *std::max_element(distances.begin(), distances.end());
    if (max_dist <= 0.0) max_dist = 1.0;  // Safety check
    max_dist *= dist_normalization_factor;

    std::vector<double> normalized_distances(distances);
    for (size_t i = 0; i < n_vertices; ++i) {
        normalized_distances[i] /= max_dist;
    }

    // Calculate kernel weights
    std::vector<double> kernel_weights(n_vertices);
    kernel_fn(normalized_distances.data(), n_vertices, kernel_weights.data());

    // Combine kernel weights with external weights if provided
    std::vector<double> combined_weights = kernel_weights;
    if (external_weights && external_weights->size() == n_vertices) {
        for (size_t i = 0; i < n_vertices; ++i) {
            combined_weights[i] *= (*external_weights)[i];
        }
    }

    // Normalize weights to sum to 1
    double total_weight = std::accumulate(combined_weights.begin(), combined_weights.end(), 0.0);
    if (total_weight <= 0.0) {
        // If all weights are zero or negative, use uniform weights
        std::fill(combined_weights.begin(), combined_weights.end(), 1.0 / n_vertices);
    } else {
        for (size_t i = 0; i < n_vertices; ++i) {
            combined_weights[i] /= total_weight;
        }
    }

    // Extract response values for vertices in the same order
    Eigen::VectorXd response(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        response(i) = y[vertex_indices[i]];
    }

    // Convert weights to diagonal matrix
    Eigen::VectorXd weight_vector = Eigen::Map<Eigen::VectorXd>(combined_weights.data(), n_vertices);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W = weight_vector.asDiagonal();

    // Compute matrices for weighted least squares
    Eigen::MatrixXd XtW = embedding.transpose() * W;
    Eigen::MatrixXd XtWX = XtW * embedding;
    Eigen::VectorXd XtWy = XtW * response;

    // Identity matrix for inversion
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(XtWX.rows(), XtWX.cols());

    // Attempt to solve the system robustly
    Eigen::VectorXd beta;
    std::optional<Eigen::MatrixXd> XtWX_inv = safe_matrix_inverse(XtWX, identity);

    if (XtWX_inv) {
        // Use the computed inverse to get coefficients
        beta = (*XtWX_inv) * XtWy;
    } else {
        // Last resort: Heavy regularization and direct solve
        XtWX.diagonal().array() += 0.1;
        beta = XtWX.ldlt().solve(XtWy);
    }

    // Calculate predictions
    Eigen::VectorXd predictions = embedding * beta;
    result.predictions.resize(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        result.predictions[i] = predictions(i);
    }

    // Store vertex indices and weights for later reference
    result.vertices = vertex_indices;
    result.weights = combined_weights;

    // Calculate LOOCV errors using the stable method
    result.errors = compute_loocv_errors(embedding, XtWX_inv, XtW, response, predictions);

    // Calculate mean error
    result.mean_error = 0.0;
    size_t valid_errors = 0;
    for (size_t i = 0; i < n_vertices; ++i) {
        if (std::isfinite(result.errors[i])) {
            result.mean_error += result.errors[i];
            valid_errors++;
        }
    }
    result.mean_error = valid_errors > 0 ? result.mean_error / valid_errors : std::numeric_limits<double>::infinity();

    // Store coefficients
    result.intercept = beta(0);
    result.coefficients = beta;

    return result;
}

/**
 * @brief Fits a robust weighted linear model using Cleveland's iterative reweighting
 *
 * This function extends the standard weighted linear model fitting by incorporating
 * Cleveland's robust iterative reweighting approach. It performs the following steps:
 * 1. Fits an initial weighted linear model
 * 2. Calculates residuals and robust weights
 * 3. Iteratively refits the model with updated weights
 * 4. Calculates final predictions and LOOCV errors
 *
 * The robustness iterations help reduce the influence of outliers by downweighting
 * observations with large residuals using bisquare weighting.
 *
 * @param embedding Matrix where each row is a vertex's coordinates in the embedding space
 * @param y Vector of response values for each vertex
 * @param vertex_map Map of vertex indices to distances from reference vertex
 * @param kernel_type Type of kernel function for distance-based weighting
 * @param dist_normalization_factor Factor for normalizing distances in kernel weights
 * @param n_iterations Number of robustness iterations (default: 3)
 * @param robust_scale Scale for bisquare weight calculation (default: 6.0)
 *
 * @return lm_t Structure containing the fitted robust model and diagnostics
 */
lm_t cleveland_fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    double dist_normalization_factor,
    size_t n_iterations,
    double robust_scale
    ) {

    // Perform initial fit without any robust weights
    lm_t result = fit_weighted_linear_model(embedding, y, vertex_map, dist_normalization_factor);

    // If no iterations requested, return standard fit
    if (n_iterations == 0) {
        return result;
    }

    // Extract vertex indices for convenience
    const std::vector<size_t>& vertices = result.vertices;
    size_t n_vertices = vertices.size();

    // Create working copy of weights for robustness iterations
    std::vector<double> robust_weights(n_vertices, 1.0);  // Start with neutral weights

    // Cleveland's robust fitting with iterative reweighting
    for (size_t iter = 0; iter < n_iterations; ++iter) {
        // Get response values for current vertices
        std::vector<double> y_current(n_vertices);
        for (size_t i = 0; i < n_vertices; ++i) {
            y_current[i] = y[vertices[i]];
        }

        // Calculate residuals from previous fit
        std::vector<double> residuals(n_vertices);
        for (size_t i = 0; i < n_vertices; ++i) {
            residuals[i] = std::abs(y_current[i] - result.predictions[i]);
        }

        // Compute median absolute residual (MAD)
        std::vector<double> sorted_residuals = residuals;
        std::nth_element(sorted_residuals.begin(),
                        sorted_residuals.begin() + sorted_residuals.size()/2,
                        sorted_residuals.end());
        double median_abs_residual = sorted_residuals[sorted_residuals.size()/2];

        // Avoid division by zero
        if (median_abs_residual < 1e-10) {
            break; // No need for more iterations, fit is already good
        }

        // Scale residuals using robust_scale * MAD
        double scale = robust_scale * median_abs_residual;
        for (auto& r : residuals) {
            r /= scale;
        }

        // Apply bisquare weights
        for (size_t i = 0; i < n_vertices; ++i) {
            double u = residuals[i];
            if (u >= 1.0) {
                robust_weights[i] = 0.0;
            } else {
                double tmp = 1.0 - u*u;
                robust_weights[i] = tmp*tmp;
            }
        }

        // Fit model with updated weights
        result = fit_weighted_linear_model(embedding, y, vertex_map, dist_normalization_factor, robust_weights);
    }

    return result;
}

/**
 * @brief Fits a weighted linear model to data in the spectral embedding space
 *
 * @details This function provides compatibility with the original API while
 * delegating to the more robust implementation in fit_weighted_linear_model().
 * It fits a linear model with weights determined by kernel-weighted distances
 * in the original graph.
 *
 * @param embedding Matrix of spectral coordinates where each row is a vertex and first column is 1.0
 * @param y Vector of response values for all vertices in the graph
 * @param vertex_map Map of vertex indices to their distances from the reference vertex
 * @param dist_normalization_factor Factor to adjust the normalization of distances
 *
 * @return lm_t Structure containing fitted model, predictions, and LOOCV errors
 *
 * @see fit_weighted_linear_model For the implementation details and numerical stability features
 */
lm_t fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    double dist_normalization_factor) {

    return fit_weighted_linear_model(embedding, y, vertex_map, dist_normalization_factor);
}
