#ifndef ULOGIT_HPP
#define ULOGIT_HPP

#include <Eigen/Dense>
#include <vector>

/**
 * @struct ulogit_eigen_t
 * @brief Structure containing results from Eigen-based logistic regression fit
 *
 * @details This structure holds the fitted model parameters and diagnostic information
 * from logistic regression using Eigen linear algebra library
 */
struct eigen_ulogit_t {
    std::vector<double> predictions; ///< Fitted probabilities for each observation
    std::vector<double> errors;      ///< LOOCV prediction errors (when computed)
    Eigen::VectorXd beta;            ///< Model coefficients (intercept, linear [, quadratic])
    bool converged;                  ///< Whether the fitting algorithm converged
    int iterations;                  ///< Number of iterations used in fitting
};

eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic = false,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8,
    bool with_errors = false);

struct ulogit_t {
    std::vector<double> w;        // weight values for the x_values
    int x_min_index;             // index of smallest x value in dataset
    int x_max_index;             // index of largest x value in dataset

    // Model evaluation
    std::vector<double> predictions;   // predicted probabilities
    std::vector<double> errors;        // Leave-One-Out Cross-Validation log-loss errors

    // debugging info
    int iteration_count;
    bool converged;
};

ulogit_t ulogit(const double* x,
                const double* y,
                const std::vector<double>& w,
                int max_iterations = 100,
                double ridge_lambda = 0.002,
                double max_beta = 100.0,
                double tolerance = 1e-8,
                bool verbose = false);

std::vector<double> ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double max_beta = 100.0,
    double tolerance = 1e-8);

#endif // ULOGIT_HPP
