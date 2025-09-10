#if 0

struct eigen_ulogit_t {
    Eigen::VectorXd beta;
    std::vector<double> predictions;
    std::vector<double> loocv_brier_errors;
    std::vector<std::string> warnings;
    bool converged;
    int iterations;
};

// Helper function to compute orthogonal polynomials
Eigen::MatrixXd compute_orthogonal_basis(const double* x, int n, bool fit_quadratic) {
    int p = fit_quadratic ? 3 : 2;
    Eigen::MatrixXd basis(n, p);

    // First basis is constant
    basis.col(0).setOnes();

    // Compute mean and standard deviation of x for centering and scaling
    double mean_x = 0.0;
    double var_x = 0.0;
    for (int i = 0; i < n; ++i) {
        mean_x += x[i];
    }
    mean_x /= n;

    for (int i = 0; i < n; ++i) {
        var_x += (x[i] - mean_x) * (x[i] - mean_x);
    }
    var_x /= (n - 1);
    double sd_x = std::sqrt(var_x);

    // Second basis is linear term (centered and scaled)
    for (int i = 0; i < n; ++i) {
        basis(i, 1) = (x[i] - mean_x) / sd_x;
    }

    if (fit_quadratic) {
        // Third basis is quadratic term (orthogonalized)
        for (int i = 0; i < n; ++i) {
            double z = (x[i] - mean_x) / sd_x;
            basis(i, 2) = (z * z - 1.0) / std::sqrt(2.0);
        }
    }

    return basis;
}

eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Compute orthogonal polynomial basis
    Eigen::MatrixXd X_ortho = compute_orthogonal_basis(x, n, fit_quadratic);

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);

    // Rest of the initialization code remains the same...

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Compute current probabilities using orthogonal basis
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        for (int i = 0; i < n; ++i) {
            // Compute linear predictor using orthogonal basis
            double eta = X_ortho.row(i).dot(result.beta);

            // Rest of the probability computation remains the same...
        }

        // The weighted least squares step now uses X_ortho directly
        Eigen::MatrixXd X_w = weights.asDiagonal() * X_ortho;

        // Rest of the IWLS algorithm remains the same...
    }

    // Final predictions using orthogonal basis
    for (int i = 0; i < n; ++i) {
        double eta = X_ortho.row(i).dot(result.beta);
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    // Standard Rf_error computation would also use X_ortho instead of raw polynomials
    if (with_errors) {
        // Use X_ortho for computing the hat matrix...
    }

    return result;
}



// Using orthogonal polynomials would be beneficial in this logistic regression implementation, especially when fit_quadratic is true.

// The benefits of using orthogonal polynomials in this case include:

// 1. Better numerical stability: Orthogonal polynomials reduce multicollinearity between the linear and quadratic terms, which can be especially problematic in logistic regression.

// 2. Improved conditioning: The design matrix will have better condition number, making the IWLS (Iteratively Weighted Least Squares) solutions more stable.

// 3. Easier interpretation: The coefficients for each polynomial term are independent of each other, making it easier to assess their individual contributions.

// 4. More efficient convergence: The orthogonalization typically leads to faster convergence of the IWLS algorithm.

// The key changes I made in the implementation:

// 1. Added a `compute_orthogonal_basis` function that creates an orthogonal polynomial basis using:
//    - Constant term (1's)
//    - Centered and scaled linear term
//    - Orthogonalized quadratic term (when fit_quadratic is true)

// 2. Modified the main fitting function to use this orthogonal basis instead of raw polynomials.

// 3. The basis computation uses proper centering and scaling to ensure numerical stability.

// Would you like me to add any additional features or modify the implementation further?

#endif


#if 0

// Using orthogonal polynomials in logistic regression is more complex than in linear regression because:

// The loss function is non-linear (log-likelihood instead of squared Rf_error)
// The IWLS algorithm uses iteratively updated weights
// The hat matrix computation is more involved

// However, we can still benefit from orthogonalization by:

// Using orthogonal basis functions for the design matrix
// Updating the orthogonalization with respect to the current weights in each iteration

// Here's how I would modify the code:

struct orthogonal_logistic_fit_t {
    std::vector<double> predictions;
    std::vector<double> loocv_brier_errors;
    Eigen::VectorXd beta;
    bool converged;
    int iterations;
    std::vector<std::string> warnings;
    bool fit_quadratic;

    // Store orthogonalization information
    struct orthogonal_basis_t {
        std::vector<double> norms;
        std::vector<std::vector<double>> proj_coeffs;
    } basis;

    std::vector<double> predict(const std::vector<double>& x) const {
        if (x.empty()) return std::vector<double>();

        int n = x.size();
        int p = fit_quadratic ? 3 : 2;
        std::vector<std::vector<double>> basis_values(p, std::vector<double>(n));

        // Compute orthogonal basis values
        for (int i = 0; i < n; ++i) {
            basis_values[0][i] = 1.0 / basis.norms[0];
        }

        for (int j = 1; j < p; ++j) {
            for (int i = 0; i < n; ++i) {
                basis_values[j][i] = x[i] * basis_values[j-1][i];
                for (int k = 0; k < j; ++k) {
                    basis_values[j][i] -= basis.proj_coeffs[j][k] * basis_values[k][i];
                }
                basis_values[j][i] /= basis.norms[j];
            }
        }

        // Compute predictions
        std::vector<double> predictions(n);
        for (int i = 0; i < n; ++i) {
            double eta = 0.0;
            for (int j = 0; j < p; ++j) {
                eta += beta(j) * basis_values[j][i];
            }
            predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }

        return predictions;
    }
};

orthogonal_logistic_fit_t eigen_orthogonal_logit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    orthogonal_logistic_fit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    result.fit_quadratic = fit_quadratic;
    result.basis.norms.resize(p);
    result.basis.proj_coeffs.resize(p);

    // Initialize as before...

    double prev_deviance = std::numeric_limits<double>::max();

    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);

        // Compute current probabilities and working weights...

        // Create orthogonal basis with respect to current working weights
        Eigen::MatrixXd X(n, p);
        X.col(0).setOnes();
        result.basis.norms[0] = std::sqrt(
            (X.col(0).array().square() * Eigen::Map<const Eigen::ArrayXd>(working_w.data(), n)).sum()
        );
        X.col(0) /= result.basis.norms[0];

        for (int j = 1; j < p; ++j) {
            result.basis.proj_coeffs[j].resize(j);

            // Compute next polynomial
            Eigen::VectorXd new_col = X.col(j-1).cwiseProduct(Eigen::Map<const Eigen::VectorXd>(x, n));

            // Gram-Schmidt with respect to working weights
            for (int k = 0; k < j; ++k) {
                double numerator = (new_col.array() * X.col(k).array() *
                    Eigen::Map<const Eigen::ArrayXd>(working_w.data(), n)).sum();
                double denominator = (X.col(k).array().square() *
                    Eigen::Map<const Eigen::ArrayXd>(working_w.data(), n)).sum();
                result.basis.proj_coeffs[j][k] = numerator / denominator;
                new_col -= result.basis.proj_coeffs[j][k] * X.col(k);
            }

            result.basis.norms[j] = std::sqrt(
                (new_col.array().square() * Eigen::Map<const Eigen::ArrayXd>(working_w.data(), n)).sum()
            );
            X.col(j) = new_col / result.basis.norms[j];
        }

        // Weight the matrices using working weights
        Eigen::VectorXd sqrt_w = working_w.array().sqrt();
        Eigen::MatrixXd X_w = sqrt_w.asDiagonal() * X;
        Eigen::VectorXd z_w = sqrt_w.cwiseProduct(Eigen::Map<const Eigen::VectorXd>(working_y.data(), n));

        // Since X is orthogonal with respect to working weights,
        // XtWX is diagonal plus ridge penalty
        Eigen::VectorXd XtWX_diag = X_w.cwiseProduct(X_w).colwise().sum();
        XtWX_diag.array() += ridge_lambda;

        // Solve system directly
        Eigen::VectorXd new_beta = (X_w.transpose() * z_w).cwiseQuotient(XtWX_diag);

        // Continue with step halving and convergence checks...
    }

    // Modified hat matrix computation for LOOCV if with_errors is true...

    return result;
}

// Key modifications and benefits:

// Orthogonalization is performed with respect to current working weights in each iteration
// The normal equations become diagonal (before ridge penalty) due to orthogonality
// Improved numerical stability in the IWLS iterations
// More efficient hat matrix computation for LOOCV
// Better handling of ill-conditioning through orthogonalization

// The main tradeoff is increased computational cost per iteration due to the orthogonalization process. However, this may be offset by:

// Faster convergence due to better conditioning
// More stable iterations requiring less step halving
// Simpler matrix operations due to diagonal structure

// Would you like me to elaborate on any of these aspects or show more details of the implementation?

#endif
