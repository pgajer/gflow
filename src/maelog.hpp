#ifndef MAELOG_H_
#define MAELOG_H_

#include <vector>

struct point_t {
// int model_i; // model index
// int i; // index
// double x; // x-value
// double y; // y-value
	double w; // weight
	double p; // predicted value

	double beta1; // beta[1] coeff
	double beta2; // beta[2] coeff

	double e; // LOOCV error at the given point
	//double deviance_error;
	double brier_error;
	//double abs_error;
};

struct local_logit_t {
	std::vector<double> preds;

	std::vector<double> beta1s;
	std::vector<double> beta2s;

	// std::vector<double> errs;
	//std::vector<double> deviance_errors;
	std::vector<double> brier_errors;
	//std::vector<double> abs_errors;

};

/**
 * @brief Results structure for model-averaged bandwidth logistic regression
 *
 * This structure contains both the results of the model fitting process and
 * the parameters used in the fit. When bandwidth selection is performed
 * (pilot_bandwidth <= 0), it includes the candidate bandwidths and their
 * associated errors.
 */
struct maelog_t {
	// Grid search results
	std::vector<double> candidate_bandwidths;             ///< Grid of bandwidths tested during optimization
	std::vector<std::vector<double>> bw_predictions;      ///< Predictions for each bandwidth in LOOCV estimation (cv_folds = 0)

	// Mean errors and optimal indices
	std::vector<double> mean_brier_errors;                      ///< Mean Brier error for each candidate bandwidth
	int opt_brier_bw_idx;                                       ///< Index of bandwidth with minimal mean Brier error

	// Model coefficients only for the case of pilot_bandwidth > 0
	std::vector<double> beta1s;                           ///< beta[1] coefficient of the local linear or quadratic model
	std::vector<double> beta2s;                           ///< beta[2] coefficient of the local quadratic model

	// Input parameters (unchanged)
	bool fit_quadratic;      ///< Whether quadratic term was included in local models
	double pilot_bandwidth;  ///< Fixed bandwidth if > 0, otherwise bandwidth is selected by CV
	int kernel_type;         ///< Type of kernel function used for local weighting
	int cv_folds;            ///< Number of CV folds (0 for LOOCV approximation)
	double min_bw_factor;    ///< Lower bound factor for bandwidth grid relative to h_rot
	double max_bw_factor;    ///< Upper bound factor for bandwidth grid relative to h_rot
	int max_iterations;      ///< Maximum iterations for logistic regression fitting
	double ridge_lambda;     ///< Ridge regularization parameter
	double tolerance;
};

maelog_t maelog(
	const std::vector<double>& x,
	const std::vector<double>& y,
	bool fit_quadratic,
	double pilot_bandwidth,
	int kernel_type,
	int min_points = 6,
	int cv_folds = 0,
	int n_bws = 50,
	double min_bw_factor = 0.05,
	double max_bw_factor = 0.9,
	int max_iterations = 100,
	double ridge_lambda = 1e-6,
	double tolerance = 1e-8,
	bool with_errors = false,
	bool with_bw_preditions = true);

#endif // MAELOG_H_
