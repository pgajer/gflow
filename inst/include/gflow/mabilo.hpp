#ifndef MABILO_HPP
#define MABILO_HPP

#include <vector>

struct mabilo_t {
    // k values
    int opt_k;     // optimal model averaging k value - the one with the smallest mean LOOCV error
    int opt_k_idx; // optimal model averaging k value index

    // Errors
    std::vector<double> k_mean_errors;   // mean LOOCV squared errors for each k for model averaged predictions
    std::vector<double> smoothed_k_mean_errors;
    std::vector<double> k_mean_true_errors; // mean absolute error between predictions and y_true

    // The best (over all k) model evaluation
    std::vector<double> predictions; // optimal k model averaged predictions

    std::vector<std::vector<double>> k_predictions; // for each k model averaged predictions

    // Bayesian bootstrap creadible intervals
    std::vector<double> bb_predictions; // central location of the Bayesian bootstrap estimates
    std::vector<double> cri_L; // credible intervals lower limit
    std::vector<double> cri_U; // credible intervals upper limit
};

mabilo_t uwmabilo(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const std::vector<double>& y_true,
                  int k_min,
                  int k_max,
                  int distance_kernel,
                  double dist_normalization_factor,
                  double epsilon,
                  bool verbose);

mabilo_t wmabilo(const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& y_true,
                 const std::vector<double>& w,
                 int k_min,
                 int k_max,
                 int distance_kernel,
                 double dist_normalization_factor,
                 double epsilon,
                 bool verbose);

#endif // MABILO_HPP
