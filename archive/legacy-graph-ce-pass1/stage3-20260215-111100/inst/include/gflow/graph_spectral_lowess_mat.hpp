#ifndef GRAPH_SPECTRAL_LOWESS_MAT_HPP
#define GRAPH_SPECTRAL_LOWESS_MAT_HPP

#include <vector>

struct graph_spectral_lowess_mat_t {
    std::vector<std::vector<double>> predictions;  ///< predictions[j][i] is an estimate of E(Y[,j]|G) at the i-th vertex
    std::vector<std::vector<double>> errors;       ///< errors[j][i] is an estimate of Y[,j] prediction error at the i-th vertex
    std::vector<std::vector<double>> scale;        ///< scale[j][i] is a local scale (radius/bandwidth) of Y[,j] prediction at the i-th vertex
};

#endif // GRAPH_SPECTRAL_LOWESS_MAT_HPP
