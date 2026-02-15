#ifndef NADA_GRAPH_SPECTRAL_LOWESS_H_
#define NADA_GRAPH_SPECTRAL_LOWESS_H_

#include <vector>
#include <cstddef>
using std::size_t;

struct nada_graph_spectral_lowess_t {
	std::vector<double> predictions;      ///< predictions[i] is an estimate of E(Y|G) at the i-th vertex
    std::vector<double> errors;           ///< errors[i] is an estimate of the prediction error at the i-th vertex
    std::vector<double> global_bws;       ///< global_bw[i] is the i-th global radius/bandwidth for all vertices
	double opt_bw;                        ///< optimal bandwidth value
	size_t opt_bw_idx;                    ///< the index of opt_bw in global_bws
	std::vector<std::vector<double>> bw_predictions; ///< bw_predictions[bw_idx] are predictions for the bandwidth index 'bw_idx'
};


#endif // NADA_GRAPH_SPECTRAL_LOWESS_H_
