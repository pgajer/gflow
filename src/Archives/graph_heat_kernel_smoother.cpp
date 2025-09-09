#include "graph_heat_kernel_smoother.h"

graph_heat_kernel_smoother_t
set_wgraph_t::graph_heat_kernel_smoother(
	const std::vector<double>& y,
	size_t n_evectors_to_compute,
	size_t laplacian_power,
	double tau_factor,
	size_t n_candidates,
	bool   log_grid,
	bool   with_t_predictions,
	bool   verbose
	) const
{
	graph_heat_kernel_smoother_t result;
	size_t n_vertices = adjacency_list.size();

	// 1) Eigen‐decomposition of (I - L)^k
	double tau = tau_factor * graph_diameter;
	auto [evals, V] = compute_graph_shifted_kernel_laplacian_spectrum(
		n_evectors_to_compute,
		tau,
		laplacian_power,
		verbose
	);

	size_t m = evals.size();

	// 2) Sort eigenpairs ascending by eigenvalue
	std::vector<std::pair<double, Eigen::VectorXd>> pairs;
	pairs.reserve(m);
	for (size_t i = 0; i < m; ++i)
		pairs.emplace_back(evals[i], V.col(i));
	std::sort(pairs.begin(), pairs.end(),
			  [](auto &a, auto &b){ return a.first < b.first; });
	for (size_t i = 0; i < m; ++i) {
		evals[i] = pairs[i].first;
		V.col(i) = pairs[i].second;
	}
	result.evalues.assign(evals.data(), evals.data() + m);
	result.evectors = std::move(V);

	// 3) Graph Fourier coefficients
	Eigen::VectorXd y_ev = Eigen::Map<const Eigen::VectorXd>(y.data(), n_vertices);
	Eigen::VectorXd gft  = result.evectors.transpose() * y_ev;

	// 4) Build diffusion‐time grid [t_min, t_max], excluding t=0
	double eps = 1e-11;
	double eval_max = result.evalues.back();
	double t_max = (eval_max > 0) ? -std::log(eps) / eval_max : 1.0;
	double t_min = t_max / 1e3; // small but positive minimum time

	result.candidate_ts.reserve(n_candidates);
	for (size_t j = 0; j < n_candidates; ++j) {
		double frac = static_cast<double>(j) / static_cast<double>(n_candidates - 1);
		double t = log_grid
			? std::exp(std::log(t_min) + frac * (std::log(t_max) - std::log(t_min)))
			: t_min + frac * (t_max - t_min);
		result.candidate_ts.push_back(t);
	}

	// 5) Loop over t_j, compute GCV and (optionally) store y_{t_j}
	result.gcv_scores.resize(n_candidates);
	if (with_t_predictions)
		result.t_predictions.assign(n_candidates, std::vector<double>(n_vertices));

	for (size_t idx = 0; idx < n_candidates; ++idx) {
		double t = result.candidate_ts[idx];

		// Spectral weights w_i = exp(-t * λ_i), but clip λ_i < 0
		Eigen::ArrayXd w(m);
		for (size_t i = 0; i < m; ++i) {
			if (result.evalues[i] >= 0.0)
				w[i] = std::exp(-t * result.evalues[i]);
			else
				w[i] = 0.0;
		}

		// Smoothed signal y_t
		Eigen::VectorXd y_t = result.evectors * (w * gft.array()).matrix();

		// GCV: ||y - y_t||² / (n - trace(S_t))²
		double trS = w.sum();
		double norm2 = (y_ev - y_t).squaredNorm();
		double denom = static_cast<double>(n_vertices) - trS;
		result.gcv_scores[idx] = norm2 / (denom * denom);

		if (with_t_predictions) {
			auto &buf = result.t_predictions[idx];
			for (size_t v = 0; v < n_vertices; ++v)
				buf[v] = y_t[v];
		} else if (idx == 0) {
			result.predictions.resize(n_vertices);  // preallocate
		}
	}

	// 6) Pick optimal t (lowest GCV)
	result.opt_t_idx = std::min_element(
		result.gcv_scores.begin(), result.gcv_scores.end()
	) - result.gcv_scores.begin();

	// 7) Fill final prediction
	if (with_t_predictions) {
		result.predictions = result.t_predictions[result.opt_t_idx];
	} else {
		double t = result.candidate_ts[result.opt_t_idx];
		Eigen::ArrayXd w(m);
		for (size_t i = 0; i < m; ++i) {
			if (result.evalues[i] >= 0.0)
				w[i] = std::exp(-t * result.evalues[i]);
			else
				w[i] = 0.0;
		}
		Eigen::VectorXd y_t = result.evectors * (w * gft.array()).matrix();
		for (size_t v = 0; v < n_vertices; ++v)
			result.predictions[v] = y_t[v];
	}

	return result;
}
