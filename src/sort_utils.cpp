#include "sort_utils.hpp"
#include <algorithm>
#include <numeric>

void sort_by_x_keep_y_and_order(
	const std::vector<double>& x,
	const std::vector<double>& y,
	std::vector<double>& x_sorted,
	std::vector<double>& y_sorted,
	std::vector<std::size_t>& order
	) {
	const std::size_t n = x.size();
	order.resize(n);
	std::iota(order.begin(), order.end(), 0u);

	std::sort(order.begin(), order.end(),
			  [&](std::size_t i, std::size_t j){ return x[i] < x[j]; });

	x_sorted.resize(n);
	y_sorted.resize(n);
	for (std::size_t k = 0; k < n; ++k) {
		x_sorted[k] = x[ order[k] ];
		y_sorted[k] = y[ order[k] ];
	}

	std::vector<std::size_t> inv(n);
	for (std::size_t k = 0; k < n; ++k) {
		inv[ order[k] ] = k;
	}

	order = std::move(inv);
}
