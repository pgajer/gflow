#ifndef VECT_WGRAPH_H_
#define VECT_WGRAPH_H_

#include <vector>
#include <string>

#include "edge_info.hpp"

struct vect_wgraph_t {
	// Core graph structure
	std::vector<std::vector<edge_info_t>> adjacency_list;

	// Default constructor
	vect_wgraph_t() = default;

	void print(const std::string& name = "",
			   size_t vertex_index_shift = 0) const;
};

#if 0
struct vect_wgraph_t {
	std::vector<std::vector<int>> adj_list;
	std::vector<std::vector<double>> weight_list;

	// Core graph interface
	struct edge_t {
		int target;
		double weight;
	};

	// Efficient accessors that don't copy data
	std::span<const int> neighbors(int vertex) const;
	std::span<const double> weights(int vertex) const;

	// Iterator support for modern C++ integration
	auto edges(int vertex) const {
		return std::views::iota(0u, adj_list[vertex].size())
			| std::views::transform([this, vertex](size_t i) {
				return edge_t{
					adj_list[vertex][i],
					weight_list[vertex][i]
				};
			});
	}
};
#endif

#endif // VECT_WGRAPH_H_
