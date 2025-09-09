#ifndef UNIFIED_GRAPH_H_
#define UNIFIED_GRAPH_H_


class weighted_graph {
private:
	std::vector<std::vector<int>> adj_list;
	std::vector<std::vector<double>> weight_list;

public:
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
		// Returns a view that combines adj_list[vertex] and weight_list[vertex]
		// without copying data
		return std::views::zip(adj_list[vertex], weight_list[vertex])
			| std::views::transform([](const auto& p) {
				return edge_t{p.first, p.second};
			});
	}

	// Analysis methods can be organized into categories using namespaces
	namespace algorithms {
		// These free functions take the graph as parameter
		// but are still logically grouped
		bfs_result_t compute_tau_bfs(
			const weighted_graph& g,
			int start,
			double tau);

		gradient_flow_t compute_gradient_flow(
			const weighted_graph& g,
			const std::vector<double>& y,
			double tau);
	}
};

// Implementation can be split across multiple files for better organization
namespace weighted_graph::algorithms {
	namespace detail {
		// Internal implementation details
	}

	bfs_result_t compute_tau_bfs(
		const weighted_graph& g,
		int start,
		double tau) {
		// Implementation
	}
}

struct unified_graph {
	// Internal representation type
	enum class repr_type {
		VECTOR_PAIR,
		SET_EDGES
	};

	// Constructor specifying the representation type
	explicit unified_graph(repr_type type) : representation(type) {}

	// Edge info structure (used for both representations)
	struct edge_info_t {
		int vertex;
		double weight;

		bool operator<(const edge_info_t& other) const {
			return vertex < other.vertex;
		}
	};

	// Add edge method that handles both representations
	void add_edge(int from, int to, double weight) {
		ensure_size(std::max(from, to) + 1);

		if (representation == repr_type::VECTOR_PAIR) {
			adj_list[from].push_back(to);
			weight_list[from].push_back(weight);
		} else {
			set_list[from].insert({to, weight});
		}
	}

	// Get neighbors and weights - returns a vector of edge_info_t for consistency
	std::vector<edge_info_t> get_edges(int vertex) const {
		std::vector<edge_info_t> result;

		if (representation == repr_type::VECTOR_PAIR) {
			result.reserve(adj_list[vertex].size());
			for (size_t i = 0; i < adj_list[vertex].size(); ++i) {
				result.push_back({adj_list[vertex][i], weight_list[vertex][i]});
			}
		} else {
			result.reserve(set_list[vertex].size());
			for (const auto& edge : set_list[vertex]) {
				result.push_back(edge);
			}
		}

		return result;
	}

	// Get number of vertices
	size_t size() const {
		if (representation == repr_type::VECTOR_PAIR) {
			return adj_list.size();
		}
		return set_list.size();
	}

private:
	repr_type representation;

	// Storage for vector pair representation
	std::vector<std::vector<int>> adj_list;
	std::vector<std::vector<double>> weight_list;

	// Storage for set representation
	std::vector<std::set<edge_info_t>> set_list;

	// Helper to ensure the graph has enough space
	void ensure_size(size_t new_size) {
		if (representation == repr_type::VECTOR_PAIR) {
			if (adj_list.size() < new_size) {
				adj_list.resize(new_size);
				weight_list.resize(new_size);
			}
		} else {
			if (set_list.size() < new_size) {
				set_list.resize(new_size);
			}
		}
	}
};


#endif   // UNIFIED_GRAPH_H_
