#ifndef EDGE_INFO_H_
#define EDGE_INFO_H_

#include <cstddef>  // For size_t

/**
 * @struct edge_info_t
 * @brief Represents an edge in the graph with its target vertex and weight
 *
 * This structure is used within the adjacency lists to store edge information.
 * The comparison operator enables its use in ordered containers like std::set.
 */
struct edge_info_t {
	size_t vertex;   ///< Target vertex of the edge
	double weight;   ///< Weight (length) of the edge

	/**
	 * @brief Comparison operator for ordering in sets
	 * @param other The other edge_info_t to compare with
	 * @return true if this edge's vertex is less than the other's vertex
	 */
	bool operator<(const edge_info_t& other) const {
		return vertex < other.vertex;
	}
};

#endif // EDGE_INFO_H_
