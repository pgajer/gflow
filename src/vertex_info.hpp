#ifndef VERTEX_INFO_H_
#define VERTEX_INFO_H_

/**
 * used for sorting shortest paths staring at a given vertex by the distance from that vertex
 */
struct vertex_info_t {
	size_t vertex;
	double distance;
};

/**
 * @brief Stores information about an endpoint vertex and its complete path from a reference vertex
 */
struct vertex_shortest_path_info_t {
	size_t vertex;            ///< The endpoint vertex
	double distance;          ///< Distance from reference vertex to this endpoint
	std::vector<size_t> path; ///< The complete shortest path from reference vertex to endpoint
};

#endif // VERTEX_INFO_H_
