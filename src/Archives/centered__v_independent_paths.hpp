#ifndef CENTERED_PATHS_H_
#define CENTERED_PATHS_H_

#include <vector>

/**
 * @brief Structure representing a path in a graph centered around a specific vertex
 *
 * Stores information about a path including its vertices, the target vertex it's meant
 * to be centered around, total weight (length), and how far the target vertex is from
 * the center of the path.
 */
struct path_t {
    std::vector<int> vertices;    ///< Sequence of vertices in the path
    int target_vertex;            ///< The vertex this path is meant to be centered around
    double total_weight;          ///< Total length of the path
    double center_offset;         ///< Distance of target vertex from path center (0 = perfectly centered)
    std::vector<double> dist_to_target; ///< dist_to_target[i] is the distance from vertices[i] to target_vertex (using edge weights)

    path_t() : target_vertex(-1), total_weight(0), center_offset(INFINITY) {}

    // For priority queue ordering
    bool operator<(const path_t& other) const {
        if (std::abs(total_weight - other.total_weight) > 1e-10) {
            return total_weight < other.total_weight;  // prefer paths closer to target length
        }
        return center_offset > other.center_offset;    // prefer more centered paths
    }
};

#endif // CENTERED_PATHS_H_
