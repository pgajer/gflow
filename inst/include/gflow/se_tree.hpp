#ifndef SE_TREE_HPP
#define SE_TREE_HPP

#include <cstddef>
#include <vector>
#include <unordered_set>
#include <unordered_map>

/**
 * @brief Node in a Spurious Extrema tree
 */
struct se_node_t {
    std::size_t vertex;                      ///< Vertex index of this extremum
    bool is_maximum;                    ///< True if maximum, false if minimum
    bool is_spurious;                   ///< True if spurious, false if non-spurious (leaf)
    std::size_t parent;                      ///< Parent node vertex (SIZE_MAX if root)
    std::vector<std::size_t> children;       ///< Child node vertices

    // Basin information for this node
    std::vector<std::size_t> basin_vertices; ///< Vertices in this extremum's basin
    std::vector<std::size_t> basin_boundary; ///< Boundary of this extremum's basin
};

/**
 * @brief Classification of SE tree based on terminal (non-spurious) extrema
 */
enum class se_tree_class_t {
    BOTH_TYPES = 0,     ///< Tree has both NSMIN and NSMAX terminals
    ONLY_MIN = 1,       ///< Tree has only NSMIN terminals
    ONLY_MAX = 2,       ///< Tree has only NSMAX terminals
    NO_TERMINALS = 3    ///< Tree has no non-spurious terminals (isolated SE cluster)
};

/**
 * @brief Result of SE tree computation for a single spurious extremum
 *
 * The SE tree captures the hierarchical structure of spurious extrema
 * connected through their basins of attraction. Starting from a root
 * spurious extremum, the tree alternates between minima and maxima,
 * terminating at non-spurious extrema (leaves).
 */
struct se_tree_t {
    std::size_t root_vertex;                             ///< Root spurious extremum
    bool root_is_maximum;                           ///< Type of root extremum
    std::unordered_map<std::size_t, se_node_t> nodes;    ///< All nodes indexed by vertex
    std::vector<std::size_t> ns_min_terminals;           ///< Non-spurious minimum terminals
    std::vector<std::size_t> ns_max_terminals;           ///< Non-spurious maximum terminals
    se_tree_class_t classification;                 ///< Tree classification

    // Precomputed support for harmonic repair
    std::vector<std::size_t> hr_support_vertices;        ///< Union of basins for HR
    std::vector<std::size_t> hr_support_boundary;        ///< Boundary of HR support region
};

#endif // SE_TREE_HPP
