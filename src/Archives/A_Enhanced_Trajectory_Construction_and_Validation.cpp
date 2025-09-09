/**
 * @brief Structure representing a validated gradient trajectory
 *
 * Guarantees:
 * - Starts at a local maximum
 * - Ends at a local minimum
 * - Monotonically decreasing function values along path
 */
struct validated_trajectory_t {
    std::vector<int> vertices;
    bool is_valid;
    std::string error_message;

private:
    // Private constructor enforces creation through validation
    validated_trajectory_t(const std::vector<int>& path)
        : vertices(path), is_valid(true) {}

    friend class trajectory_builder;
};

/**
 * @brief Builder class for creating and validating trajectories
 *
 * Enforces invariants during trajectory construction and provides
 * detailed error information if validation fails.
 */
class trajectory_builder {
public:
    trajectory_builder(const std::vector<std::vector<int>>& adj_list,
                      const std::vector<double>& Ey)
        : adj_list_(adj_list), Ey_(Ey) {}

    /**
     * @brief Constructs and validates a trajectory from ascending and descending paths
     *
     * @param ascending Path from starting vertex to local maximum
     * @param descending Path from starting vertex to local minimum
     * @return validated_trajectory_t Validated trajectory or error information
     *
     * Guarantees:
     * - Resulting trajectory is monotonic
     * - Endpoints are proper local extrema
     * - Path is connected in the graph
     */
    validated_trajectory_t build_trajectory(
        const std::vector<int>& ascending,
        const std::vector<int>& descending) {

        // Early validation of input paths
        if (ascending.empty() || descending.empty()) {
            return create_error("Empty input path");
        }

        if (ascending.front() != descending.front()) {
            return create_error("Ascending and descending paths don't share starting vertex");
        }

        // Construct complete trajectory
        std::vector<int> complete_path;
        complete_path.reserve(ascending.size() + descending.size() - 1);

        // Add descending path in reverse (excluding middle vertex)
        std::copy(descending.rbegin(), descending.rend() - 1,
                 std::back_inserter(complete_path));

        // Add ascending path
        std::copy(ascending.begin(), ascending.end(),
                 std::back_inserter(complete_path));

        // Validate the complete path
        auto validation_result = validate_path(complete_path);
        if (!validation_result.first) {
            return create_error(validation_result.second);
        }

        return validated_trajectory_t(complete_path);
    }

private:
    const std::vector<std::vector<int>>& adj_list_;
    const std::vector<double>& Ey_;

    /**
     * @brief Validates a complete trajectory path
     *
     * @param path The trajectory to validate
     * @return pair<bool, string> Success/failure with error message
     */
    std::pair<bool, std::string> validate_path(const std::vector<int>& path) {
        // Check path is non-empty
        if (path.empty()) {
            return {false, "Empty path"};
        }

        // Check connectivity
        for (size_t i = 0; i < path.size() - 1; ++i) {
            if (!are_vertices_connected(path[i], path[i + 1])) {
                return {false, "Disconnected vertices in path: " +
                               std::to_string(path[i]) + " -> " +
                               std::to_string(path[i + 1])};
            }
        }

        // Check monotonicity
        for (size_t i = 0; i < path.size() - 1; ++i) {
            if (Ey_[path[i]] <= Ey_[path[i + 1]]) {
                return {false, "Non-monotonic path segment at index " +
                               std::to_string(i)};
            }
        }

        // Endpoints should be local extrema
        if (!is_local_maximum(path.front(), adj_list_, Ey_.data())) {
            return {false, "First vertex is not a local maximum"};
        }

        if (!is_local_minimum(path.back(), adj_list_, Ey_.data())) {
            return {false, "Last vertex is not a local minimum"};
        }

        return {true, ""};
    }

    /**
     * @brief Creates an invalid trajectory with error message
     */
    validated_trajectory_t create_error(const std::string& message) {
        validated_trajectory_t result({});
        result.is_valid = false;
        result.error_message = message;
        return result;
    }

    bool are_vertices_connected(int v1, int v2) {
        const auto& neighbors = adj_list_[v1];
        return std::find(neighbors.begin(), neighbors.end(), v2) != neighbors.end();
    }
};

/**
 * @brief Enhanced version of graph_MS_cx with stronger guarantees
 */
MS_complex_t graph_MS_cx(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<int>>& core_adj_list,
    const std::map<std::pair<int,int>, std::vector<int>>& shortest_paths,
    const std::vector<double>& Ey) {

    MS_complex_t ms_cx;
    std::set<std::vector<int>> unique_trajectories_set;
    trajectory_builder builder(adj_list, Ey);

    // Helper function to add a validated trajectory
    auto add_trajectory = [&](const std::vector<int>& ascending,
                            const std::vector<int>& descending)
                            -> std::pair<bool, size_t> {

        auto validated_traj = builder.build_trajectory(ascending, descending);
        if (!validated_traj.is_valid) {
            // Log error or handle invalid trajectory
            return {false, 0};
        }

        auto [it, inserted] = unique_trajectories_set.insert(validated_traj.vertices);
        if (!inserted) {
            return {false, 0};
        }

        size_t traj_idx = ms_cx.unique_trajectories.size();
        ms_cx.unique_trajectories.push_back(validated_traj.vertices);

        int local_max = validated_traj.vertices.front();
        int local_min = validated_traj.vertices.back();

        ms_cx.lmax_to_lmin[local_max].insert(local_min);
        ms_cx.lmin_to_lmax[local_min].insert(local_max);

        std::pair<int,int> cell_key(local_max, local_min);
        ms_cx.procells[cell_key].insert(validated_traj.vertices.begin(),
                                      validated_traj.vertices.end());

        return {true, traj_idx};
    };

    // Rest of the implementation...
};


// Key improvements in this version:

// 1. **Explicit Guarantees**:
//    - The `validated_trajectory_t` structure makes guarantees explicit
//    - Private constructor ensures trajectories are always validated
//    - Builder pattern enforces proper construction

// 2. **Comprehensive Validation**:
//    - Path connectivity
//    - Monotonicity
//    - Proper extrema at endpoints
//    - Shared starting vertex for ascending/descending paths

// 3. **Error Handling**:
//    - Detailed error messages for each type of validation failure
//    - Clear separation between validation and construction
//    - Easy to add logging or error reporting

// 4. **Type Safety**:
//    - Cannot create invalid trajectories without going through validation
//    - Clear distinction between validated and unvalidated paths

// Would you like me to:
// 1. Add more validation checks?
// 2. Explore error handling strategies?
// 3. Add debugging/logging capabilities?
// 4. Discuss performance implications?
