#include "cpp_utils.hpp"

// Convert from map to vector representation
std::vector<double> map_to_vector(
    const std::unordered_map<size_t, double>& map,
    double default_value) {
    // Handle empty map case
    if (map.empty()) {
        return {};
    }

    // Find the maximum key in the map
    size_t max_key = 0;
    for (const auto& [key, _] : map) {
        if (key > max_key) {
            max_key = key;
        }
    }

    // Create vector with size max_key + 1 to accommodate all indices
    std::vector<double> vec(max_key + 1, default_value);

    // Fill in the values
    for (const auto& [idx, val] : map) {
        vec[idx] = val;
    }

    return vec;
}

// Convert from vector to map representation (skipping default values)
std::unordered_map<size_t, double> vector_to_map(
    const std::vector<double>& vec,
    double default_value) {
    std::unordered_map<size_t, double> map;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] != default_value) {
            map[i] = vec[i];
        }
    }
    return map;
}
