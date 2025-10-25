#ifndef GRADIENT_BASIN_H_
#define GRADIENT_BASIN_H_

#include <unordered_map>
using std::size_t;

struct gradient_basin_t {
    size_t vertex;                                   // extremum vertex (origin of basin)
    double value;                                    // y value at extremum
    bool is_maximum;                                 // true if descending basin, false if ascending
    size_t hop_idx;                                  // maximum hop distance in basin
    std::unordered_map<size_t, size_t> hop_dist_map; // vertex -> hop distance from origin
    std::unordered_map<size_t, size_t> predecessors; // vertex -> predecessor (for trajectory reconstruction)
    std::unordered_map<size_t, double> y_nbhd_bd_map; // boundary vertices -> y value
    std::vector<size_t> terminal_extrema;            // indices of terminal extrema reachable from origin
};

#endif // GRADIENT_BASIN_H_
