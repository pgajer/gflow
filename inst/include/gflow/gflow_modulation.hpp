#ifndef GFLOW_MODULATION_HPP
#define GFLOW_MODULATION_HPP

#include <string>

// ============================================================================
// Modulation Types
// ============================================================================

/**
 * @brief Gradient flow modulation strategies
 *
 * Controls how the "steepest" direction is determined when following
 * gradient trajectories. The score for moving from vertex v to neighbor u
 * is computed differently under each modulation:
 *
 * - NONE: score = Δy = y[u] - y[v]
 *   Standard gradient flow following raw function differences.
 *
 * - DENSITY: score = ρ(u) · Δy
 *   Prefers directions toward higher-density regions, useful when
 *   the function should follow population concentrations.
 *
 * - EDGELEN: score = dl([v,u]) · Δy
 *   Multiplies by the edge length density weight dl, which is the
 *   max-normalized KDE of edge lengths evaluated at this edge's length.
 *   Edges near the mode of the length distribution receive weight ≈ 1,
 *   while very short or very long edges are down-weighted. This prevents
 *   "basin jumping" through atypical edges.
 *
 * - DENSITY_EDGELEN: score = ρ(u) · dl([v,u]) · Δy
 *   Combined modulation incorporating both effects.
 */
enum class gflow_modulation_t {
    NONE = 0,           ///< Standard gradient flow: Δy
    DENSITY = 1,        ///< Density-modulated: ρ(u) · Δy
    EDGELEN = 2,        ///< Edge-length-modulated: dl([v,u]) · Δy
    DENSITY_EDGELEN = 3 ///< Combined: ρ(u) · dl([v,u]) · Δy
};

/**
 * @brief Convert modulation enum to string for display
 */
inline std::string gflow_modulation_to_string(gflow_modulation_t mod) {
    switch (mod) {
        case gflow_modulation_t::NONE: return "NONE";
        case gflow_modulation_t::DENSITY: return "DENSITY";
        case gflow_modulation_t::EDGELEN: return "EDGELEN";
        case gflow_modulation_t::DENSITY_EDGELEN: return "DENSITY_EDGELEN";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Convert string to modulation enum
 */
inline gflow_modulation_t string_to_gflow_modulation(const std::string& s) {
    if (s == "NONE" || s == "none") return gflow_modulation_t::NONE;
    if (s == "DENSITY" || s == "density") return gflow_modulation_t::DENSITY;
    if (s == "EDGELEN" || s == "edgelen") return gflow_modulation_t::EDGELEN;
    if (s == "DENSITY_EDGELEN" || s == "density_edgelen") return gflow_modulation_t::DENSITY_EDGELEN;
    return gflow_modulation_t::NONE;  // Default
}

#endif // GFLOW_MODULATION_HPP
