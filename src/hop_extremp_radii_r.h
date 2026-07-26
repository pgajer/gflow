#ifndef GFLOW_HOP_EXTREMP_RADII_R_H
#define GFLOW_HOP_EXTREMP_RADII_R_H

#include <Rinternals.h>

SEXP S_compute_hop_extremp_radii_batch(
    SEXP s_adj_list,
    SEXP s_edge_densities,
    SEXP s_vertex_densities,
    SEXP s_candidates,
    SEXP s_y,
    SEXP s_p_threshold,
    SEXP s_detect_maxima,
    SEXP s_max_hop
);

#endif
