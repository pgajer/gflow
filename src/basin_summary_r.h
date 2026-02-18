/**
 * @file basin_summary_r.h
 * @brief R interface declarations for basin summary acceleration helpers.
 */

#ifndef BASIN_SUMMARY_R_H
#define BASIN_SUMMARY_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_precompute_basin_vertex_metrics(
    SEXP s_adj_list,
    SEXP s_edgelen_list,
    SEXP s_hop_k
);

SEXP S_summary_basins_of_attraction_cpp(
    SEXP s_object,
    SEXP s_adj_list,
    SEXP s_edgelen_list,
    SEXP s_hop_k,
    SEXP s_vertex_metrics
);

#ifdef __cplusplus
}
#endif

#endif // BASIN_SUMMARY_R_H
