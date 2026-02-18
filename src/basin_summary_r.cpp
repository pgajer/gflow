#include "SEXP_cpp_conversion_utils.hpp"

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

struct vertex_metrics_t {
    std::vector<double> mean_nbrs_dist;
    std::vector<double> mean_hopk_dist;
    std::vector<int> deg;
    std::vector<double> p_mean_nbrs_dist;
    std::vector<double> p_mean_hopk_dist;
    std::vector<double> p_deg;
    int hop_k = 0;
};

struct basin_row_t {
    int vertex = NA_INTEGER; // 1-based
    double value = NA_REAL;
    double rel_value = NA_REAL;
    int hop_idx = NA_INTEGER;
    int basin_size = 0;
    double p_mean_nbrs_dist = NA_REAL;
    double p_mean_hopk_dist = NA_REAL;
    double deg = NA_REAL;
    double p_deg = NA_REAL;
};

static SEXP get_list_element(SEXP list, const char* name) {
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (Rf_isNull(names) || TYPEOF(names) != STRSXP) {
        return R_NilValue;
    }
    const int n = LENGTH(list);
    for (int i = 0; i < n; ++i) {
        if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

static int get_nrows_from_dim(SEXP x) {
    SEXP dim = Rf_getAttrib(x, R_DimSymbol);
    if (!Rf_isNull(dim) && LENGTH(dim) >= 1) {
        return INTEGER(dim)[0];
    }
    return LENGTH(x);
}

static void validate_graph_inputs(
    const std::vector<std::vector<int>>& adj,
    const std::vector<std::vector<double>>& edgelen
) {
    if (adj.size() != edgelen.size()) {
        Rf_error("adj.list and edgelen.list must have the same length");
    }
    const size_t n = adj.size();
    for (size_t i = 0; i < n; ++i) {
        if (adj[i].size() != edgelen[i].size()) {
            Rf_error("Mismatch at vertex %zu: adj.list has %zu neighbors but edgelen.list has %zu lengths",
                     i + 1, adj[i].size(), edgelen[i].size());
        }
    }
}

static vertex_metrics_t compute_vertex_metrics_cpp(
    const std::vector<std::vector<int>>& adj,
    const std::vector<std::vector<double>>& edgelen,
    int hop_k
) {
    const size_t n = adj.size();
    vertex_metrics_t out;
    out.mean_nbrs_dist.resize(n, R_PosInf);
    out.mean_hopk_dist.resize(n, R_PosInf);
    out.deg.resize(n, 0);
    out.p_mean_nbrs_dist.resize(n, NA_REAL);
    out.p_mean_hopk_dist.resize(n, NA_REAL);
    out.p_deg.resize(n, NA_REAL);
    out.hop_k = hop_k;

    if (hop_k < 1) {
        Rf_error("hop.k must be a positive integer");
    }

    // Per-vertex immediate neighbor metrics.
    for (size_t v = 0; v < n; ++v) {
        const size_t d = edgelen[v].size();
        out.deg[v] = static_cast<int>(d);
        if (d > 0) {
            double sum = std::accumulate(edgelen[v].begin(), edgelen[v].end(), 0.0);
            out.mean_nbrs_dist[v] = sum / static_cast<double>(d);
        } else {
            out.mean_nbrs_dist[v] = R_PosInf;
        }
    }

    // Hop-k weighted-distance metric with BFS on hop count.
    // This matches the R implementation's shortest-hop discovery rule.
    std::vector<int> mark(n, 0);
    std::vector<int> hop_dist(n, -1);
    std::vector<double> wdist(n, 0.0);
    std::vector<int> queue;
    queue.reserve(1024);
    int cur_mark = 0;

    for (size_t src = 0; src < n; ++src) {
        if (cur_mark >= INT_MAX - 2) {
            std::fill(mark.begin(), mark.end(), 0);
            cur_mark = 0;
        }
        ++cur_mark;

        queue.clear();
        queue.push_back(static_cast<int>(src));
        mark[src] = cur_mark;
        hop_dist[src] = 0;
        wdist[src] = 0.0;

        size_t head = 0;
        while (head < queue.size()) {
            const int current = queue[head++];
            const int current_hop = hop_dist[static_cast<size_t>(current)];
            if (current_hop >= hop_k) {
                continue;
            }

            const auto& nbrs = adj[static_cast<size_t>(current)];
            const auto& w = edgelen[static_cast<size_t>(current)];
            const size_t nn = nbrs.size();
            for (size_t j = 0; j < nn; ++j) {
                const int nb = nbrs[j];
                if (nb < 0 || static_cast<size_t>(nb) >= n) {
                    Rf_error("Neighbor index out of range at vertex %d", current + 1);
                }
                const size_t nb_u = static_cast<size_t>(nb);
                if (mark[nb_u] != cur_mark) {
                    mark[nb_u] = cur_mark;
                    hop_dist[nb_u] = current_hop + 1;
                    wdist[nb_u] = wdist[static_cast<size_t>(current)] + w[j];
                    queue.push_back(nb);
                }
            }
        }

        double sum = 0.0;
        int count = 0;
        for (const int v : queue) {
            const size_t vu = static_cast<size_t>(v);
            if (hop_dist[vu] == hop_k) {
                sum += wdist[vu];
                ++count;
            }
        }
        out.mean_hopk_dist[src] = (count > 0) ? (sum / static_cast<double>(count)) : R_PosInf;
    }

    // Percentiles matching R semantics:
    // p.mean.* -> sum(metric <= x) / n
    // p.deg    -> sum(deg >= x) / n
    const double denom = static_cast<double>(n);
    std::vector<double> sorted_nbrs = out.mean_nbrs_dist;
    std::sort(sorted_nbrs.begin(), sorted_nbrs.end());

    std::vector<double> sorted_hopk = out.mean_hopk_dist;
    std::sort(sorted_hopk.begin(), sorted_hopk.end());

    std::vector<int> sorted_deg = out.deg;
    std::sort(sorted_deg.begin(), sorted_deg.end());

    for (size_t i = 0; i < n; ++i) {
        const double x_nbr = out.mean_nbrs_dist[i];
        const auto it_nbr = std::upper_bound(sorted_nbrs.begin(), sorted_nbrs.end(), x_nbr);
        out.p_mean_nbrs_dist[i] = static_cast<double>(std::distance(sorted_nbrs.begin(), it_nbr)) / denom;

        const double x_hop = out.mean_hopk_dist[i];
        const auto it_hop = std::upper_bound(sorted_hopk.begin(), sorted_hopk.end(), x_hop);
        out.p_mean_hopk_dist[i] = static_cast<double>(std::distance(sorted_hopk.begin(), it_hop)) / denom;

        const int x_deg = out.deg[i];
        const auto it_deg = std::lower_bound(sorted_deg.begin(), sorted_deg.end(), x_deg);
        out.p_deg[i] = static_cast<double>(std::distance(it_deg, sorted_deg.end())) / denom;
    }

    return out;
}

static SEXP vertex_metrics_to_R(const vertex_metrics_t& metrics) {
    const R_xlen_t n = static_cast<R_xlen_t>(metrics.mean_nbrs_dist.size());
    SEXP r_list = PROTECT(Rf_allocVector(VECSXP, 8));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 8));

    SEXP r_mean_nbrs = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_mean_hopk = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_deg = PROTECT(Rf_allocVector(INTSXP, n));
    SEXP r_p_mean_nbrs = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_p_mean_hopk = PROTECT(Rf_allocVector(REALSXP, n));
    SEXP r_p_deg = PROTECT(Rf_allocVector(REALSXP, n));

    std::copy(metrics.mean_nbrs_dist.begin(), metrics.mean_nbrs_dist.end(), REAL(r_mean_nbrs));
    std::copy(metrics.mean_hopk_dist.begin(), metrics.mean_hopk_dist.end(), REAL(r_mean_hopk));
    for (R_xlen_t i = 0; i < n; ++i) {
        INTEGER(r_deg)[i] = metrics.deg[static_cast<size_t>(i)];
    }
    std::copy(metrics.p_mean_nbrs_dist.begin(), metrics.p_mean_nbrs_dist.end(), REAL(r_p_mean_nbrs));
    std::copy(metrics.p_mean_hopk_dist.begin(), metrics.p_mean_hopk_dist.end(), REAL(r_p_mean_hopk));
    std::copy(metrics.p_deg.begin(), metrics.p_deg.end(), REAL(r_p_deg));

    SET_VECTOR_ELT(r_list, 0, r_mean_nbrs);
    SET_STRING_ELT(r_names, 0, Rf_mkChar("mean.nbrs.dist"));

    SET_VECTOR_ELT(r_list, 1, r_mean_hopk);
    SET_STRING_ELT(r_names, 1, Rf_mkChar("mean.hopk.dist"));

    SET_VECTOR_ELT(r_list, 2, r_deg);
    SET_STRING_ELT(r_names, 2, Rf_mkChar("deg"));

    SET_VECTOR_ELT(r_list, 3, r_p_mean_nbrs);
    SET_STRING_ELT(r_names, 3, Rf_mkChar("p.mean.nbrs.dist"));

    SET_VECTOR_ELT(r_list, 4, r_p_mean_hopk);
    SET_STRING_ELT(r_names, 4, Rf_mkChar("p.mean.hopk.dist"));

    SET_VECTOR_ELT(r_list, 5, r_p_deg);
    SET_STRING_ELT(r_names, 5, Rf_mkChar("p.deg"));

    SET_VECTOR_ELT(r_list, 6, Rf_ScalarInteger(static_cast<int>(n)));
    SET_STRING_ELT(r_names, 6, Rf_mkChar("n.vertices"));

    SET_VECTOR_ELT(r_list, 7, Rf_ScalarInteger(metrics.hop_k));
    SET_STRING_ELT(r_names, 7, Rf_mkChar("hop.k"));

    Rf_setAttrib(r_list, R_NamesSymbol, r_names);

    SEXP r_class = PROTECT(Rf_mkString("basin_vertex_metrics"));
    Rf_setAttrib(r_list, R_ClassSymbol, r_class);

    UNPROTECT(9);
    return r_list;
}

static vertex_metrics_t vertex_metrics_from_R(SEXP s_metrics, int hop_k_expected, size_t n_expected) {
    if (TYPEOF(s_metrics) != VECSXP) {
        Rf_error("vertex.metrics must be a list");
    }

    SEXP s_mean_nbrs = get_list_element(s_metrics, "mean.nbrs.dist");
    SEXP s_mean_hopk = get_list_element(s_metrics, "mean.hopk.dist");
    SEXP s_deg = get_list_element(s_metrics, "deg");
    SEXP s_p_mean_nbrs = get_list_element(s_metrics, "p.mean.nbrs.dist");
    SEXP s_p_mean_hopk = get_list_element(s_metrics, "p.mean.hopk.dist");
    SEXP s_p_deg = get_list_element(s_metrics, "p.deg");

    if (Rf_isNull(s_mean_nbrs) || Rf_isNull(s_mean_hopk) || Rf_isNull(s_deg) ||
        Rf_isNull(s_p_mean_nbrs) || Rf_isNull(s_p_mean_hopk) || Rf_isNull(s_p_deg)) {
        Rf_error("vertex.metrics is missing required fields");
    }

    const R_xlen_t n = XLENGTH(s_mean_nbrs);
    if (n <= 0 || static_cast<size_t>(n) != n_expected) {
        Rf_error("vertex.metrics length does not match number of vertices");
    }

    if (XLENGTH(s_mean_hopk) != n || XLENGTH(s_deg) != n ||
        XLENGTH(s_p_mean_nbrs) != n || XLENGTH(s_p_mean_hopk) != n ||
        XLENGTH(s_p_deg) != n) {
        Rf_error("vertex.metrics fields must have equal lengths");
    }

    SEXP s_hop_k = get_list_element(s_metrics, "hop.k");
    if (!Rf_isNull(s_hop_k) && LENGTH(s_hop_k) == 1) {
        const int hk = Rf_asInteger(s_hop_k);
        if (hk != hop_k_expected) {
            Rf_error("vertex.metrics$hop.k (%d) does not match requested hop.k (%d)", hk, hop_k_expected);
        }
    }

    vertex_metrics_t out;
    out.hop_k = hop_k_expected;
    out.mean_nbrs_dist.resize(static_cast<size_t>(n));
    out.mean_hopk_dist.resize(static_cast<size_t>(n));
    out.deg.resize(static_cast<size_t>(n));
    out.p_mean_nbrs_dist.resize(static_cast<size_t>(n));
    out.p_mean_hopk_dist.resize(static_cast<size_t>(n));
    out.p_deg.resize(static_cast<size_t>(n));

    for (R_xlen_t i = 0; i < n; ++i) {
        out.mean_nbrs_dist[static_cast<size_t>(i)] = REAL(s_mean_nbrs)[i];
        out.mean_hopk_dist[static_cast<size_t>(i)] = REAL(s_mean_hopk)[i];
        out.p_mean_nbrs_dist[static_cast<size_t>(i)] = REAL(s_p_mean_nbrs)[i];
        out.p_mean_hopk_dist[static_cast<size_t>(i)] = REAL(s_p_mean_hopk)[i];
        out.p_deg[static_cast<size_t>(i)] = REAL(s_p_deg)[i];
    }

    if (TYPEOF(s_deg) == INTSXP) {
        for (R_xlen_t i = 0; i < n; ++i) {
            out.deg[static_cast<size_t>(i)] = INTEGER(s_deg)[i];
        }
    } else if (TYPEOF(s_deg) == REALSXP) {
        for (R_xlen_t i = 0; i < n; ++i) {
            out.deg[static_cast<size_t>(i)] = static_cast<int>(REAL(s_deg)[i]);
        }
    } else {
        Rf_error("vertex.metrics$deg must be integer or numeric");
    }

    return out;
}

static std::vector<basin_row_t> extract_basin_rows(
    SEXP s_basins,
    const vertex_metrics_t& metrics,
    double mean_y
) {
    std::vector<basin_row_t> rows;
    const int n_basins = LENGTH(s_basins);
    rows.reserve(static_cast<size_t>(n_basins));
    const size_t n_vertices = metrics.mean_nbrs_dist.size();

    for (int i = 0; i < n_basins; ++i) {
        SEXP s_basin = VECTOR_ELT(s_basins, i);
        if (TYPEOF(s_basin) != VECSXP || LENGTH(s_basin) < 4) {
            Rf_error("Invalid basin structure");
        }

        const int vertex = Rf_asInteger(VECTOR_ELT(s_basin, 0));
        if (vertex == NA_INTEGER || vertex < 1 || static_cast<size_t>(vertex) > n_vertices) {
            Rf_error("Basin vertex index out of bounds");
        }
        const size_t idx0 = static_cast<size_t>(vertex - 1);

        const double value = Rf_asReal(VECTOR_ELT(s_basin, 1));
        const int hop_idx = Rf_asInteger(VECTOR_ELT(s_basin, 2));
        const int basin_size = get_nrows_from_dim(VECTOR_ELT(s_basin, 3));

        basin_row_t row;
        row.vertex = vertex;
        row.value = value;
        row.rel_value = value / mean_y;
        row.hop_idx = hop_idx;
        row.basin_size = basin_size;
        row.p_mean_nbrs_dist = metrics.p_mean_nbrs_dist[idx0];
        row.p_mean_hopk_dist = metrics.p_mean_hopk_dist[idx0];
        row.deg = static_cast<double>(metrics.deg[idx0]);
        row.p_deg = metrics.p_deg[idx0];
        rows.push_back(row);
    }

    return rows;
}

static SEXP basin_rows_to_df(
    std::vector<basin_row_t>& min_rows,
    std::vector<basin_row_t>& max_rows
) {
    std::stable_sort(min_rows.begin(), min_rows.end(),
                     [](const basin_row_t& a, const basin_row_t& b) {
                         return a.value < b.value;
                     });
    std::stable_sort(max_rows.begin(), max_rows.end(),
                     [](const basin_row_t& a, const basin_row_t& b) {
                         return a.value > b.value;
                     });

    const R_xlen_t n_min = static_cast<R_xlen_t>(min_rows.size());
    const R_xlen_t n_max = static_cast<R_xlen_t>(max_rows.size());
    const R_xlen_t n_total = n_min + n_max;

    SEXP r_label = PROTECT(Rf_allocVector(STRSXP, n_total));
    SEXP r_vertex = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP r_value = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_rel_value = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_type = PROTECT(Rf_allocVector(STRSXP, n_total));
    SEXP r_hop_idx = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_basin_size = PROTECT(Rf_allocVector(INTSXP, n_total));
    SEXP r_p_mean_nbrs = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_p_mean_hopk = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_deg = PROTECT(Rf_allocVector(REALSXP, n_total));
    SEXP r_p_deg = PROTECT(Rf_allocVector(REALSXP, n_total));

    R_xlen_t row = 0;
    for (R_xlen_t i = 0; i < n_min; ++i, ++row) {
        const basin_row_t& x = min_rows[static_cast<size_t>(i)];
        SET_STRING_ELT(r_label, row, Rf_mkChar(("m" + std::to_string(i + 1)).c_str()));
        INTEGER(r_vertex)[row] = x.vertex;
        REAL(r_value)[row] = x.value;
        REAL(r_rel_value)[row] = x.rel_value;
        SET_STRING_ELT(r_type, row, Rf_mkChar("min"));
        REAL(r_hop_idx)[row] = (x.hop_idx == NA_INTEGER) ? NA_REAL : static_cast<double>(x.hop_idx);
        INTEGER(r_basin_size)[row] = x.basin_size;
        REAL(r_p_mean_nbrs)[row] = x.p_mean_nbrs_dist;
        REAL(r_p_mean_hopk)[row] = x.p_mean_hopk_dist;
        REAL(r_deg)[row] = x.deg;
        REAL(r_p_deg)[row] = x.p_deg;
    }

    for (R_xlen_t i = 0; i < n_max; ++i, ++row) {
        const basin_row_t& x = max_rows[static_cast<size_t>(i)];
        SET_STRING_ELT(r_label, row, Rf_mkChar(("M" + std::to_string(i + 1)).c_str()));
        INTEGER(r_vertex)[row] = x.vertex;
        REAL(r_value)[row] = x.value;
        REAL(r_rel_value)[row] = x.rel_value;
        SET_STRING_ELT(r_type, row, Rf_mkChar("max"));
        REAL(r_hop_idx)[row] = (x.hop_idx == NA_INTEGER) ? NA_REAL : static_cast<double>(x.hop_idx);
        INTEGER(r_basin_size)[row] = x.basin_size;
        REAL(r_p_mean_nbrs)[row] = x.p_mean_nbrs_dist;
        REAL(r_p_mean_hopk)[row] = x.p_mean_hopk_dist;
        REAL(r_deg)[row] = x.deg;
        REAL(r_p_deg)[row] = x.p_deg;
    }

    SEXP r_df = PROTECT(Rf_allocVector(VECSXP, 11));
    SET_VECTOR_ELT(r_df, 0, r_label);
    SET_VECTOR_ELT(r_df, 1, r_vertex);
    SET_VECTOR_ELT(r_df, 2, r_value);
    SET_VECTOR_ELT(r_df, 3, r_rel_value);
    SET_VECTOR_ELT(r_df, 4, r_type);
    SET_VECTOR_ELT(r_df, 5, r_hop_idx);
    SET_VECTOR_ELT(r_df, 6, r_basin_size);
    SET_VECTOR_ELT(r_df, 7, r_p_mean_nbrs);
    SET_VECTOR_ELT(r_df, 8, r_p_mean_hopk);
    SET_VECTOR_ELT(r_df, 9, r_deg);
    SET_VECTOR_ELT(r_df, 10, r_p_deg);

    SEXP r_col_names = PROTECT(Rf_allocVector(STRSXP, 11));
    SET_STRING_ELT(r_col_names, 0, Rf_mkChar("label"));
    SET_STRING_ELT(r_col_names, 1, Rf_mkChar("vertex"));
    SET_STRING_ELT(r_col_names, 2, Rf_mkChar("value"));
    SET_STRING_ELT(r_col_names, 3, Rf_mkChar("rel.value"));
    SET_STRING_ELT(r_col_names, 4, Rf_mkChar("type"));
    SET_STRING_ELT(r_col_names, 5, Rf_mkChar("hop.idx"));
    SET_STRING_ELT(r_col_names, 6, Rf_mkChar("basin.size"));
    SET_STRING_ELT(r_col_names, 7, Rf_mkChar("p.mean.nbrs.dist"));
    SET_STRING_ELT(r_col_names, 8, Rf_mkChar("p.mean.hopk.dist"));
    SET_STRING_ELT(r_col_names, 9, Rf_mkChar("deg"));
    SET_STRING_ELT(r_col_names, 10, Rf_mkChar("p.deg"));
    Rf_setAttrib(r_df, R_NamesSymbol, r_col_names);

    SEXP r_rownames = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(r_rownames)[0] = NA_INTEGER;
    INTEGER(r_rownames)[1] = -static_cast<int>(n_total);
    Rf_setAttrib(r_df, R_RowNamesSymbol, r_rownames);

    SEXP r_class = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_class, 0, Rf_mkChar("basin_summary"));
    SET_STRING_ELT(r_class, 1, Rf_mkChar("data.frame"));
    Rf_setAttrib(r_df, R_ClassSymbol, r_class);

    UNPROTECT(15);
    return r_df;
}

} // namespace

extern "C" SEXP S_precompute_basin_vertex_metrics(
    SEXP s_adj_list,
    SEXP s_edgelen_list,
    SEXP s_hop_k
) {
    const int hop_k = Rf_asInteger(s_hop_k);
    if (hop_k == NA_INTEGER || hop_k < 1) {
        Rf_error("hop.k must be a positive integer");
    }

    auto adj = convert_adj_list_from_R(s_adj_list);
    auto edgelen = convert_weight_list_from_R(s_edgelen_list);
    validate_graph_inputs(adj, edgelen);

    vertex_metrics_t metrics = compute_vertex_metrics_cpp(adj, edgelen, hop_k);
    return vertex_metrics_to_R(metrics);
}

extern "C" SEXP S_summary_basins_of_attraction_cpp(
    SEXP s_object,
    SEXP s_adj_list,
    SEXP s_edgelen_list,
    SEXP s_hop_k,
    SEXP s_vertex_metrics
) {
    if (TYPEOF(s_object) != VECSXP) {
        Rf_error("object must be of class 'basins_of_attraction'");
    }

    const int hop_k = Rf_asInteger(s_hop_k);
    if (hop_k == NA_INTEGER || hop_k < 1) {
        Rf_error("hop.k must be a positive integer");
    }

    SEXP s_lmin = get_list_element(s_object, "lmin_basins");
    SEXP s_lmax = get_list_element(s_object, "lmax_basins");
    SEXP s_y = get_list_element(s_object, "y");
    if (Rf_isNull(s_lmin) || Rf_isNull(s_lmax) || Rf_isNull(s_y)) {
        Rf_error("object must include lmin_basins, lmax_basins, and y");
    }
    if (!Rf_isReal(s_y)) {
        Rf_error("object$y must be numeric");
    }

    const R_xlen_t n = XLENGTH(s_y);
    if (n <= 0) {
        Rf_error("object$y must be non-empty");
    }

    double sum_y = 0.0;
    for (R_xlen_t i = 0; i < n; ++i) {
        sum_y += REAL(s_y)[i];
    }
    const double mean_y = sum_y / static_cast<double>(n);

    vertex_metrics_t metrics;
    if (Rf_isNull(s_vertex_metrics)) {
        if (Rf_isNull(s_adj_list) || Rf_isNull(s_edgelen_list)) {
            Rf_error("adj.list and edgelen.list must be provided when vertex.metrics is NULL");
        }
        auto adj = convert_adj_list_from_R(s_adj_list);
        auto edgelen = convert_weight_list_from_R(s_edgelen_list);
        validate_graph_inputs(adj, edgelen);
        if (adj.size() != static_cast<size_t>(n)) {
            Rf_error("Length of graph inputs must match object$y length");
        }
        metrics = compute_vertex_metrics_cpp(adj, edgelen, hop_k);
    } else {
        metrics = vertex_metrics_from_R(s_vertex_metrics, hop_k, static_cast<size_t>(n));
    }

    std::vector<basin_row_t> min_rows = extract_basin_rows(s_lmin, metrics, mean_y);
    std::vector<basin_row_t> max_rows = extract_basin_rows(s_lmax, metrics, mean_y);
    return basin_rows_to_df(min_rows, max_rows);
}
