#include "local_geodesic_pruning_r.h"

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace {

struct edge_t {
    int u;
    int v;
    double weight;
};

struct pruned_edge_t {
    int u;
    int v;
    double edge_weight;
    double alt_path_length;
    double path_edge_ratio;
};

std::uint64_t edge_key(int u, int v) {
    if (u > v) {
        std::swap(u, v);
    }
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32U) |
           static_cast<std::uint32_t>(v);
}

bool is_integerish(double x) {
    return std::isfinite(x) && x == std::floor(x);
}

int read_vertex_index(SEXP row, R_xlen_t pos, const char* name) {
    if (TYPEOF(row) == INTSXP) {
        const int value = INTEGER(row)[pos];
        if (value == NA_INTEGER) {
            Rf_error("%s contains NA vertex IDs.", name);
        }
        return value;
    }
    if (TYPEOF(row) == REALSXP) {
        const double value = REAL(row)[pos];
        if (!is_integerish(value)) {
            Rf_error("%s contains non-integer vertex IDs.", name);
        }
        if (value < static_cast<double>(std::numeric_limits<int>::min()) ||
            value > static_cast<double>(std::numeric_limits<int>::max())) {
            Rf_error("%s contains vertex IDs outside the integer range.", name);
        }
        return static_cast<int>(value);
    }
    Rf_error("%s must contain integer or numeric vectors.", name);
}

double read_weight(SEXP row, R_xlen_t pos) {
    if (TYPEOF(row) == REALSXP) {
        return REAL(row)[pos];
    }
    if (TYPEOF(row) == INTSXP) {
        const int value = INTEGER(row)[pos];
        if (value == NA_INTEGER) {
            return NA_REAL;
        }
        return static_cast<double>(value);
    }
    Rf_error("weight.list must contain numeric vectors.");
}

SEXP make_pruned_edge_stats(const std::vector<pruned_edge_t>& edges) {
    const int nrows = static_cast<int>(edges.size());
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 5));
    const char* name_values[] = {
        "u", "v", "edge_length", "alt_path_length", "path_edge_ratio"
    };
    for (int i = 0; i < 5; ++i) {
        SET_STRING_ELT(names, i, Rf_mkChar(name_values[i]));
    }
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(1);

    SEXP u = PROTECT(Rf_allocVector(INTSXP, nrows));
    SEXP v = PROTECT(Rf_allocVector(INTSXP, nrows));
    SEXP edge_length = PROTECT(Rf_allocVector(REALSXP, nrows));
    SEXP alt_path_length = PROTECT(Rf_allocVector(REALSXP, nrows));
    SEXP path_edge_ratio = PROTECT(Rf_allocVector(REALSXP, nrows));
    for (int i = 0; i < nrows; ++i) {
        INTEGER(u)[i] = edges[static_cast<size_t>(i)].u + 1;
        INTEGER(v)[i] = edges[static_cast<size_t>(i)].v + 1;
        REAL(edge_length)[i] = edges[static_cast<size_t>(i)].edge_weight;
        REAL(alt_path_length)[i] = edges[static_cast<size_t>(i)].alt_path_length;
        REAL(path_edge_ratio)[i] = edges[static_cast<size_t>(i)].path_edge_ratio;
    }
    SET_VECTOR_ELT(out, 0, u);
    SET_VECTOR_ELT(out, 1, v);
    SET_VECTOR_ELT(out, 2, edge_length);
    SET_VECTOR_ELT(out, 3, alt_path_length);
    SET_VECTOR_ELT(out, 4, path_edge_ratio);
    UNPROTECT(5);

    SEXP row_names = PROTECT(Rf_allocVector(INTSXP, nrows));
    for (int i = 0; i < nrows; ++i) {
        INTEGER(row_names)[i] = i + 1;
    }
    Rf_setAttrib(out, R_RowNamesSymbol, row_names);
    UNPROTECT(1);

    SEXP klass = PROTECT(Rf_mkString("data.frame"));
    Rf_setAttrib(out, R_ClassSymbol, klass);
    UNPROTECT(1);

    UNPROTECT(1);
    return out;
}

std::vector<std::vector<int>> exact_knn_index(const double* X, int n, int p, int k) {
    std::vector<std::vector<int>> out(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double, int>> distances;
        distances.reserve(static_cast<size_t>(n - 1));
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            double d2 = 0.0;
            for (int col = 0; col < p; ++col) {
                const double diff = X[i + n * col] - X[j + n * col];
                d2 += diff * diff;
            }
            distances.push_back({d2, j});
        }
        std::sort(distances.begin(), distances.end(), [](const auto& a, const auto& b) {
            if (a.first < b.first) return true;
            if (a.first > b.first) return false;
            return a.second < b.second;
        });
        out[static_cast<size_t>(i)].reserve(static_cast<size_t>(k));
        for (int r = 0; r < k; ++r) {
            out[static_cast<size_t>(i)].push_back(distances[static_cast<size_t>(r)].second);
        }
    }
    return out;
}

double local_dijkstra_distance(const std::vector<std::vector<std::pair<int, double>>>& adj,
                               const std::vector<int>& local_vertices,
                               int source,
                               int target,
                               std::uint64_t excluded_key,
                               double cutoff) {
    const int n = static_cast<int>(adj.size());
    std::vector<char> in_local(static_cast<size_t>(n), 0);
    for (int v : local_vertices) {
        in_local[static_cast<size_t>(v)] = 1;
    }
    in_local[static_cast<size_t>(source)] = 1;
    in_local[static_cast<size_t>(target)] = 1;

    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(static_cast<size_t>(n), inf);
    using queue_item_t = std::pair<double, int>;
    std::priority_queue<queue_item_t, std::vector<queue_item_t>, std::greater<queue_item_t>> pq;

    dist[static_cast<size_t>(source)] = 0.0;
    pq.push({0.0, source});

    while (!pq.empty()) {
        const double du = pq.top().first;
        const int u = pq.top().second;
        pq.pop();
        if (du != dist[static_cast<size_t>(u)]) {
            continue;
        }
        if (du > cutoff) {
            break;
        }
        if (u == target) {
            return du;
        }
        for (const auto& edge : adj[static_cast<size_t>(u)]) {
            const int v = edge.first;
            if (!in_local[static_cast<size_t>(v)]) {
                continue;
            }
            if (edge_key(u, v) == excluded_key) {
                continue;
            }
            const double nd = du + edge.second;
            if (nd < dist[static_cast<size_t>(v)] && nd <= cutoff) {
                dist[static_cast<size_t>(v)] = nd;
                pq.push({nd, v});
            }
        }
    }
    return inf;
}

void remove_undirected_edge(std::vector<std::vector<std::pair<int, double>>>& adj,
                            int u,
                            int v) {
    auto remove_one = [](std::vector<std::pair<int, double>>& row, int target) {
        const auto it = std::find_if(row.begin(), row.end(), [&](const auto& edge) {
            return edge.first == target;
        });
        if (it != row.end()) {
            row.erase(it);
        }
    };
    remove_one(adj[static_cast<size_t>(u)], v);
    remove_one(adj[static_cast<size_t>(v)], u);
}

std::vector<edge_t> validate_and_convert_graph(SEXP s_adj_list,
                                               SEXP s_weight_list,
                                               int n,
                                               std::vector<std::vector<std::pair<int, double>>>& adj) {
    if (!Rf_isVectorList(s_adj_list) || !Rf_isVectorList(s_weight_list)) {
        Rf_error("adj.list and weight.list must be lists.");
    }
    if (Rf_length(s_adj_list) != n || Rf_length(s_weight_list) != n) {
        Rf_error("adj.list and weight.list must have length nrow(X).");
    }

    adj.assign(static_cast<size_t>(n), {});
    std::vector<edge_t> edges;

    for (int i = 0; i < n; ++i) {
        SEXP nbrs = VECTOR_ELT(s_adj_list, i);
        SEXP weights = VECTOR_ELT(s_weight_list, i);
        if ((TYPEOF(nbrs) != INTSXP && TYPEOF(nbrs) != REALSXP) ||
            (TYPEOF(weights) != REALSXP && TYPEOF(weights) != INTSXP)) {
            Rf_error("adj.list and weight.list entries must be numeric vectors.");
        }
        const R_xlen_t row_len = Rf_xlength(nbrs);
        if (Rf_xlength(weights) != row_len) {
            Rf_error("adj.list and weight.list entries must have matching lengths.");
        }
        std::set<int> seen_neighbors;
        adj[static_cast<size_t>(i)].reserve(static_cast<size_t>(row_len));
        for (R_xlen_t pos = 0; pos < row_len; ++pos) {
            const int j1 = read_vertex_index(nbrs, pos, "adj.list");
            if (j1 < 1 || j1 > n) {
                Rf_error("adj.list contains vertex IDs outside 1..nrow(X).");
            }
            const int j = j1 - 1;
            if (j == i) {
                Rf_error("adj.list cannot contain self edges.");
            }
            if (!seen_neighbors.insert(j).second) {
                Rf_error("adj.list contains duplicate neighbors within a row.");
            }
            const double w = read_weight(weights, pos);
            if (!std::isfinite(w) || w <= 0.0) {
                Rf_error("weight.list must contain finite positive edge weights.");
            }
            adj[static_cast<size_t>(i)].push_back({j, w});
            if (i < j) {
                edges.push_back(edge_t{i, j, w});
            }
        }
    }

    const double weight_tol = 1e-12;
    for (int i = 0; i < n; ++i) {
        for (const auto& edge : adj[static_cast<size_t>(i)]) {
            const int j = edge.first;
            const double w = edge.second;
            const auto& reciprocal_row = adj[static_cast<size_t>(j)];
            const auto it = std::find_if(reciprocal_row.begin(), reciprocal_row.end(),
                                         [&](const auto& reciprocal) {
                                             return reciprocal.first == i;
                                         });
            if (it == reciprocal_row.end()) {
                Rf_error("adj.list must describe an undirected graph with reciprocal edges.");
            }
            const double scale = std::max({1.0, std::fabs(w), std::fabs(it->second)});
            if (std::fabs(w - it->second) > weight_tol * scale) {
                Rf_error("weight.list has mismatched reciprocal edge weights.");
            }
        }
    }

    return edges;
}

SEXP make_result(const std::vector<std::vector<std::pair<int, double>>>& adj,
                 int n_edges_before,
                 const std::vector<pruned_edge_t>& pruned_edges,
                 double prune_tau,
                 int prune_local_k,
                 bool with_pruned_edge_stats) {
    const int n = static_cast<int>(adj.size());
    int directed_edges_after = 0;
    for (const auto& row : adj) {
        directed_edges_after += static_cast<int>(row.size());
    }
    const int n_edges_after = directed_edges_after / 2;

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 9));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 9));
    const char* name_values[] = {
        "adj_list", "weight_list",
        "n_edges_before_pruning", "n_edges_after_pruning", "n_pruned_edges",
        "pruned_edge_stats", "prune_tau", "prune_local_k",
        "with_pruned_edge_stats"
    };
    for (int i = 0; i < 9; ++i) {
        SET_STRING_ELT(names, i, Rf_mkChar(name_values[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(1);

    SEXP out_adj = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP out_weight = PROTECT(Rf_allocVector(VECSXP, n));
    for (int i = 0; i < n; ++i) {
        const auto& row = adj[static_cast<size_t>(i)];
        SEXP a = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(row.size())));
        SEXP w = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(row.size())));
        for (R_xlen_t pos = 0; pos < static_cast<R_xlen_t>(row.size()); ++pos) {
            INTEGER(a)[pos] = row[static_cast<size_t>(pos)].first + 1;
            REAL(w)[pos] = row[static_cast<size_t>(pos)].second;
        }
        SET_VECTOR_ELT(out_adj, i, a);
        SET_VECTOR_ELT(out_weight, i, w);
        UNPROTECT(2);
    }
    SET_VECTOR_ELT(result, 0, out_adj);
    SET_VECTOR_ELT(result, 1, out_weight);
    UNPROTECT(2);

    SET_VECTOR_ELT(result, 2, Rf_ScalarInteger(n_edges_before));
    SET_VECTOR_ELT(result, 3, Rf_ScalarInteger(n_edges_after));
    SET_VECTOR_ELT(result, 4, Rf_ScalarInteger(n_edges_before - n_edges_after));

    SEXP stats = PROTECT(make_pruned_edge_stats(pruned_edges));
    SET_VECTOR_ELT(result, 5, stats);
    UNPROTECT(1);

    SET_VECTOR_ELT(result, 6, Rf_ScalarReal(prune_tau));
    SET_VECTOR_ELT(result, 7, Rf_ScalarInteger(prune_local_k));
    SET_VECTOR_ELT(result, 8, Rf_ScalarLogical(with_pruned_edge_stats ? TRUE : FALSE));

    UNPROTECT(1);
    return result;
}

} // namespace

extern "C" SEXP S_prune_graph_local_geodesic(SEXP s_X,
                                             SEXP s_adj_list,
                                             SEXP s_weight_list,
                                             SEXP s_prune_tau,
                                             SEXP s_prune_local_k,
                                             SEXP s_with_pruned_edge_stats) {
    if (!Rf_isMatrix(s_X) || TYPEOF(s_X) != REALSXP) {
        Rf_error("X must be a numeric matrix.");
    }
    SEXP s_dim = PROTECT(Rf_getAttrib(s_X, R_DimSymbol));
    if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) < 2) {
        UNPROTECT(1);
        Rf_error("X must have a valid dim attribute.");
    }
    const int n = INTEGER(s_dim)[0];
    const int p = INTEGER(s_dim)[1];
    UNPROTECT(1);
    if (n < 2 || p < 1) {
        Rf_error("X must have at least two rows and one column.");
    }

    const double prune_tau = Rf_asReal(s_prune_tau);
    if (!std::isfinite(prune_tau) || prune_tau <= 1.0) {
        Rf_error("prune.tau must be a finite number greater than 1.");
    }
    const int prune_local_k = Rf_asInteger(s_prune_local_k);
    if (prune_local_k < 1 || prune_local_k >= n) {
        Rf_error("prune.local.k must be a positive integer smaller than nrow(X).");
    }
    const int stats_flag = Rf_asLogical(s_with_pruned_edge_stats);
    if (stats_flag == NA_LOGICAL) {
        Rf_error("with.pruned.edge.stats must be TRUE or FALSE.");
    }
    const bool with_pruned_edge_stats = (stats_flag == TRUE);

    std::vector<std::vector<std::pair<int, double>>> adj;
    std::vector<edge_t> candidates = validate_and_convert_graph(s_adj_list, s_weight_list, n, adj);
    const int n_edges_before = static_cast<int>(candidates.size());
    if (n_edges_before == 0) {
        std::vector<pruned_edge_t> empty_stats;
        return make_result(adj, n_edges_before, empty_stats, prune_tau, prune_local_k,
                           with_pruned_edge_stats);
    }

    std::sort(candidates.begin(), candidates.end(), [](const edge_t& a, const edge_t& b) {
        if (a.weight > b.weight) return true;
        if (a.weight < b.weight) return false;
        if (a.u < b.u) return true;
        if (a.u > b.u) return false;
        return a.v < b.v;
    });

    const std::vector<std::vector<int>> local_nn =
        exact_knn_index(REAL(s_X), n, p, prune_local_k);

    std::set<std::uint64_t> alive_edges;
    for (const auto& edge : candidates) {
        alive_edges.insert(edge_key(edge.u, edge.v));
    }

    std::vector<pruned_edge_t> pruned_edges;
    if (with_pruned_edge_stats) {
        pruned_edges.reserve(candidates.size());
    }

    const double tol = 1e-12;
    for (const auto& candidate : candidates) {
        const std::uint64_t key = edge_key(candidate.u, candidate.v);
        if (alive_edges.find(key) == alive_edges.end()) {
            continue;
        }

        double edge_length = std::numeric_limits<double>::quiet_NaN();
        for (const auto& edge : adj[static_cast<size_t>(candidate.u)]) {
            if (edge.first == candidate.v) {
                edge_length = edge.second;
                break;
            }
        }
        if (!std::isfinite(edge_length)) {
            continue;
        }

        std::vector<int> local_vertices;
        local_vertices.reserve(static_cast<size_t>(2 + 2 * prune_local_k));
        local_vertices.push_back(candidate.u);
        local_vertices.push_back(candidate.v);
        for (int x : local_nn[static_cast<size_t>(candidate.u)]) {
            local_vertices.push_back(x);
        }
        for (int x : local_nn[static_cast<size_t>(candidate.v)]) {
            local_vertices.push_back(x);
        }
        std::sort(local_vertices.begin(), local_vertices.end());
        local_vertices.erase(std::unique(local_vertices.begin(), local_vertices.end()),
                             local_vertices.end());

        const double cutoff = prune_tau * edge_length;
        const double alt = local_dijkstra_distance(adj, local_vertices,
                                                   candidate.u, candidate.v,
                                                   key, cutoff + tol);
        if (std::isfinite(alt) && alt <= cutoff + tol) {
            alive_edges.erase(key);
            remove_undirected_edge(adj, candidate.u, candidate.v);
            if (with_pruned_edge_stats) {
                pruned_edges.push_back(pruned_edge_t{
                    candidate.u,
                    candidate.v,
                    edge_length,
                    alt,
                    alt / edge_length
                });
            }
        }
    }

    return make_result(adj, n_edges_before, pruned_edges, prune_tau, prune_local_k,
                       with_pruned_edge_stats);
}
