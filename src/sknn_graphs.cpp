#include "sknn_graphs_r.h"

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

struct edge_value_t {
    double weight;
};

struct mst_edge_t {
    int u;
    int v;
    double weight;
};

struct dsu_t {
    std::vector<int> parent;
    std::vector<int> rank;

    explicit dsu_t(int n) : parent(static_cast<size_t>(n)), rank(static_cast<size_t>(n), 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[static_cast<size_t>(x)] != x) {
            parent[static_cast<size_t>(x)] = find(parent[static_cast<size_t>(x)]);
        }
        return parent[static_cast<size_t>(x)];
    }

    bool unite(int a, int b) {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) {
            return false;
        }
        if (rank[static_cast<size_t>(ra)] < rank[static_cast<size_t>(rb)]) {
            std::swap(ra, rb);
        }
        parent[static_cast<size_t>(rb)] = ra;
        if (rank[static_cast<size_t>(ra)] == rank[static_cast<size_t>(rb)]) {
            rank[static_cast<size_t>(ra)] += 1;
        }
        return true;
    }
};

std::uint64_t edge_key(int u, int v) {
    if (u > v) {
        std::swap(u, v);
    }
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(u)) << 32U) |
           static_cast<std::uint32_t>(v);
}

std::pair<int, int> unpack_edge_key(std::uint64_t key) {
    const int u = static_cast<int>(key >> 32U);
    const int v = static_cast<int>(key & 0xffffffffU);
    return {u, v};
}

double euclidean_distance(const double* X, int n, int p, int i, int j) {
    double acc = 0.0;
    for (int col = 0; col < p; ++col) {
        const double diff = X[i + n * col] - X[j + n * col];
        acc += diff * diff;
    }
    return std::sqrt(acc);
}

std::vector<double> pairwise_distances(const double* X, int n, int p) {
    std::vector<double> D(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d = euclidean_distance(X, n, p, i, j);
            D[static_cast<size_t>(i) * n + j] = d;
            D[static_cast<size_t>(j) * n + i] = d;
        }
    }
    return D;
}

std::vector<std::vector<int>> knn_sets_from_dist(const std::vector<double>& D, int n, int k) {
    std::vector<std::vector<int>> nn(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        std::vector<int> ord;
        ord.reserve(static_cast<size_t>(n - 1));
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                ord.push_back(j);
            }
        }
        std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            const double da = D[static_cast<size_t>(i) * n + a];
            const double db = D[static_cast<size_t>(i) * n + b];
            if (da < db) {
                return true;
            }
            if (da > db) {
                return false;
            }
            return a < b;
        });
        nn[static_cast<size_t>(i)].assign(ord.begin(), ord.begin() + k);
    }
    return nn;
}

std::vector<std::vector<int>> knn_sets_from_precomputed(SEXP s_knn_index, int n, int k) {
    if (!Rf_isMatrix(s_knn_index) || TYPEOF(s_knn_index) != INTSXP) {
        Rf_error("knn_index must be an integer matrix.");
    }
    SEXP s_dim = PROTECT(Rf_getAttrib(s_knn_index, R_DimSymbol));
    if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) < 2) {
        UNPROTECT(1);
        Rf_error("knn_index must have a valid dim attribute.");
    }
    const int nr = INTEGER(s_dim)[0];
    const int nc = INTEGER(s_dim)[1];
    UNPROTECT(1);
    if (nr != n || nc != k) {
        Rf_error("knn_index must have nrow(X) rows and k columns.");
    }

    const int* idx = INTEGER(s_knn_index);
    std::vector<std::vector<int>> nn(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        nn[static_cast<size_t>(i)].reserve(static_cast<size_t>(k));
        for (int col = 0; col < k; ++col) {
            const int j = idx[i + n * col] - 1;
            if (j < 0 || j >= n || j == i) {
                Rf_error("knn_index contains an invalid neighbor index.");
            }
            if (std::find(nn[static_cast<size_t>(i)].begin(),
                          nn[static_cast<size_t>(i)].end(),
                          j) != nn[static_cast<size_t>(i)].end()) {
                Rf_error("knn_index contains duplicate neighbors within a row.");
            }
            nn[static_cast<size_t>(i)].push_back(j);
        }
    }
    return nn;
}

void add_edge(std::map<std::uint64_t, edge_value_t>& edges, int u, int v, double weight) {
    if (u == v) {
        return;
    }
    const std::uint64_t key = edge_key(u, v);
    if (edges.find(key) == edges.end()) {
        edges.emplace(key, edge_value_t{weight});
    }
}

std::vector<std::vector<int>> adjacency_from_edges(const std::map<std::uint64_t, edge_value_t>& edges,
                                                   int n) {
    std::vector<std::vector<int>> adj(static_cast<size_t>(n));
    for (const auto& kv : edges) {
        const auto [u, v] = unpack_edge_key(kv.first);
        adj[static_cast<size_t>(u)].push_back(v);
        adj[static_cast<size_t>(v)].push_back(u);
    }
    for (auto& row : adj) {
        std::sort(row.begin(), row.end());
    }
    return adj;
}

std::vector<int> components_from_edges(const std::map<std::uint64_t, edge_value_t>& edges,
                                       int n,
                                       int& n_components) {
    std::vector<std::vector<int>> adj = adjacency_from_edges(edges, n);
    std::vector<int> comp(static_cast<size_t>(n), -1);
    n_components = 0;
    std::vector<int> stack;
    for (int s = 0; s < n; ++s) {
        if (comp[static_cast<size_t>(s)] >= 0) {
            continue;
        }
        comp[static_cast<size_t>(s)] = n_components;
        stack.assign(1, s);
        while (!stack.empty()) {
            const int u = stack.back();
            stack.pop_back();
            for (int v : adj[static_cast<size_t>(u)]) {
                if (comp[static_cast<size_t>(v)] < 0) {
                    comp[static_cast<size_t>(v)] = n_components;
                    stack.push_back(v);
                }
            }
        }
        n_components += 1;
    }
    return comp;
}

void add_component_mst_edges(std::map<std::uint64_t, edge_value_t>& edges,
                             const std::vector<double>& D,
                             int n,
                             const std::vector<int>& comp,
                             int n_components,
                             std::vector<mst_edge_t>& added_edges) {
    if (n_components <= 1) {
        return;
    }

    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> best(static_cast<size_t>(n_components) * n_components, inf);
    std::vector<int> best_u(static_cast<size_t>(n_components) * n_components, -1);
    std::vector<int> best_v(static_cast<size_t>(n_components) * n_components, -1);

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const int ci = comp[static_cast<size_t>(i)];
            const int cj = comp[static_cast<size_t>(j)];
            if (ci == cj) {
                continue;
            }
            const int a = std::min(ci, cj);
            const int b = std::max(ci, cj);
            const size_t idx = static_cast<size_t>(a) * n_components + b;
            const double d = D[static_cast<size_t>(i) * n + j];
            if (d < best[idx] || (d == best[idx] &&
                                  (i < best_u[idx] || (i == best_u[idx] && j < best_v[idx])))) {
                best[idx] = d;
                best_u[idx] = i;
                best_v[idx] = j;
            }
        }
    }

    std::vector<mst_edge_t> candidates;
    for (int a = 0; a < n_components - 1; ++a) {
        for (int b = a + 1; b < n_components; ++b) {
            const size_t idx = static_cast<size_t>(a) * n_components + b;
            if (std::isfinite(best[idx])) {
                candidates.push_back(mst_edge_t{best_u[idx], best_v[idx], best[idx]});
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const mst_edge_t& x, const mst_edge_t& y) {
        if (x.weight < y.weight) return true;
        if (x.weight > y.weight) return false;
        if (x.u < y.u) return true;
        if (x.u > y.u) return false;
        return x.v < y.v;
    });

    dsu_t dsu(n_components);
    for (const auto& e : candidates) {
        const int cu = comp[static_cast<size_t>(e.u)];
        const int cv = comp[static_cast<size_t>(e.v)];
        if (dsu.unite(cu, cv)) {
            add_edge(edges, e.u, e.v, e.weight);
            added_edges.push_back(e);
            if (static_cast<int>(added_edges.size()) == n_components - 1) {
                break;
            }
        }
    }
}

bool try_component_mst_from_bridge_knn(std::map<std::uint64_t, edge_value_t>& edges,
                                       const double* X,
                                       int n,
                                       int p,
                                       const std::vector<std::vector<int>>& bridge_nn,
                                       int bridge_k,
                                       const std::vector<int>& comp,
                                       int n_components,
                                       std::vector<mst_edge_t>& added_edges) {
    if (n_components <= 1) {
        return true;
    }

    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> best(static_cast<size_t>(n_components) * n_components, inf);
    std::vector<int> best_u(static_cast<size_t>(n_components) * n_components, -1);
    std::vector<int> best_v(static_cast<size_t>(n_components) * n_components, -1);

    for (int i = 0; i < n; ++i) {
        const int ci = comp[static_cast<size_t>(i)];
        const int row_k = std::min(bridge_k, static_cast<int>(bridge_nn[static_cast<size_t>(i)].size()));
        for (int pos = 0; pos < row_k; ++pos) {
            const int j = bridge_nn[static_cast<size_t>(i)][static_cast<size_t>(pos)];
            const int cj = comp[static_cast<size_t>(j)];
            if (ci == cj) {
                continue;
            }
            const int a = std::min(ci, cj);
            const int b = std::max(ci, cj);
            const size_t idx = static_cast<size_t>(a) * n_components + b;
            const double d = euclidean_distance(X, n, p, i, j);
            const int u = std::min(i, j);
            const int v = std::max(i, j);
            if (d < best[idx] || (d == best[idx] &&
                                  (u < best_u[idx] || (u == best_u[idx] && v < best_v[idx])))) {
                best[idx] = d;
                best_u[idx] = u;
                best_v[idx] = v;
            }
        }
    }

    std::vector<mst_edge_t> candidates;
    for (int a = 0; a < n_components - 1; ++a) {
        for (int b = a + 1; b < n_components; ++b) {
            const size_t idx = static_cast<size_t>(a) * n_components + b;
            if (std::isfinite(best[idx])) {
                candidates.push_back(mst_edge_t{best_u[idx], best_v[idx], best[idx]});
            }
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const mst_edge_t& x, const mst_edge_t& y) {
        if (x.weight < y.weight) return true;
        if (x.weight > y.weight) return false;
        if (x.u < y.u) return true;
        if (x.u > y.u) return false;
        return x.v < y.v;
    });

    dsu_t dsu(n_components);
    std::vector<mst_edge_t> proposed;
    proposed.reserve(static_cast<size_t>(n_components - 1));
    for (const auto& e : candidates) {
        const int cu = comp[static_cast<size_t>(e.u)];
        const int cv = comp[static_cast<size_t>(e.v)];
        if (dsu.unite(cu, cv)) {
            proposed.push_back(e);
            if (static_cast<int>(proposed.size()) == n_components - 1) {
                break;
            }
        }
    }
    if (static_cast<int>(proposed.size()) != n_components - 1) {
        return false;
    }

    for (const auto& e : proposed) {
        add_edge(edges, e.u, e.v, e.weight);
        added_edges.push_back(e);
    }
    return true;
}

void add_global_mst_edges(std::map<std::uint64_t, edge_value_t>& edges,
                          const std::vector<double>& D,
                          int n,
                          std::vector<mst_edge_t>& added_edges) {
    std::vector<mst_edge_t> candidates;
    candidates.reserve(static_cast<size_t>(n) * static_cast<size_t>(n - 1) / 2U);
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            candidates.push_back(mst_edge_t{i, j, D[static_cast<size_t>(i) * n + j]});
        }
    }
    std::sort(candidates.begin(), candidates.end(), [](const mst_edge_t& x, const mst_edge_t& y) {
        if (x.weight < y.weight) return true;
        if (x.weight > y.weight) return false;
        if (x.u < y.u) return true;
        if (x.u > y.u) return false;
        return x.v < y.v;
    });

    dsu_t dsu(n);
    int n_mst_edges = 0;
    for (const auto& e : candidates) {
        if (!dsu.unite(e.u, e.v)) {
            continue;
        }
        n_mst_edges += 1;
        const std::uint64_t key = edge_key(e.u, e.v);
        if (edges.find(key) == edges.end()) {
            add_edge(edges, e.u, e.v, e.weight);
            added_edges.push_back(e);
        }
        if (n_mst_edges == n - 1) {
            break;
        }
    }
}

SEXP make_integer_vector_1based(const std::vector<int>& values) {
    SEXP out = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(values.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(values.size()); ++i) {
        INTEGER(out)[i] = values[static_cast<size_t>(i)] + 1;
    }
    UNPROTECT(1);
    return out;
}

SEXP make_edge_matrix(const std::vector<mst_edge_t>& edges) {
    SEXP out = PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(edges.size()), 2));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(edges.size()); ++i) {
        INTEGER(out)[i] = edges[static_cast<size_t>(i)].u + 1;
        INTEGER(out)[i + static_cast<R_xlen_t>(edges.size())] = edges[static_cast<size_t>(i)].v + 1;
    }
    UNPROTECT(1);
    return out;
}

} // namespace

extern "C" SEXP S_create_sknn_graph(SEXP s_X,
                                    SEXP s_k,
                                    SEXP s_connect_components,
                                    SEXP s_connect_method,
                                    SEXP s_neighbor_method,
                                    SEXP s_ann_eps,
                                    SEXP s_knn_index,
                                    SEXP s_bridge_knn_index,
                                    SEXP s_bridge_k,
                                    SEXP s_bridge_k_max,
                                    SEXP s_bridge_growth) {
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

    const int k = Rf_asInteger(s_k);
    if (n < 2 || p < 1) {
        Rf_error("X must have at least two rows and one column.");
    }
    if (k < 1 || k >= n) {
        Rf_error("k must be a positive integer smaller than nrow(X).");
    }
    const bool connect_components = (Rf_asLogical(s_connect_components) == TRUE);
    const int connect_method = Rf_asInteger(s_connect_method);
    if (connect_method < 0 || connect_method > 2) {
        Rf_error("connect_method must be 0 (component.mst), 1 (global.mst), or 2 (component.mst.ann).");
    }
    const int neighbor_method = Rf_asInteger(s_neighbor_method);
    if (neighbor_method < 0 || neighbor_method > 1) {
        Rf_error("neighbor_method must be 0 (exact) or 1 (ann).");
    }
    const double ann_eps = Rf_asReal(s_ann_eps);
    if (!std::isfinite(ann_eps) || ann_eps < 0.0) {
        Rf_error("ann_eps must be a finite non-negative number.");
    }
    const int bridge_k = Rf_asInteger(s_bridge_k);
    const int bridge_k_max = Rf_asInteger(s_bridge_k_max);
    const double bridge_growth = Rf_asReal(s_bridge_growth);
    if (bridge_k < 1 || bridge_k >= n) {
        Rf_error("bridge_k must be a positive integer smaller than nrow(X).");
    }
    if (bridge_k_max < bridge_k || bridge_k_max >= n) {
        Rf_error("bridge_k_max must be an integer between bridge_k and nrow(X) - 1.");
    }
    if (!std::isfinite(bridge_growth) || bridge_growth <= 1.0) {
        Rf_error("bridge_growth must be a finite number greater than 1.");
    }

    const double* X = REAL(s_X);
    std::vector<double> D;
    std::vector<std::vector<int>> nn;
    if (neighbor_method == 0) {
        D = pairwise_distances(X, n, p);
        nn = knn_sets_from_dist(D, n, k);
    } else {
        nn = knn_sets_from_precomputed(s_knn_index, n, k);
    }

    auto get_distance = [&](int i, int j) {
        if (!D.empty()) {
            return D[static_cast<size_t>(i) * n + j];
        }
        return euclidean_distance(X, n, p, i, j);
    };

    std::map<std::uint64_t, edge_value_t> edges;
    for (int i = 0; i < n; ++i) {
        for (int j : nn[static_cast<size_t>(i)]) {
            add_edge(edges, i, j, get_distance(i, j));
        }
    }

    int n_components_before = 0;
    std::vector<int> component_before = components_from_edges(edges, n, n_components_before);

    std::vector<mst_edge_t> added_edges;
    bool bridge_exact_fallback_used = false;
    std::string bridge_method = "none";
    int bridge_k_used = NA_INTEGER;
    if (connect_components && n_components_before > 1) {
        if (connect_method == 0) {
            if (D.empty()) {
                D = pairwise_distances(X, n, p);
            }
            add_component_mst_edges(edges, D, n, component_before, n_components_before, added_edges);
            bridge_method = "exact";
        } else if (connect_method == 1) {
            if (D.empty()) {
                D = pairwise_distances(X, n, p);
            }
            add_global_mst_edges(edges, D, n, added_edges);
            bridge_method = "global.mst";
        } else {
            std::vector<std::vector<int>> bridge_nn =
                knn_sets_from_precomputed(s_bridge_knn_index, n, bridge_k_max);
            int current_k = bridge_k;
            while (true) {
                std::map<std::uint64_t, edge_value_t> candidate_edges = edges;
                std::vector<mst_edge_t> candidate_added;
                if (try_component_mst_from_bridge_knn(candidate_edges, X, n, p,
                                                      bridge_nn, current_k,
                                                      component_before,
                                                      n_components_before,
                                                      candidate_added)) {
                    edges.swap(candidate_edges);
                    added_edges.swap(candidate_added);
                    bridge_method = "ann";
                    bridge_k_used = current_k;
                    break;
                }
                if (current_k >= bridge_k_max) {
                    if (D.empty()) {
                        D = pairwise_distances(X, n, p);
                    }
                    add_component_mst_edges(edges, D, n, component_before,
                                            n_components_before, added_edges);
                    bridge_method = "ann_then_exact";
                    bridge_exact_fallback_used = true;
                    bridge_k_used = current_k;
                    break;
                }
                const int grown = static_cast<int>(std::ceil(static_cast<double>(current_k) * bridge_growth));
                current_k = std::min(bridge_k_max, std::max(current_k + 1, grown));
            }
        }
    }

    int n_components_after = 0;
    std::vector<int> component_after = components_from_edges(edges, n, n_components_after);
    std::vector<std::vector<std::pair<int, double>>> adj(static_cast<size_t>(n));
    std::vector<mst_edge_t> edge_list;
    edge_list.reserve(edges.size());
    for (const auto& kv : edges) {
        const auto [u, v] = unpack_edge_key(kv.first);
        const double w = kv.second.weight;
        adj[static_cast<size_t>(u)].push_back({v, w});
        adj[static_cast<size_t>(v)].push_back({u, w});
        edge_list.push_back(mst_edge_t{u, v, w});
    }
    for (auto& row : adj) {
        std::sort(row.begin(), row.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
    }

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 26));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 26));
    const char* name_values[] = {
        "adj_list", "weight_list", "edge_matrix", "edge_weight",
        "knn_index", "n_vertices", "n_edges", "k",
        "connect_components", "connect_method", "edge_weight_type",
        "n_components_before", "n_components_after",
        "component_id_before", "component_id_after",
        "mst_edge_matrix", "mst_edge_weight", "n_mst_edges_added",
        "neighbor_method", "ann_eps",
        "bridge_method", "bridge_k", "bridge_k_max", "bridge_growth",
        "bridge_k_used", "bridge_exact_fallback_used"
    };
    for (int i = 0; i < 26; ++i) {
        SET_STRING_ELT(names, i, Rf_mkChar(name_values[i]));
    }
    Rf_setAttrib(result, R_NamesSymbol, names);
    UNPROTECT(1);

    SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n));
    SEXP weight_list = PROTECT(Rf_allocVector(VECSXP, n));
    for (int i = 0; i < n; ++i) {
        const auto& row = adj[static_cast<size_t>(i)];
        SEXP a = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(row.size())));
        SEXP w = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(row.size())));
        for (R_xlen_t j = 0; j < static_cast<R_xlen_t>(row.size()); ++j) {
            INTEGER(a)[j] = row[static_cast<size_t>(j)].first + 1;
            REAL(w)[j] = row[static_cast<size_t>(j)].second;
        }
        SET_VECTOR_ELT(adj_list, i, a);
        SET_VECTOR_ELT(weight_list, i, w);
        UNPROTECT(2);
    }
    SET_VECTOR_ELT(result, 0, adj_list);
    SET_VECTOR_ELT(result, 1, weight_list);
    UNPROTECT(2);

    SEXP edge_matrix = PROTECT(make_edge_matrix(edge_list));
    SET_VECTOR_ELT(result, 2, edge_matrix);
    UNPROTECT(1);

    SEXP edge_weight = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(edge_list.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(edge_list.size()); ++i) {
        REAL(edge_weight)[i] = edge_list[static_cast<size_t>(i)].weight;
    }
    SET_VECTOR_ELT(result, 3, edge_weight);
    UNPROTECT(1);

    SEXP knn_index = PROTECT(Rf_allocMatrix(INTSXP, n, k));
    for (int col = 0; col < k; ++col) {
        for (int i = 0; i < n; ++i) {
            INTEGER(knn_index)[i + n * col] = nn[static_cast<size_t>(i)][static_cast<size_t>(col)] + 1;
        }
    }
    SET_VECTOR_ELT(result, 4, knn_index);
    UNPROTECT(1);

    SET_VECTOR_ELT(result, 5, Rf_ScalarInteger(n));
    SET_VECTOR_ELT(result, 6, Rf_ScalarInteger(static_cast<int>(edges.size())));
    SET_VECTOR_ELT(result, 7, Rf_ScalarInteger(k));
    SET_VECTOR_ELT(result, 8, Rf_ScalarLogical(connect_components ? TRUE : FALSE));
    SET_VECTOR_ELT(result, 9, Rf_mkString(connect_method == 0 ? "component.mst" :
                                          (connect_method == 1 ? "global.mst" : "component.mst.ann")));
    SET_VECTOR_ELT(result, 10, Rf_mkString("distance"));
    SET_VECTOR_ELT(result, 11, Rf_ScalarInteger(n_components_before));
    SET_VECTOR_ELT(result, 12, Rf_ScalarInteger(n_components_after));

    SEXP comp_before = PROTECT(make_integer_vector_1based(component_before));
    SEXP comp_after = PROTECT(make_integer_vector_1based(component_after));
    SET_VECTOR_ELT(result, 13, comp_before);
    SET_VECTOR_ELT(result, 14, comp_after);
    UNPROTECT(2);

    SEXP mst_edge_matrix = PROTECT(make_edge_matrix(added_edges));
    SET_VECTOR_ELT(result, 15, mst_edge_matrix);
    UNPROTECT(1);

    SEXP mst_edge_weight = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(added_edges.size())));
    for (R_xlen_t i = 0; i < static_cast<R_xlen_t>(added_edges.size()); ++i) {
        REAL(mst_edge_weight)[i] = added_edges[static_cast<size_t>(i)].weight;
    }
    SET_VECTOR_ELT(result, 16, mst_edge_weight);
    UNPROTECT(1);

    SET_VECTOR_ELT(result, 17, Rf_ScalarInteger(static_cast<int>(added_edges.size())));
    SET_VECTOR_ELT(result, 18, Rf_mkString(neighbor_method == 0 ? "exact" : "ann"));
    SET_VECTOR_ELT(result, 19, Rf_ScalarReal(ann_eps));
    SET_VECTOR_ELT(result, 20, Rf_mkString(bridge_method.c_str()));
    SET_VECTOR_ELT(result, 21, Rf_ScalarInteger(bridge_k));
    SET_VECTOR_ELT(result, 22, Rf_ScalarInteger(bridge_k_max));
    SET_VECTOR_ELT(result, 23, Rf_ScalarReal(bridge_growth));
    SET_VECTOR_ELT(result, 24, Rf_ScalarInteger(bridge_k_used));
    SET_VECTOR_ELT(result, 25, Rf_ScalarLogical(bridge_exact_fallback_used ? TRUE : FALSE));

    UNPROTECT(1);
    return result;
}
