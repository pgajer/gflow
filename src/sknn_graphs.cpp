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
                                    SEXP s_connect_method) {
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
    if (connect_method < 0 || connect_method > 1) {
        Rf_error("connect_method must be 0 (component.mst) or 1 (global.mst).");
    }

    const double* X = REAL(s_X);
    std::vector<double> D = pairwise_distances(X, n, p);
    std::vector<std::vector<int>> nn = knn_sets_from_dist(D, n, k);

    std::map<std::uint64_t, edge_value_t> edges;
    for (int i = 0; i < n; ++i) {
        for (int j : nn[static_cast<size_t>(i)]) {
            add_edge(edges, i, j, D[static_cast<size_t>(i) * n + j]);
        }
    }

    int n_components_before = 0;
    std::vector<int> component_before = components_from_edges(edges, n, n_components_before);

    std::vector<mst_edge_t> added_edges;
    if (connect_components && n_components_before > 1) {
        if (connect_method == 0) {
            add_component_mst_edges(edges, D, n, component_before, n_components_before, added_edges);
        } else {
            add_global_mst_edges(edges, D, n, added_edges);
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

    SEXP result = PROTECT(Rf_allocVector(VECSXP, 18));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 18));
    const char* name_values[] = {
        "adj_list", "weight_list", "edge_matrix", "edge_weight",
        "knn_index", "n_vertices", "n_edges", "k",
        "connect_components", "connect_method", "edge_weight_type",
        "n_components_before", "n_components_after",
        "component_id_before", "component_id_after",
        "mst_edge_matrix", "mst_edge_weight", "n_mst_edges_added"
    };
    for (int i = 0; i < 18; ++i) {
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
    SET_VECTOR_ELT(result, 9, Rf_mkString(connect_method == 0 ? "component.mst" : "global.mst"));
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

    UNPROTECT(1);
    return result;
}
