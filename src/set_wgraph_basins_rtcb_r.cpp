#include "set_wgraph.hpp"
#include "SEXP_cpp_conversion_utils.hpp"

#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <cstring>

#include <R.h>
#include <Rinternals.h>

namespace {

struct rtcb_params_t {
    size_t n_min = 50;
    double m_min = 0.7;
    double q_min = -1.0;          // if in (0,1], takes precedence over m_min
    int run_max = 2;
    double tau0 = 0.1;
    double kappa = 1.5;
    int k_max = 24;
    int h_max = 1000;
    double d_path_max = std::numeric_limits<double>::infinity();
    double eta_step = 0.0;
    double epsilon_d = 1e-12;
    double eps_num = 1e-12;
    bool sink_prune = false;
    int sink_prune_max_iter = 5;
    int max_paths_per_sink = 64;
};

struct arc_t {
    size_t to = 0;
    double w = 0.0;
};

struct rtcb_label_t {
    size_t vertex = 0;
    double d = 0.0;
    double G = 0.0;
    double L = 0.0;
    double base = 0.0; // lambda * G - L
    int run = 0;
    int hops = 0;
    int pred = -1;
    bool active = true;
};

struct rtcb_search_result_t {
    std::vector<rtcb_label_t> labels;
    std::vector<std::vector<int>> labels_at_vertex;
    std::unordered_map<size_t, std::vector<int>> accepted_sink_labels;
};

SEXP get_list_element(SEXP list, const char* name) {
    if (Rf_isNull(list) || !Rf_isNewList(list)) {
        return R_NilValue;
    }
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    if (Rf_isNull(names)) {
        return R_NilValue;
    }
    R_xlen_t n = Rf_xlength(list);
    for (R_xlen_t i = 0; i < n; ++i) {
        if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list, i);
        }
    }
    return R_NilValue;
}

double get_list_double(SEXP list, const char* name, double def) {
    SEXP elt = get_list_element(list, name);
    if (Rf_isNull(elt) || Rf_length(elt) == 0) {
        return def;
    }
    return Rf_asReal(elt);
}

int get_list_int(SEXP list, const char* name, int def) {
    SEXP elt = get_list_element(list, name);
    if (Rf_isNull(elt) || Rf_length(elt) == 0) {
        return def;
    }
    return Rf_asInteger(elt);
}

bool get_list_bool(SEXP list, const char* name, bool def) {
    SEXP elt = get_list_element(list, name);
    if (Rf_isNull(elt) || Rf_length(elt) == 0) {
        return def;
    }
    return Rf_asLogical(elt) != 0;
}

rtcb_params_t parse_rtcb_params(SEXP s_params, size_t n_vertices) {
    rtcb_params_t p;

    // Scale default n_min with graph size while keeping it conservative.
    size_t n_min_default = static_cast<size_t>(std::max(20.0, std::sqrt(static_cast<double>(std::max<size_t>(1, n_vertices)))));
    p.n_min = std::min(n_vertices, n_min_default);

    if (Rf_isNull(s_params) || !Rf_isNewList(s_params)) {
        return p;
    }

    p.n_min = static_cast<size_t>(std::max(1, get_list_int(s_params, "n_min", static_cast<int>(p.n_min))));
    if (n_vertices > 0) {
        p.n_min = std::min(p.n_min, n_vertices);
    }

    p.m_min = get_list_double(s_params, "m_min", p.m_min);
    p.q_min = get_list_double(s_params, "q_min", p.q_min);
    p.run_max = std::max(0, get_list_int(s_params, "run_max", p.run_max));
    p.tau0 = std::max(0.0, get_list_double(s_params, "tau0", p.tau0));
    p.kappa = std::max(1e-9, get_list_double(s_params, "kappa", p.kappa));
    p.k_max = std::max(1, get_list_int(s_params, "k_max", p.k_max));
    p.h_max = std::max(1, get_list_int(s_params, "h_max", p.h_max));
    p.d_path_max = get_list_double(s_params, "d_path_max", p.d_path_max);
    p.eta_step = std::max(0.0, get_list_double(s_params, "eta_step", p.eta_step));
    p.epsilon_d = std::max(1e-18, get_list_double(s_params, "epsilon_d", p.epsilon_d));
    p.eps_num = std::max(1e-18, get_list_double(s_params, "eps_num", p.eps_num));
    p.sink_prune = get_list_bool(s_params, "sink_prune", p.sink_prune);
    p.sink_prune_max_iter = std::max(1, get_list_int(s_params, "sink_prune_max_iter", p.sink_prune_max_iter));
    p.max_paths_per_sink = std::max(1, get_list_int(s_params, "max_paths_per_sink", p.max_paths_per_sink));

    return p;
}

double compute_lambda(const rtcb_params_t& params) {
    if (std::isfinite(params.q_min) && params.q_min > 0.0 && params.q_min <= 1.0) {
        return (1.0 - params.q_min) / std::max(1e-12, params.q_min);
    }
    double m = std::max(-0.999999, std::min(0.999999, params.m_min));
    return (1.0 - m) / std::max(1e-12, 1.0 + m);
}

bool dominates_label(
    const rtcb_label_t& a,
    const rtcb_label_t& b,
    double eps
    ) {
    bool no_worse = (a.d <= b.d + eps) &&
        (a.base >= b.base - eps) &&
        (a.run <= b.run) &&
        (a.hops <= b.hops);

    bool strictly_better = (a.d < b.d - eps) ||
        (a.base > b.base + eps) ||
        (a.run < b.run) ||
        (a.hops < b.hops);

    return no_worse && strictly_better;
}

void deduplicate_predecessors(std::unordered_map<size_t, std::vector<size_t>>& all_predecessors) {
    for (auto& kv : all_predecessors) {
        auto& preds = kv.second;
        std::sort(preds.begin(), preds.end());
        preds.erase(std::unique(preds.begin(), preds.end()), preds.end());
    }
}

SEXP basin_to_r(const gradient_basin_t& basin, bool with_trajectories, size_t n_vertices) {
    bool is_empty_basin = (basin.hop_idx == std::numeric_limits<size_t>::max());
    int n_fields = with_trajectories ? 8 : 7;

    SEXP r_basin = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP r_basin_names = PROTECT(Rf_allocVector(STRSXP, n_fields));

    int idx = 0;
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("vertex"));
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("value"));
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("hop_idx"));
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("basin_df"));
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("basin_bd_df"));
    SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("terminal_extrema"));
    if (with_trajectories) {
        SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("all_predecessors"));
        SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("trajectory_sets"));
    } else {
        SET_STRING_ELT(r_basin_names, idx++, Rf_mkChar("all_predecessors"));
    }
    Rf_setAttrib(r_basin, R_NamesSymbol, r_basin_names);

    idx = 0;
    SET_VECTOR_ELT(r_basin, idx++, Rf_ScalarInteger(static_cast<int>(basin.vertex + 1)));
    SET_VECTOR_ELT(r_basin, idx++, Rf_ScalarReal(basin.value));

    if (is_empty_basin) {
        SET_VECTOR_ELT(r_basin, idx++, Rf_ScalarInteger(NA_INTEGER));
    } else {
        SET_VECTOR_ELT(r_basin, idx++, Rf_ScalarInteger(static_cast<int>(basin.hop_idx)));
    }

    // basin_df
    {
        SEXP r_basin_df;
        if (is_empty_basin) {
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.hop_dist_map.size();
            r_basin_df = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(m), 2));
            double* pr = REAL(r_basin_df);
            size_t i = 0;
            for (const auto& kv : basin.hop_dist_map) {
                pr[i] = static_cast<double>(kv.first + 1);
                pr[i + m] = static_cast<double>(kv.second);
                i++;
            }
        }
        SEXP cn = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(cn, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(cn, 1, Rf_mkChar("hop_distance"));
        SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dn, 0, R_NilValue);
        SET_VECTOR_ELT(dn, 1, cn);
        Rf_setAttrib(r_basin_df, R_DimNamesSymbol, dn);
        SET_VECTOR_ELT(r_basin, idx++, r_basin_df);
        UNPROTECT(3);
    }

    // basin_bd_df
    {
        SEXP r_basin_bd_df;
        if (is_empty_basin) {
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, 0, 2));
        } else {
            size_t m = basin.y_nbhd_bd_map.size();
            r_basin_bd_df = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(m), 2));
            double* pr = REAL(r_basin_bd_df);
            size_t i = 0;
            for (const auto& kv : basin.y_nbhd_bd_map) {
                pr[i] = static_cast<double>(kv.first + 1);
                pr[i + m] = kv.second;
                i++;
            }
        }
        SEXP cn = PROTECT(Rf_allocVector(STRSXP, 2));
        SET_STRING_ELT(cn, 0, Rf_mkChar("vertex"));
        SET_STRING_ELT(cn, 1, Rf_mkChar("y_value"));
        SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dn, 0, R_NilValue);
        SET_VECTOR_ELT(dn, 1, cn);
        Rf_setAttrib(r_basin_bd_df, R_DimNamesSymbol, dn);
        SET_VECTOR_ELT(r_basin, idx++, r_basin_bd_df);
        UNPROTECT(3);
    }

    // terminal_extrema
    {
        SEXP r_term = PROTECT(Rf_allocVector(INTSXP, static_cast<int>(basin.terminal_extrema.size())));
        int* p = INTEGER(r_term);
        for (size_t i = 0; i < basin.terminal_extrema.size(); ++i) {
            p[i] = static_cast<int>(basin.terminal_extrema[i] + 1);
        }
        SET_VECTOR_ELT(r_basin, idx++, r_term);
        UNPROTECT(1);
    }

    // all_predecessors
    {
        SEXP r_all_pred;
        if (with_trajectories) {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, static_cast<int>(n_vertices)));
            for (size_t v = 0; v < n_vertices; ++v) {
                auto it = basin.all_predecessors.find(v);
                if (it != basin.all_predecessors.end() && !it->second.empty()) {
                    const auto& preds = it->second;
                    SEXP r_preds = PROTECT(Rf_allocVector(INTSXP, static_cast<int>(preds.size())));
                    int* p = INTEGER(r_preds);
                    for (size_t j = 0; j < preds.size(); ++j) {
                        p[j] = static_cast<int>(preds[j] + 1);
                    }
                    SET_VECTOR_ELT(r_all_pred, static_cast<int>(v), r_preds);
                    UNPROTECT(1);
                } else {
                    SET_VECTOR_ELT(r_all_pred, static_cast<int>(v), Rf_allocVector(INTSXP, 0));
                }
            }
        } else {
            r_all_pred = PROTECT(Rf_allocVector(VECSXP, 0));
        }
        SET_VECTOR_ELT(r_basin, idx++, r_all_pred);
        UNPROTECT(1);
    }

    if (with_trajectories) {
        SEXP r_traj_sets = PROTECT(Rf_allocVector(VECSXP, static_cast<int>(basin.trajectory_sets.size())));
        for (size_t i = 0; i < basin.trajectory_sets.size(); ++i) {
            const auto& ts = basin.trajectory_sets[i];
            SEXP r_tset = PROTECT(Rf_allocVector(VECSXP, 2));
            SEXP r_tnames = PROTECT(Rf_allocVector(STRSXP, 2));
            SET_STRING_ELT(r_tnames, 0, Rf_mkChar("terminal_vertex"));
            SET_STRING_ELT(r_tnames, 1, Rf_mkChar("trajectories"));
            Rf_setAttrib(r_tset, R_NamesSymbol, r_tnames);

            SET_VECTOR_ELT(r_tset, 0, Rf_ScalarInteger(static_cast<int>(ts.terminal_vertex + 1)));

            SEXP r_trajs = PROTECT(Rf_allocVector(VECSXP, static_cast<int>(ts.trajectories.size())));
            for (size_t j = 0; j < ts.trajectories.size(); ++j) {
                const auto& traj = ts.trajectories[j];
                SEXP r_traj = PROTECT(Rf_allocVector(INTSXP, static_cast<int>(traj.size())));
                int* p = INTEGER(r_traj);
                for (size_t k = 0; k < traj.size(); ++k) {
                    p[k] = static_cast<int>(traj[k] + 1);
                }
                SET_VECTOR_ELT(r_trajs, static_cast<int>(j), r_traj);
                UNPROTECT(1);
            }
            SET_VECTOR_ELT(r_tset, 1, r_trajs);

            SET_VECTOR_ELT(r_traj_sets, static_cast<int>(i), r_tset);
            UNPROTECT(3);
        }
        SET_VECTOR_ELT(r_basin, idx++, r_traj_sets);
        UNPROTECT(1);
    }

    UNPROTECT(2);
    return r_basin;
}

rtcb_search_result_t run_rtcb_search(
    const set_wgraph_t& graph,
    const std::vector<std::vector<arc_t>>& out_arcs,
    const std::vector<bool>& candidate_mask,
    const std::vector<size_t>& sinks,
    const std::vector<double>& y,
    const std::vector<double>& tau,
    size_t source,
    bool detect_maxima,
    const rtcb_params_t& params,
    double lambda
    ) {
    const size_t n = graph.adjacency_list.size();
    rtcb_search_result_t out;
    out.labels_at_vertex.resize(n);

    // Priority queue over labels by shortest path length first.
    using pq_entry_t = std::pair<double, int>;
    std::priority_queue<pq_entry_t, std::vector<pq_entry_t>, std::greater<pq_entry_t>> pq;

    rtcb_label_t init;
    init.vertex = source;
    init.d = 0.0;
    init.G = 0.0;
    init.L = 0.0;
    init.base = 0.0;
    init.run = 0;
    init.hops = 0;
    init.pred = -1;
    init.active = true;

    out.labels.push_back(init);
    out.labels_at_vertex[source].push_back(0);
    pq.push({0.0, 0});

    auto keep_best_k = [&](size_t v) {
        std::vector<int> active_ids;
        active_ids.reserve(out.labels_at_vertex[v].size());
        for (int id : out.labels_at_vertex[v]) {
            if (id >= 0 && static_cast<size_t>(id) < out.labels.size() && out.labels[id].active) {
                active_ids.push_back(id);
            }
        }
        if (static_cast<int>(active_ids.size()) <= params.k_max) {
            return;
        }
        std::sort(active_ids.begin(), active_ids.end(), [&](int a, int b) {
            const auto& A = out.labels[a];
            const auto& B = out.labels[b];
            if (A.base != B.base) return A.base > B.base;
            if (A.d != B.d) return A.d < B.d;
            if (A.run != B.run) return A.run < B.run;
            return A.hops < B.hops;
        });
        for (size_t i = static_cast<size_t>(params.k_max); i < active_ids.size(); ++i) {
            out.labels[active_ids[i]].active = false;
        }
    };

    while (!pq.empty()) {
        auto [queued_d, label_id] = pq.top();
        pq.pop();

        if (label_id < 0 || static_cast<size_t>(label_id) >= out.labels.size()) {
            continue;
        }
        const auto lab = out.labels[static_cast<size_t>(label_id)];
        if (!lab.active) {
            continue;
        }
        if (queued_d > lab.d + params.eps_num) {
            continue;
        }

        const size_t u = lab.vertex;
        if (!candidate_mask[u]) {
            continue;
        }

        for (const auto& arc : out_arcs[u]) {
            const size_t v = arc.to;
            if (!candidate_mask[v]) {
                continue;
            }

            const double sigma = detect_maxima ? -1.0 : 1.0;
            const double delta = sigma * (y[v] - y[u]);

            double g = 0.0;
            double loss = 0.0;
            int run_next = 0;

            if (std::abs(delta) <= params.eta_step) {
                run_next = 0;
            } else if (delta < -params.eta_step) {
                g = 0.0;
                loss = -delta;
                run_next = lab.run + 1;
            } else {
                g = delta;
                loss = 0.0;
                run_next = 0;
            }

            const double G_next = lab.G + g;
            const double L_next = lab.L + loss;
            const double base_next = lambda * G_next - L_next;
            const double d_next = lab.d + arc.w;
            const int hops_next = lab.hops + 1;
            const double slack_next = base_next + tau[v];

            if (slack_next < -params.eps_num) {
                continue;
            }
            if (run_next > params.run_max) {
                continue;
            }
            if (hops_next > params.h_max) {
                continue;
            }
            if (d_next > params.d_path_max) {
                continue;
            }

            rtcb_label_t cand;
            cand.vertex = v;
            cand.d = d_next;
            cand.G = G_next;
            cand.L = L_next;
            cand.base = base_next;
            cand.run = run_next;
            cand.hops = hops_next;
            cand.pred = label_id;
            cand.active = true;

            bool rejected = false;
            for (int existing_id : out.labels_at_vertex[v]) {
                if (existing_id < 0 || static_cast<size_t>(existing_id) >= out.labels.size()) {
                    continue;
                }
                const auto& ex = out.labels[static_cast<size_t>(existing_id)];
                if (!ex.active) {
                    continue;
                }
                if (dominates_label(ex, cand, params.eps_num)) {
                    rejected = true;
                    break;
                }
            }
            if (rejected) {
                continue;
            }

            for (int existing_id : out.labels_at_vertex[v]) {
                if (existing_id < 0 || static_cast<size_t>(existing_id) >= out.labels.size()) {
                    continue;
                }
                auto& ex = out.labels[static_cast<size_t>(existing_id)];
                if (!ex.active) {
                    continue;
                }
                if (dominates_label(cand, ex, params.eps_num)) {
                    ex.active = false;
                }
            }

            const int new_id = static_cast<int>(out.labels.size());
            out.labels.push_back(cand);
            out.labels_at_vertex[v].push_back(new_id);
            keep_best_k(v);
            if (out.labels[static_cast<size_t>(new_id)].active) {
                pq.push({d_next, new_id});
            }
        }
    }

    for (size_t sink_v : sinks) {
        auto& ids = out.labels_at_vertex[sink_v];
        std::vector<int> accepted;
        accepted.reserve(ids.size());
        for (int id : ids) {
            if (id < 0 || static_cast<size_t>(id) >= out.labels.size()) {
                continue;
            }
            const auto& lab = out.labels[static_cast<size_t>(id)];
            if (lab.active && lab.base >= -params.eps_num) {
                accepted.push_back(id);
            }
        }
        if (accepted.empty()) {
            continue;
        }
        std::sort(accepted.begin(), accepted.end(), [&](int a, int b) {
            const auto& A = out.labels[static_cast<size_t>(a)];
            const auto& B = out.labels[static_cast<size_t>(b)];
            if (A.base != B.base) return A.base > B.base;
            if (A.d != B.d) return A.d < B.d;
            return A.hops < B.hops;
        });
        if (static_cast<int>(accepted.size()) > params.max_paths_per_sink) {
            accepted.resize(static_cast<size_t>(params.max_paths_per_sink));
        }
        out.accepted_sink_labels[sink_v] = std::move(accepted);
    }

    return out;
}

gradient_basin_t compute_rtcb_basin(
    const set_wgraph_t& graph,
    size_t source,
    const std::vector<double>& y,
    bool detect_maxima,
    double edge_length_thld,
    bool with_trajectories,
    const rtcb_params_t& params
    ) {
    const size_t n = graph.adjacency_list.size();

    gradient_basin_t basin;
    basin.vertex = source;
    basin.value = y[source];
    basin.is_maximum = detect_maxima;
    basin.hop_idx = std::numeric_limits<size_t>::max();

    if (n == 0 || source >= n) {
        return basin;
    }

    // -------------------------------------------------------------------------
    // Stage A: Candidate region via Dijkstra until n_min settled vertices.
    // -------------------------------------------------------------------------
    std::vector<double> d_in(n, std::numeric_limits<double>::infinity());
    std::vector<bool> settled(n, false);

    using dijkstra_entry_t = std::pair<double, size_t>;
    std::priority_queue<dijkstra_entry_t, std::vector<dijkstra_entry_t>, std::greater<dijkstra_entry_t>> pq0;
    d_in[source] = 0.0;
    pq0.push({0.0, source});

    size_t settled_count = 0;
    while (!pq0.empty() && settled_count < params.n_min) {
        auto [du, u] = pq0.top();
        pq0.pop();

        if (settled[u]) {
            continue;
        }
        if (du > d_in[u] + params.eps_num) {
            continue;
        }

        settled[u] = true;
        settled_count++;

        for (const auto& edge : graph.adjacency_list[u]) {
            if (edge.weight > edge_length_thld + params.eps_num) {
                continue;
            }
            const size_t v = edge.vertex;
            const double nd = du + edge.weight;
            if (nd + params.eps_num < d_in[v]) {
                d_in[v] = nd;
                pq0.push({nd, v});
            }
        }
    }

    // If the reachable component is smaller than n_min, include all reachable nodes.
    if (settled_count < params.n_min) {
        while (!pq0.empty()) {
            auto [du, u] = pq0.top();
            pq0.pop();
            if (settled[u]) {
                continue;
            }
            if (du > d_in[u] + params.eps_num) {
                continue;
            }
            settled[u] = true;
            settled_count++;
            for (const auto& edge : graph.adjacency_list[u]) {
                if (edge.weight > edge_length_thld + params.eps_num) {
                    continue;
                }
                const size_t v = edge.vertex;
                const double nd = du + edge.weight;
                if (nd + params.eps_num < d_in[v]) {
                    d_in[v] = nd;
                    pq0.push({nd, v});
                }
            }
        }
    }

    if (!settled[source]) {
        return basin;
    }

    std::vector<bool> candidate_mask = settled;
    std::vector<size_t> candidate_vertices;
    candidate_vertices.reserve(n);
    for (size_t v = 0; v < n; ++v) {
        if (candidate_mask[v]) {
            candidate_vertices.push_back(v);
        }
    }
    if (candidate_vertices.empty()) {
        return basin;
    }

    const double lambda = compute_lambda(params);

    rtcb_search_result_t search;
    std::vector<size_t> sinks;
    std::vector<std::vector<arc_t>> out_arcs(n), in_arcs(n);
    std::vector<double> d_out(n, std::numeric_limits<double>::infinity());
    std::vector<double> tau(n, 0.0);

    const int n_iters = params.sink_prune ? params.sink_prune_max_iter : 1;
    for (int iter = 0; iter < n_iters; ++iter) {
        // Build DAG and reverse DAG on current candidate set.
        for (auto& vec : out_arcs) vec.clear();
        for (auto& vec : in_arcs) vec.clear();
        sinks.clear();
        std::fill(d_out.begin(), d_out.end(), std::numeric_limits<double>::infinity());
        std::fill(tau.begin(), tau.end(), 0.0);

        for (size_t u : candidate_vertices) {
            if (!candidate_mask[u]) continue;
            for (const auto& edge : graph.adjacency_list[u]) {
                const size_t v = edge.vertex;
                if (!candidate_mask[v]) continue;
                if (edge.weight > edge_length_thld + params.eps_num) continue;

                if (d_in[v] > d_in[u] + params.epsilon_d) {
                    out_arcs[u].push_back({v, edge.weight});
                    in_arcs[v].push_back({u, edge.weight});
                }
            }
        }

        for (size_t u : candidate_vertices) {
            if (!candidate_mask[u]) continue;
            if (out_arcs[u].empty()) {
                sinks.push_back(u);
            }
        }
        if (sinks.empty()) {
            sinks.push_back(source);
        }

        // Distance-to-sink via reverse Dijkstra on DAG.
        using out_entry_t = std::pair<double, size_t>;
        std::priority_queue<out_entry_t, std::vector<out_entry_t>, std::greater<out_entry_t>> pq_out;
        for (size_t s : sinks) {
            d_out[s] = 0.0;
            pq_out.push({0.0, s});
        }
        while (!pq_out.empty()) {
            auto [du, u] = pq_out.top();
            pq_out.pop();
            if (du > d_out[u] + params.eps_num) continue;
            for (const auto& in_edge : in_arcs[u]) {
                const size_t p = in_edge.to;
                const double nd = du + in_edge.w;
                if (nd + params.eps_num < d_out[p]) {
                    d_out[p] = nd;
                    pq_out.push({nd, p});
                }
            }
        }

        // Local tolerance schedule tau(v) = tau0 * (1 - p(v))^kappa.
        for (size_t v : candidate_vertices) {
            if (!candidate_mask[v]) continue;
            double progress = 0.0;
            if (std::isfinite(d_out[v])) {
                const double denom = d_in[v] + d_out[v] + params.epsilon_d;
                progress = (denom > 0.0) ? (d_in[v] / denom) : 0.0;
            } else {
                progress = 0.0;
            }
            progress = std::max(0.0, std::min(1.0, progress));
            tau[v] = params.tau0 * std::pow(std::max(0.0, 1.0 - progress), params.kappa);
        }

        search = run_rtcb_search(
            graph, out_arcs, candidate_mask, sinks, y, tau, source, detect_maxima, params, lambda
        );

        if (!params.sink_prune) {
            break;
        }

        // Optional fixed-point sink pruning.
        bool changed = false;
        for (size_t sink_v : sinks) {
            if (search.accepted_sink_labels.find(sink_v) == search.accepted_sink_labels.end()) {
                if (sink_v != source && candidate_mask[sink_v]) {
                    candidate_mask[sink_v] = false;
                    changed = true;
                }
            }
        }
        if (!changed) {
            break;
        }
    }

    // -------------------------------------------------------------------------
    // Stage C: Basin extraction from accepted sink labels.
    // -------------------------------------------------------------------------
    std::unordered_set<size_t> basin_vertices;
    std::unordered_map<size_t, size_t> hop_dist_map;
    std::unordered_map<size_t, std::vector<size_t>> all_predecessors;
    std::unordered_set<size_t> terminal_set;
    std::unordered_map<size_t, std::vector<std::vector<size_t>>> trajectories_by_sink;

    for (const auto& kv : search.accepted_sink_labels) {
        const size_t sink_v = kv.first;
        terminal_set.insert(sink_v);
        for (int label_id : kv.second) {
            if (label_id < 0 || static_cast<size_t>(label_id) >= search.labels.size()) {
                continue;
            }
            int curr = label_id;
            std::vector<size_t> trajectory;
            while (curr >= 0) {
                const auto& lab = search.labels[static_cast<size_t>(curr)];
                const size_t v = lab.vertex;
                basin_vertices.insert(v);
                auto it_h = hop_dist_map.find(v);
                if (it_h == hop_dist_map.end() || static_cast<size_t>(lab.hops) < it_h->second) {
                    hop_dist_map[v] = static_cast<size_t>(lab.hops);
                }
                trajectory.push_back(v);

                if (lab.pred >= 0) {
                    const size_t pv = search.labels[static_cast<size_t>(lab.pred)].vertex;
                    all_predecessors[v].push_back(pv);
                }
                curr = lab.pred;
            }
            if (with_trajectories && !trajectory.empty()) {
                trajectories_by_sink[sink_v].push_back(std::move(trajectory));
            }
        }
    }

    if (basin_vertices.empty()) {
        // Fallback: keep at least the source as a valid singleton basin.
        basin_vertices.insert(source);
        hop_dist_map[source] = 0;
        terminal_set.insert(source);
    }

    deduplicate_predecessors(all_predecessors);

    basin.hop_dist_map = std::move(hop_dist_map);
    basin.all_predecessors = std::move(all_predecessors);

    size_t max_hop = 0;
    for (const auto& kv : basin.hop_dist_map) {
        if (kv.second > max_hop) max_hop = kv.second;
    }
    basin.hop_idx = max_hop;

    basin.terminal_extrema.assign(terminal_set.begin(), terminal_set.end());
    std::sort(basin.terminal_extrema.begin(), basin.terminal_extrema.end());

    if (with_trajectories) {
        basin.trajectory_sets.clear();
        basin.trajectory_sets.reserve(trajectories_by_sink.size());
        for (auto& kv : trajectories_by_sink) {
            trajectory_set_t tset;
            tset.terminal_vertex = kv.first;
            tset.trajectories = std::move(kv.second);
            basin.trajectory_sets.push_back(std::move(tset));
        }
        std::sort(basin.trajectory_sets.begin(), basin.trajectory_sets.end(),
                  [](const trajectory_set_t& a, const trajectory_set_t& b) {
                      return a.terminal_vertex < b.terminal_vertex;
                  });
    }

    // Boundary map: neighbors of basin vertices that are outside the basin.
    for (size_t u : basin_vertices) {
        for (const auto& edge : graph.adjacency_list[u]) {
            const size_t v = edge.vertex;
            if (basin_vertices.find(v) == basin_vertices.end()) {
                basin.y_nbhd_bd_map[v] = y[v];
            }
        }
    }

    return basin;
}

} // anonymous namespace

extern "C" SEXP S_compute_basins_of_attraction_rtcb(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_edge_length_quantile_thld,
    SEXP s_with_trajectories,
    SEXP s_params
    ) {
    std::vector<std::vector<int>> adj_list = convert_adj_list_from_R(s_adj_list);
    std::vector<std::vector<double>> weight_list = convert_weight_list_from_R(s_weight_list);

    std::vector<double> y(REAL(s_y), REAL(s_y) + LENGTH(s_y));
    const size_t n = y.size();

    if (!Rf_isReal(s_edge_length_quantile_thld) || LENGTH(s_edge_length_quantile_thld) != 1) {
        Rf_error("edge_length_quantile_thld must be a single numeric value");
    }
    const double edge_length_quantile_thld = REAL(s_edge_length_quantile_thld)[0];
    if (std::isnan(edge_length_quantile_thld) || edge_length_quantile_thld < 0.0) {
        Rf_error("edge_length_quantile_thld must be non-negative and finite");
    }

    const bool with_trajectories = (Rf_asLogical(s_with_trajectories) != 0);
    const rtcb_params_t params = parse_rtcb_params(s_params, n);

    set_wgraph_t graph(adj_list, weight_list);
    const double edge_length_thld = graph.compute_quantile_edge_length(edge_length_quantile_thld);

    // Identify strict local extrema (same convention as geodesic backend).
    std::vector<size_t> local_minima;
    std::vector<size_t> local_maxima;
    local_minima.reserve(n);
    local_maxima.reserve(n);

    for (size_t v = 0; v < n; ++v) {
        bool is_local_min = true;
        bool is_local_max = true;
        for (const auto& edge : graph.adjacency_list[v]) {
            const size_t u = edge.vertex;
            if (y[u] >= y[v]) is_local_max = false;
            if (y[u] <= y[v]) is_local_min = false;
        }
        if (is_local_min) local_minima.push_back(v);
        if (is_local_max) local_maxima.push_back(v);
    }

    std::vector<gradient_basin_t> min_basins;
    std::vector<gradient_basin_t> max_basins;
    min_basins.reserve(local_minima.size());
    max_basins.reserve(local_maxima.size());

    for (size_t m : local_minima) {
        gradient_basin_t basin = compute_rtcb_basin(
            graph, m, y, false, edge_length_thld, with_trajectories, params
        );
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            min_basins.push_back(std::move(basin));
        }
    }
    for (size_t M : local_maxima) {
        gradient_basin_t basin = compute_rtcb_basin(
            graph, M, y, true, edge_length_thld, with_trajectories, params
        );
        if (basin.hop_idx != std::numeric_limits<size_t>::max()) {
            max_basins.push_back(std::move(basin));
        }
    }

    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP r_names = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(r_names, 0, Rf_mkChar("lmin_basins"));
    SET_STRING_ELT(r_names, 1, Rf_mkChar("lmax_basins"));
    Rf_setAttrib(r_result, R_NamesSymbol, r_names);

    SEXP r_min = PROTECT(Rf_allocVector(VECSXP, static_cast<int>(min_basins.size())));
    for (size_t i = 0; i < min_basins.size(); ++i) {
        SET_VECTOR_ELT(r_min, static_cast<int>(i), basin_to_r(min_basins[i], with_trajectories, n));
    }
    SET_VECTOR_ELT(r_result, 0, r_min);

    SEXP r_max = PROTECT(Rf_allocVector(VECSXP, static_cast<int>(max_basins.size())));
    for (size_t i = 0; i < max_basins.size(); ++i) {
        SET_VECTOR_ELT(r_max, static_cast<int>(i), basin_to_r(max_basins[i], with_trajectories, n));
    }
    SET_VECTOR_ELT(r_result, 1, r_max);

    UNPROTECT(4);
    return r_result;
}
