// traj_clustering.cpp
//
// Core C++ implementation for constructing kNN-sparsified trajectory similarity
// graphs within a single (m, M) cell.

#include "traj_clustering.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

inline std::uint64_t pair_key_32(int a, int b) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
}

// Encode a directed edge (u -> v) into a 64-bit token.
inline std::uint64_t edge_token_32(int u, int v) {
    return pair_key_32(u, v);
}

static void unique_sort(std::vector<int>& x) {
    std::sort(x.begin(), x.end());
    x.erase(std::unique(x.begin(), x.end()), x.end());
}

// Build a set-representation of each trajectory (sorted unique vertices)
// and an ordered representation used for subpath hashing.
struct traj_repr_t {
    std::vector<int> ordered;  // ordered, endpoints possibly excluded
    std::vector<int> set;      // sorted unique
};

static std::vector<traj_repr_t> build_representations(const std::vector<std::vector<int>>& traj_list,
                                                      bool exclude_endpoints) {
    std::vector<traj_repr_t> reps;
    reps.resize(traj_list.size());

    for (size_t i = 0; i < traj_list.size(); ++i) {
        const auto& t = traj_list[i];
        if (t.size() < 2) {
            reps[i].ordered.clear();
            reps[i].set.clear();
            continue;
        }

        size_t start = 0;
        size_t end = t.size();
        if (exclude_endpoints && t.size() >= 2) {
            start = 1;
            end = t.size() - 1;
        }
        if (start >= end) {
            reps[i].ordered.clear();
            reps[i].set.clear();
            continue;
        }
        reps[i].ordered.assign(t.begin() + static_cast<std::ptrdiff_t>(start),
                               t.begin() + static_cast<std::ptrdiff_t>(end));
        reps[i].set = reps[i].ordered;
        unique_sort(reps[i].set);
    }

    return reps;
}

static traj_similarity_spec_t::similarity_t parse_similarity_type(const traj_similarity_spec_t& spec) {
    return spec.similarity_type;
}

static traj_similarity_spec_t::overlap_mode_t parse_overlap_mode(const traj_similarity_spec_t& spec) {
    return spec.overlap_mode;
}

static std::unordered_map<int, std::vector<int>> build_vertex_inverted_index(
    const std::vector<traj_repr_t>& reps) {

    std::unordered_map<int, std::vector<int>> idx;
    idx.reserve(reps.size() * 8);
    for (int i = 0; i < static_cast<int>(reps.size()); ++i) {
        for (int v : reps[i].set) {
            idx[v].push_back(i);
        }
    }
    return idx;
}

static void compute_df_and_weights(
    const std::unordered_map<int, std::vector<int>>& vertex_to_traj,
    int n_traj,
    double idf_smooth,
    std::unordered_map<int, int>& df,
    std::unordered_map<int, double>& w_idf) {

    df.clear();
    w_idf.clear();
    df.reserve(vertex_to_traj.size());
    w_idf.reserve(vertex_to_traj.size());

    const double n = static_cast<double>(n_traj);
    const double numer = n + 1.0;

    for (const auto& kv : vertex_to_traj) {
        int v = kv.first;
        int d = static_cast<int>(kv.second.size());
        df[v] = d;
        const double denom = static_cast<double>(d) + idf_smooth;
        double w = 0.0;
        if (denom > 0.0) {
            w = std::log(numer / denom);
            if (w < 0.0) w = 0.0;
        }
        w_idf[v] = w;
    }
}

static std::vector<double> compute_traj_total_weights(
    const std::vector<traj_repr_t>& reps,
    const std::unordered_map<int, double>& w_idf,
    bool use_idf) {

    std::vector<double> total;
    total.resize(reps.size(), 0.0);
    for (size_t i = 0; i < reps.size(); ++i) {
        double s = 0.0;
        if (use_idf) {
            for (int v : reps[i].set) {
                auto it = w_idf.find(v);
                if (it != w_idf.end()) s += it->second;
            }
        } else {
            s = static_cast<double>(reps[i].set.size());
        }
        total[i] = s;
    }
    return total;
}

// Build candidate pairs for overlap_mode = min_subpath.
// We mark (i,j) candidate if they share at least one contiguous subpath of
// length min_subpath (in vertices), i.e. length (min_subpath-1) in edges.
static std::vector<std::unordered_set<int>> build_subpath_candidates(
    const std::vector<traj_repr_t>& reps,
    int min_subpath) {

    std::vector<std::unordered_set<int>> cand;
    cand.resize(reps.size());

    if (min_subpath <= 1) {
        // Everyone is a candidate with everyone (handled elsewhere).
        return cand;
    }

    const int s_edges = std::max(1, min_subpath - 1);
    std::unordered_map<std::uint64_t, std::vector<int>> gram_to_traj;
    gram_to_traj.reserve(reps.size() * 32);

    // For each trajectory, hash each edge-gram of length s_edges.
    for (int i = 0; i < static_cast<int>(reps.size()); ++i) {
        const auto& ord = reps[i].ordered;
        if (static_cast<int>(ord.size()) < s_edges + 1) continue;

        // Build edge tokens
        std::vector<std::uint64_t> edges;
        edges.reserve(ord.size() - 1);
        for (size_t p = 0; p + 1 < ord.size(); ++p) {
            edges.push_back(edge_token_32(ord[p], ord[p + 1]));
        }
        if (static_cast<int>(edges.size()) < s_edges) continue;

        // Sliding hash (simple combinational hash; sufficient as a filter).
        for (size_t p = 0; p + static_cast<size_t>(s_edges) <= edges.size(); ++p) {
            std::uint64_t h = 1469598103934665603ULL; // FNV offset basis
            for (int q = 0; q < s_edges; ++q) {
                h ^= edges[p + static_cast<size_t>(q)];
                h *= 1099511628211ULL; // FNV prime
            }
            gram_to_traj[h].push_back(i);
        }
    }

    // Convert grams -> candidate pairs
    for (const auto& kv : gram_to_traj) {
        const auto& lst = kv.second;
        if (lst.size() < 2) continue;
        for (size_t a = 0; a < lst.size(); ++a) {
            int i = lst[a];
            for (size_t b = 0; b < lst.size(); ++b) {
                if (a == b) continue;
                cand[i].insert(lst[b]);
            }
        }
    }
    return cand;
}

static bool is_valid_symmetrize(const std::string& s) {
    return (s == "mutual" || s == "union" || s == "none");
}

static bool is_valid_knn_select(const std::string& s) {
    return (s == "weighted" || s == "raw");
}

} // namespace

traj_graph_t build_traj_knn_graph(const std::vector<std::vector<int>>& traj_list,
                                  const traj_similarity_spec_t& spec,
                                  int k,
                                  const std::string& symmetrize,
                                  const std::string& knn_select,
                                  int n_threads) {

    if (k < 1) {
        throw std::runtime_error("k must be >= 1");
    }
    if (!is_valid_symmetrize(symmetrize)) {
        throw std::runtime_error("symmetrize must be one of: mutual, union, none");
    }
    if (!is_valid_knn_select(knn_select)) {
        throw std::runtime_error("knn_select must be one of: weighted, raw");
    }

    const int n_traj = static_cast<int>(traj_list.size());
    traj_graph_t out;
    out.n_traj = n_traj;
    if (n_traj == 0) return out;

    const auto reps = build_representations(traj_list, spec.exclude_endpoints);
    auto vertex_to_traj = build_vertex_inverted_index(reps);

    std::unordered_map<int, int> df;
    std::unordered_map<int, double> w_idf;
    compute_df_and_weights(vertex_to_traj, n_traj, spec.idf_smooth, df, w_idf);

    const bool use_idf_weight = (parse_similarity_type(spec) == traj_similarity_spec_t::similarity_t::idf_jaccard);
    const auto total_weight = compute_traj_total_weights(reps, w_idf, use_idf_weight);
    const auto total_count = compute_traj_total_weights(reps, w_idf, false); // size of set

    const auto overlap_mode = parse_overlap_mode(spec);

    std::vector<std::unordered_set<int>> subpath_candidates;
    if (overlap_mode == traj_similarity_spec_t::overlap_mode_t::min_subpath) {
        subpath_candidates = build_subpath_candidates(reps, spec.min_subpath);
    }

    // Neighbor lists (directed top-k per trajectory)
    std::vector<std::vector<std::pair<int, double>>> neigh;
    neigh.resize(n_traj);
    std::vector<std::vector<std::pair<int, double>>> neigh_raw;
    if (knn_select == "raw") neigh_raw.resize(n_traj);

    // For mutual/union weights, store directed weights for lookup
    std::vector<std::unordered_map<int, double>> neigh_w;
    neigh_w.resize(n_traj);

    // Parallel over trajectories i
    int n_threads_use = n_threads;
#ifdef _OPENMP
    if (n_threads_use < 1) n_threads_use = 1;
    omp_set_num_threads(n_threads_use);
#else
    (void)n_threads_use;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n_traj; ++i) {
        std::unordered_map<int, double> inter_w;
        std::unordered_map<int, int> inter_c;
        inter_w.reserve(256);
        inter_c.reserve(256);

        // Accumulate intersections via shared vertices
        for (int v : reps[i].set) {
            auto it_list = vertex_to_traj.find(v);
            if (it_list == vertex_to_traj.end()) continue;
            const double wv = use_idf_weight ? w_idf[v] : 1.0;
            for (int j : it_list->second) {
                if (j == i) continue;
                inter_w[j] += wv;
                inter_c[j] += 1;
            }
        }

        std::vector<std::pair<int, double>> sims;
        sims.reserve(inter_w.size());
        std::vector<std::pair<int, double>> sims_raw;
        if (knn_select == "raw") sims_raw.reserve(inter_w.size());

        const double Wi = total_weight[i];
        const double Ci = total_count[i];

        for (const auto& kv : inter_w) {
            int j = kv.first;
            const int c_ij = inter_c[j];

            // Overlap criteria
            bool ok = true;
            if (overlap_mode == traj_similarity_spec_t::overlap_mode_t::any) {
                ok = (c_ij >= 1);
            } else if (overlap_mode == traj_similarity_spec_t::overlap_mode_t::min_shared) {
                ok = (c_ij >= spec.min_shared);
            } else if (overlap_mode == traj_similarity_spec_t::overlap_mode_t::min_frac) {
                const double Cj = total_count[j];
                const double denom = std::min(Ci, Cj);
                ok = (denom > 0.0) ? (static_cast<double>(c_ij) / denom >= spec.min_frac) : false;
            } else if (overlap_mode == traj_similarity_spec_t::overlap_mode_t::min_subpath) {
                // Filter by shared subpath candidates
                if (!subpath_candidates.empty()) {
                    ok = (subpath_candidates[i].find(j) != subpath_candidates[i].end());
                }
                if (ok) {
                    ok = (c_ij >= spec.min_subpath);
                }
            }
            if (!ok) continue;

            // Weighted similarity (edge weight)
            const double Wj = total_weight[j];
            const double inter = kv.second;
            const double denom_w = Wi + Wj - inter;
            double sim_w = 0.0;
            if (denom_w > 0.0) sim_w = inter / denom_w;
            if (sim_w <= 0.0) continue;
            sims.push_back({j, sim_w});

            // Raw similarity for kNN selection if requested
            if (knn_select == "raw") {
                const double Cj = total_count[j];
                const double denom_c = Ci + Cj - static_cast<double>(c_ij);
                double sim_c = 0.0;
                if (denom_c > 0.0) sim_c = static_cast<double>(c_ij) / denom_c;
                sims_raw.push_back({j, sim_c});
            }
        }

        auto cmp = [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
            return a.second > b.second;
        };

        // Select top-k indices
        if (knn_select == "weighted") {
            std::sort(sims.begin(), sims.end(), cmp);
            if (static_cast<int>(sims.size()) > k) sims.resize(static_cast<size_t>(k));
            neigh[i] = std::move(sims);
        } else {
            // raw selection: sort by sims_raw, but keep weights from sims
            // Build lookup from j -> weighted similarity
            std::unordered_map<int, double> w_lookup;
            w_lookup.reserve(sims.size());
            for (const auto& p : sims) w_lookup[p.first] = p.second;

            std::sort(sims_raw.begin(), sims_raw.end(), cmp);
            if (static_cast<int>(sims_raw.size()) > k) sims_raw.resize(static_cast<size_t>(k));

            std::vector<std::pair<int, double>> top;
            top.reserve(sims_raw.size());
            for (const auto& p : sims_raw) {
                auto itw = w_lookup.find(p.first);
                if (itw != w_lookup.end()) top.push_back({p.first, itw->second});
            }
            neigh[i] = std::move(top);
        }

        // Store directed weights for symmetrization
        std::unordered_map<int, double> wmap;
        wmap.reserve(neigh[i].size());
        for (const auto& p : neigh[i]) {
            wmap[p.first] = p.second;
        }

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            neigh_w[i] = std::move(wmap);
        }
    }

    // Symmetrize into undirected edge list with i < j
    const std::string sym = (symmetrize == "none") ? "union" : symmetrize;
    std::vector<int> from;
    std::vector<int> to;
    std::vector<double> weight;

    from.reserve(static_cast<size_t>(n_traj) * static_cast<size_t>(k));
    to.reserve(static_cast<size_t>(n_traj) * static_cast<size_t>(k));
    weight.reserve(static_cast<size_t>(n_traj) * static_cast<size_t>(k));

    std::unordered_set<std::uint64_t> seen;
    seen.reserve(static_cast<size_t>(n_traj) * static_cast<size_t>(k) * 2);

    for (int i = 0; i < n_traj; ++i) {
        for (const auto& p : neigh[i]) {
            int j = p.first;
            if (j == i) continue;
            int a = std::min(i, j);
            int b = std::max(i, j);

            bool keep = true;
            if (sym == "mutual") {
                const auto& mj = neigh_w[j];
                keep = (mj.find(i) != mj.end());
            }
            if (!keep) continue;

            const std::uint64_t kkey = pair_key_32(a, b);
            if (seen.find(kkey) != seen.end()) continue;
            seen.insert(kkey);

            // Weight: max of available directed weights
            double w_ab = 0.0;
            auto it1 = neigh_w[a].find(b);
            if (it1 != neigh_w[a].end()) w_ab = std::max(w_ab, it1->second);
            auto it2 = neigh_w[b].find(a);
            if (it2 != neigh_w[b].end()) w_ab = std::max(w_ab, it2->second);
            if (w_ab <= 0.0) continue;

            from.push_back(a);
            to.push_back(b);
            weight.push_back(w_ab);
        }
    }

    out.from = std::move(from);
    out.to = std::move(to);
    out.weight = std::move(weight);
    return out;
}

std::vector<int> leiden_cluster(const traj_graph_t& /*graph*/,
                                double /*resolution*/,
                                int /*n_iter*/,
                                int /*seed*/) {
    throw std::runtime_error("leiden_cluster() is a placeholder; run Leiden in R (igraph) for now");
}
