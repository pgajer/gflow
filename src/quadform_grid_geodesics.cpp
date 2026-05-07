#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct GridGraph {
    int grid_size;
    int index_k;
    double c1;
    double c2;
    std::string domain_shape;
    double domain_radius;
    double step;
    int n_vertices;
    int n_edges;
    std::vector<int> cell_id;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<std::vector<int> > adj;
    std::vector<std::vector<double> > weight;
};

double quadform_edge_length_2d(double ux, double uy,
                               double vx, double vy,
                               int index_k,
                               double c1,
                               double c2) {
    const double dx = vx - ux;
    const double dy = vy - uy;
    const double A = dx * dx + dy * dy;
    if (A == 0.0) {
        return 0.0;
    }
    const double s1 = index_k >= 1 ? 1.0 : -1.0;
    const double s2 = index_k >= 2 ? 1.0 : -1.0;
    const double a = 2.0 * (ux * s1 * c1 * dx + uy * s2 * c2 * dy);
    const double b = 2.0 * (dx * s1 * c1 * dx + dy * s2 * c2 * dy);
    if (std::abs(b) <= std::sqrt(std::numeric_limits<double>::epsilon())) {
        return std::sqrt(A + a * a);
    }
    const double sqrt_A = std::sqrt(A);
    auto antiderivative = [A, sqrt_A](double z) {
        return 0.5 * (z * std::sqrt(A + z * z) + A * std::asinh(z / sqrt_A));
    };
    return (antiderivative(a + b) - antiderivative(a)) / b;
}

int cell_index(int row, int col, int grid_size) {
    return row * grid_size + col;
}

GridGraph build_quadform_grid_graph(int index_k,
                                    const Rcpp::NumericVector& coefficients,
                                    const std::string& domain_shape,
                                    double domain_radius,
                                    int grid_size) {
    if (index_k < 0 || index_k > 2) {
        throw std::runtime_error("'index.k' must be an integer between 0 and 2.");
    }
    if (!std::isfinite(domain_radius) || domain_radius <= 0.0) {
        throw std::runtime_error("'domain.radius' must be a positive finite numeric scalar.");
    }
    if (grid_size < 5) {
        throw std::runtime_error("'grid.size' must be an integer at least 5.");
    }
    if (coefficients.size() != 2 ||
        !std::isfinite(coefficients[0]) || !std::isfinite(coefficients[1]) ||
        coefficients[0] <= 0.0 || coefficients[1] <= 0.0) {
        throw std::runtime_error("'coefficients' must be a positive finite numeric vector of length 2.");
    }
    if (domain_shape != "disk" && domain_shape != "square") {
        throw std::runtime_error("'domain.shape' must be either 'disk' or 'square'.");
    }

    GridGraph g;
    g.grid_size = grid_size;
    g.index_k = index_k;
    g.c1 = coefficients[0];
    g.c2 = coefficients[1];
    g.domain_shape = domain_shape;
    g.domain_radius = domain_radius;
    g.step = 2.0 * domain_radius / static_cast<double>(grid_size - 1);
    g.cell_id.assign(static_cast<std::size_t>(grid_size) * grid_size, -1);

    const double tol_radius2 = domain_radius * domain_radius * (1.0 + 1e-12);
    for (int row = 0; row < grid_size; ++row) {
        const double x = -domain_radius + g.step * static_cast<double>(row);
        for (int col = 0; col < grid_size; ++col) {
            const double y = -domain_radius + g.step * static_cast<double>(col);
            if (domain_shape == "square" || x * x + y * y <= tol_radius2) {
                const int id = static_cast<int>(g.x.size());
                g.cell_id[static_cast<std::size_t>(cell_index(row, col, grid_size))] = id;
                g.x.push_back(x);
                g.y.push_back(y);
            }
        }
    }
    g.n_vertices = static_cast<int>(g.x.size());
    g.adj.assign(g.n_vertices, std::vector<int>());
    g.weight.assign(g.n_vertices, std::vector<double>());

    const int offsets[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
    int edge_count = 0;
    for (int row = 0; row < grid_size; ++row) {
        for (int col = 0; col < grid_size; ++col) {
            const int from = g.cell_id[static_cast<std::size_t>(cell_index(row, col, grid_size))];
            if (from < 0) {
                continue;
            }
            for (const auto& offset : offsets) {
                const int row2 = row + offset[0];
                const int col2 = col + offset[1];
                if (row2 < 0 || row2 >= grid_size || col2 < 0 || col2 >= grid_size) {
                    continue;
                }
                const int to = g.cell_id[static_cast<std::size_t>(cell_index(row2, col2, grid_size))];
                if (to < 0) {
                    continue;
                }
                const double w = quadform_edge_length_2d(
                    g.x[static_cast<std::size_t>(from)], g.y[static_cast<std::size_t>(from)],
                    g.x[static_cast<std::size_t>(to)], g.y[static_cast<std::size_t>(to)],
                    index_k,
                    g.c1,
                    g.c2
                );
                g.adj[static_cast<std::size_t>(from)].push_back(to);
                g.weight[static_cast<std::size_t>(from)].push_back(w);
                g.adj[static_cast<std::size_t>(to)].push_back(from);
                g.weight[static_cast<std::size_t>(to)].push_back(w);
                ++edge_count;
            }
        }
    }
    g.n_edges = edge_count;
    return g;
}

bool point_inside_domain(double x, double y, const GridGraph& g, double tol = 1e-8) {
    if (g.domain_shape == "disk") {
        return x * x + y * y <= g.domain_radius * g.domain_radius * (1.0 + tol);
    }
    return std::max(std::abs(x), std::abs(y)) <= g.domain_radius * (1.0 + tol);
}

int nearest_grid_vertex(const GridGraph& g, double qx, double qy, double* snap_distance) {
    if (!point_inside_domain(qx, qy, g)) {
        throw std::runtime_error("All query points must lie inside the parameter domain.");
    }
    int row0 = static_cast<int>(std::llround((qx + g.domain_radius) / g.step));
    int col0 = static_cast<int>(std::llround((qy + g.domain_radius) / g.step));
    row0 = std::max(0, std::min(g.grid_size - 1, row0));
    col0 = std::max(0, std::min(g.grid_size - 1, col0));

    double best_d2 = std::numeric_limits<double>::infinity();
    int best = -1;
    for (int radius = 0; radius < g.grid_size; ++radius) {
        bool saw_cell = false;
        for (int row = std::max(0, row0 - radius);
             row <= std::min(g.grid_size - 1, row0 + radius); ++row) {
            for (int col = std::max(0, col0 - radius);
                 col <= std::min(g.grid_size - 1, col0 + radius); ++col) {
                if (std::max(std::abs(row - row0), std::abs(col - col0)) != radius) {
                    continue;
                }
                const int id = g.cell_id[static_cast<std::size_t>(cell_index(row, col, g.grid_size))];
                if (id < 0) {
                    continue;
                }
                saw_cell = true;
                const double dx = g.x[static_cast<std::size_t>(id)] - qx;
                const double dy = g.y[static_cast<std::size_t>(id)] - qy;
                const double d2 = dx * dx + dy * dy;
                if (d2 < best_d2 || (d2 == best_d2 && id < best)) {
                    best_d2 = d2;
                    best = id;
                }
            }
        }
        if (best >= 0 && saw_cell) {
            break;
        }
    }
    if (best < 0) {
        throw std::runtime_error("Internal error: failed to find nearest grid vertex.");
    }
    if (snap_distance != nullptr) {
        *snap_distance = std::sqrt(best_d2);
    }
    return best;
}

std::vector<double> dijkstra_to_targets(const std::vector<std::vector<int> >& adj,
                                        const std::vector<std::vector<double> >& weight,
                                        int source,
                                        const std::vector<int>& targets) {
    const int n = static_cast<int>(adj.size());
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(static_cast<std::size_t>(n), inf);
    std::vector<char> settled(static_cast<std::size_t>(n), 0);
    std::vector<int> target_remaining(static_cast<std::size_t>(n), 0);
    int remaining = 0;
    for (int target : targets) {
        if (target < 0 || target >= n) {
            throw std::runtime_error("Internal error: target index is out of range.");
        }
        if (target_remaining[static_cast<std::size_t>(target)] == 0) {
            ++remaining;
        }
        ++target_remaining[static_cast<std::size_t>(target)];
    }

    typedef std::pair<double, int> QueueItem;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem> > pq;
    dist[static_cast<std::size_t>(source)] = 0.0;
    pq.push(QueueItem(0.0, source));
    while (!pq.empty() && remaining > 0) {
        const QueueItem item = pq.top();
        pq.pop();
        const double du = item.first;
        const int u = item.second;
        if (settled[static_cast<std::size_t>(u)]) {
            continue;
        }
        settled[static_cast<std::size_t>(u)] = 1;
        if (target_remaining[static_cast<std::size_t>(u)] > 0) {
            --remaining;
        }
        const std::vector<int>& nbrs = adj[static_cast<std::size_t>(u)];
        const std::vector<double>& weights = weight[static_cast<std::size_t>(u)];
        for (std::size_t pos = 0; pos < nbrs.size(); ++pos) {
            const int v = nbrs[pos];
            if (settled[static_cast<std::size_t>(v)]) {
                continue;
            }
            const double nd = du + weights[pos];
            if (nd < dist[static_cast<std::size_t>(v)]) {
                dist[static_cast<std::size_t>(v)] = nd;
                pq.push(QueueItem(nd, v));
            }
        }
    }

    std::vector<double> out;
    out.reserve(targets.size());
    for (int target : targets) {
        out.push_back(dist[static_cast<std::size_t>(target)]);
    }
    return out;
}

struct DijkstraResult {
    std::vector<double> dist;
    std::vector<int> pred;
};

DijkstraResult dijkstra_with_predecessor(const std::vector<std::vector<int> >& adj,
                                         const std::vector<std::vector<double> >& weight,
                                         int source,
                                         const std::vector<int>& targets) {
    const int n = static_cast<int>(adj.size());
    const double inf = std::numeric_limits<double>::infinity();
    DijkstraResult result;
    result.dist.assign(static_cast<std::size_t>(n), inf);
    result.pred.assign(static_cast<std::size_t>(n), -1);
    std::vector<char> settled(static_cast<std::size_t>(n), 0);
    std::vector<int> target_remaining(static_cast<std::size_t>(n), 0);
    int remaining = 0;
    for (int target : targets) {
        if (target_remaining[static_cast<std::size_t>(target)] == 0) {
            ++remaining;
        }
        ++target_remaining[static_cast<std::size_t>(target)];
    }

    typedef std::pair<double, int> QueueItem;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem> > pq;
    result.dist[static_cast<std::size_t>(source)] = 0.0;
    pq.push(QueueItem(0.0, source));
    while (!pq.empty() && remaining > 0) {
        const QueueItem item = pq.top();
        pq.pop();
        const double du = item.first;
        const int u = item.second;
        if (settled[static_cast<std::size_t>(u)]) {
            continue;
        }
        settled[static_cast<std::size_t>(u)] = 1;
        if (target_remaining[static_cast<std::size_t>(u)] > 0) {
            --remaining;
        }
        const std::vector<int>& nbrs = adj[static_cast<std::size_t>(u)];
        const std::vector<double>& weights = weight[static_cast<std::size_t>(u)];
        for (std::size_t pos = 0; pos < nbrs.size(); ++pos) {
            const int v = nbrs[pos];
            if (settled[static_cast<std::size_t>(v)]) {
                continue;
            }
            const double nd = du + weights[pos];
            if (nd < result.dist[static_cast<std::size_t>(v)]) {
                result.dist[static_cast<std::size_t>(v)] = nd;
                result.pred[static_cast<std::size_t>(v)] = u;
                pq.push(QueueItem(nd, v));
            }
        }
    }
    return result;
}

std::vector<int> reconstruct_path(int source,
                                  int target,
                                  const std::vector<int>& pred) {
    std::vector<int> path;
    if (source == target) {
        path.push_back(source);
        return path;
    }
    int current = target;
    while (current >= 0) {
        path.push_back(current);
        if (current == source) {
            std::reverse(path.begin(), path.end());
            return path;
        }
        current = pred[static_cast<std::size_t>(current)];
    }
    path.clear();
    return path;
}

std::vector<int> nearest_grid_vertices_bruteforce(const GridGraph& g,
                                                  double qx,
                                                  double qy,
                                                  int k) {
    std::priority_queue<std::pair<double, int> > heap;
    for (int id = 0; id < g.n_vertices; ++id) {
        const double dx = g.x[static_cast<std::size_t>(id)] - qx;
        const double dy = g.y[static_cast<std::size_t>(id)] - qy;
        const double d2 = dx * dx + dy * dy;
        if (static_cast<int>(heap.size()) < k) {
            heap.push(std::make_pair(d2, id));
        } else if (d2 < heap.top().first ||
                   (d2 == heap.top().first && id < heap.top().second)) {
            heap.pop();
            heap.push(std::make_pair(d2, id));
        }
    }
    std::vector<std::pair<double, int> > pairs;
    pairs.reserve(static_cast<std::size_t>(k));
    while (!heap.empty()) {
        pairs.push_back(heap.top());
        heap.pop();
    }
    std::sort(pairs.begin(), pairs.end());
    std::vector<int> out;
    out.reserve(pairs.size());
    for (const auto& p : pairs) {
        out.push_back(p.second);
    }
    return out;
}

double quadform_value_2d(double x, double y, int index_k, double c1, double c2) {
    const double s1 = index_k >= 1 ? 1.0 : -1.0;
    const double s2 = index_k >= 2 ? 1.0 : -1.0;
    return s1 * c1 * x * x + s2 * c2 * y * y;
}

void embedded_vertex_coords(int vertex,
                            const GridGraph& g,
                            const Rcpp::NumericMatrix& X,
                            int index_k,
                            double& x,
                            double& y,
                            double& z) {
    if (vertex < g.n_vertices) {
        x = g.x[static_cast<std::size_t>(vertex)];
        y = g.y[static_cast<std::size_t>(vertex)];
    } else {
        const int sample = vertex - g.n_vertices;
        x = X(sample, 0);
        y = X(sample, 1);
    }
    z = quadform_value_2d(x, y, index_k, g.c1, g.c2);
}

double embedded_distance(double ax, double ay, double az,
                         double bx, double by, double bz) {
    const double dx = bx - ax;
    const double dy = by - ay;
    const double dz = bz - az;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

struct Projection {
    double distance;
    double arclength;
};

Projection project_sample_to_embedded_path(int sample,
                                           const std::vector<int>& path,
                                           const std::vector<double>& path_s,
                                           const GridGraph& g,
                                           const Rcpp::NumericMatrix& X,
                                           int index_k) {
    double px = X(sample, 0);
    double py = X(sample, 1);
    double pz = quadform_value_2d(px, py, index_k, g.c1, g.c2);
    Projection best;
    best.distance = std::numeric_limits<double>::infinity();
    best.arclength = 0.0;
    if (path.size() == 1) {
        double ax, ay, az;
        embedded_vertex_coords(path[0], g, X, index_k, ax, ay, az);
        best.distance = embedded_distance(px, py, pz, ax, ay, az);
        best.arclength = 0.0;
        return best;
    }

    for (std::size_t r = 0; r + 1 < path.size(); ++r) {
        double ax, ay, az;
        double bx, by, bz;
        embedded_vertex_coords(path[r], g, X, index_k, ax, ay, az);
        embedded_vertex_coords(path[r + 1], g, X, index_k, bx, by, bz);
        const double vx = bx - ax;
        const double vy = by - ay;
        const double vz = bz - az;
        const double wx = px - ax;
        const double wy = py - ay;
        const double wz = pz - az;
        const double vv = vx * vx + vy * vy + vz * vz;
        double t = vv > 0.0 ? (wx * vx + wy * vy + wz * vz) / vv : 0.0;
        t = std::max(0.0, std::min(1.0, t));
        const double qx = ax + t * vx;
        const double qy = ay + t * vy;
        const double qz = az + t * vz;
        const double d = embedded_distance(px, py, pz, qx, qy, qz);
        const double arclength = path_s[r] + t * std::sqrt(vv);
        if (d < best.distance ||
            (d == best.distance && arclength < best.arclength)) {
            best.distance = d;
            best.arclength = arclength;
        }
    }
    return best;
}

struct OraclePathResult {
    double distance;
    int n_points;
    int status_code;
    std::vector<int> sample_indices;
};

OraclePathResult sample_path_oracle_distance(int source_sample,
                                             int target_sample,
                                             const std::vector<int>& path,
                                             const GridGraph& g,
                                             const Rcpp::NumericMatrix& X,
                                             int index_k,
                                             double tube_radius) {
    OraclePathResult out;
    out.distance = 0.0;
    out.n_points = 0;
    out.status_code = 3; // no path
    if (path.empty()) {
        return out;
    }
    const int n_sample = X.nrow();
    std::vector<double> path_s(path.size(), 0.0);
    for (std::size_t r = 1; r < path.size(); ++r) {
        double ax, ay, az;
        double bx, by, bz;
        embedded_vertex_coords(path[r - 1], g, X, index_k, ax, ay, az);
        embedded_vertex_coords(path[r], g, X, index_k, bx, by, bz);
        path_s[r] = path_s[r - 1] + embedded_distance(ax, ay, az, bx, by, bz);
    }

    std::vector<std::pair<double, int> > selected;
    selected.reserve(static_cast<std::size_t>(n_sample));
    for (int sample = 0; sample < n_sample; ++sample) {
        Projection p = project_sample_to_embedded_path(
            sample, path, path_s, g, X, index_k
        );
        if (p.distance <= tube_radius ||
            sample == source_sample || sample == target_sample) {
            selected.push_back(std::make_pair(p.arclength, sample));
        }
    }
    std::sort(selected.begin(), selected.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  if (a.first == b.first) {
                      return a.second < b.second;
                  }
                  return a.first < b.first;
              });

    std::vector<int> ordered;
    ordered.reserve(selected.size());
    std::vector<char> seen(static_cast<std::size_t>(n_sample), 0);
    for (const auto& item : selected) {
        const int sample = item.second;
        if (!seen[static_cast<std::size_t>(sample)]) {
            ordered.push_back(sample);
            seen[static_cast<std::size_t>(sample)] = 1;
        }
    }

    out.sample_indices = ordered;
    out.n_points = static_cast<int>(ordered.size());
    if (out.n_points < 2) {
        out.status_code = 3;
        return out;
    }
    for (std::size_t r = 0; r + 1 < ordered.size(); ++r) {
        const int a = ordered[r];
        const int b = ordered[r + 1];
        const double ax = X(a, 0);
        const double ay = X(a, 1);
        const double az = quadform_value_2d(ax, ay, index_k, g.c1, g.c2);
        const double bx = X(b, 0);
        const double by = X(b, 1);
        const double bz = quadform_value_2d(bx, by, index_k, g.c1, g.c2);
        out.distance += embedded_distance(ax, ay, az, bx, by, bz);
    }
    out.status_code = out.n_points == 2 ? 2 : 1; // direct or ok
    return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List rcpp_quadform_grid_pair_distances(int index_k,
                                             const Rcpp::NumericVector& coefficients,
                                             const std::string& domain_shape,
                                             double domain_radius,
                                             int grid_size,
                                             const Rcpp::NumericMatrix& pair_points) {
    if (pair_points.ncol() != 4) {
        Rcpp::stop("'pair.points' must have four columns: x1, y1, x2, y2.");
    }
    GridGraph g = build_quadform_grid_graph(index_k, coefficients, domain_shape,
                                            domain_radius, grid_size);
    const int n_pairs = pair_points.nrow();
    std::vector<int> source(n_pairs);
    std::vector<int> target(n_pairs);
    Rcpp::NumericVector source_snap(n_pairs);
    Rcpp::NumericVector target_snap(n_pairs);

    std::unordered_map<int, std::vector<int> > targets_by_source;
    std::unordered_map<int, std::vector<int> > rows_by_source;
    targets_by_source.reserve(static_cast<std::size_t>(n_pairs));
    rows_by_source.reserve(static_cast<std::size_t>(n_pairs));

    for (int r = 0; r < n_pairs; ++r) {
        double snap1 = 0.0;
        double snap2 = 0.0;
        source[r] = nearest_grid_vertex(g, pair_points(r, 0), pair_points(r, 1), &snap1);
        target[r] = nearest_grid_vertex(g, pair_points(r, 2), pair_points(r, 3), &snap2);
        source_snap[r] = snap1;
        target_snap[r] = snap2;
        targets_by_source[source[r]].push_back(target[r]);
        rows_by_source[source[r]].push_back(r);
    }

    Rcpp::NumericVector distances(n_pairs);
    for (const auto& item : targets_by_source) {
        const int src = item.first;
        const std::vector<int>& tgts = item.second;
        std::vector<double> d = dijkstra_to_targets(g.adj, g.weight, src, tgts);
        const std::vector<int>& rows = rows_by_source[src];
        for (std::size_t i = 0; i < rows.size(); ++i) {
            distances[rows[i]] = d[i];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("distances") = distances,
        Rcpp::Named("source_vertex") = Rcpp::wrap(source),
        Rcpp::Named("target_vertex") = Rcpp::wrap(target),
        Rcpp::Named("source_snap_distance") = source_snap,
        Rcpp::Named("target_snap_distance") = target_snap,
        Rcpp::Named("n_vertices") = g.n_vertices,
        Rcpp::Named("n_edges") = g.n_edges,
        Rcpp::Named("grid_size") = grid_size,
        Rcpp::Named("index_k") = index_k,
        Rcpp::Named("coefficients") = coefficients,
        Rcpp::Named("domain_shape") = domain_shape,
        Rcpp::Named("domain_radius") = domain_radius
    );
}

// [[Rcpp::export]]
Rcpp::List rcpp_quadform_grid_geodesic_distances(const Rcpp::NumericMatrix& X,
                                                 int index_k,
                                                 const Rcpp::NumericVector& coefficients,
                                                 const std::string& domain_shape,
                                                 double domain_radius,
                                                 int grid_size,
                                                 int sample_connection_k,
                                                 bool with_oracle,
                                                 double oracle_tube_radius,
                                                 int oracle_tube_k,
                                                 bool return_oracle_paths) {
    if (X.ncol() != 2) {
        Rcpp::stop("'X' must have two parameter-coordinate columns.");
    }
    if (sample_connection_k < 1) {
        Rcpp::stop("'sample.connection.k' must be a positive integer.");
    }
    if (with_oracle && (!std::isfinite(oracle_tube_radius) || oracle_tube_radius <= 0.0)) {
        Rcpp::stop("'oracle.tube.radius' must be a positive finite numeric scalar.");
    }
    if (with_oracle && oracle_tube_k < 1) {
        Rcpp::stop("'oracle.tube.k' must be a positive integer.");
    }
    GridGraph g = build_quadform_grid_graph(index_k, coefficients, domain_shape,
                                            domain_radius, grid_size);
    const int n_query = X.nrow();
    const int k = std::min(sample_connection_k, g.n_vertices);
    std::vector<std::vector<int> > adj = g.adj;
    std::vector<std::vector<double> > weight = g.weight;
    adj.resize(static_cast<std::size_t>(g.n_vertices + n_query));
    weight.resize(static_cast<std::size_t>(g.n_vertices + n_query));

    for (int i = 0; i < n_query; ++i) {
        const double qx = X(i, 0);
        const double qy = X(i, 1);
        if (!point_inside_domain(qx, qy, g)) {
            Rcpp::stop("All rows of 'X' must lie inside the parameter domain.");
        }
        const int qid = g.n_vertices + i;
        std::vector<int> nbrs = nearest_grid_vertices_bruteforce(g, qx, qy, k);
        for (int v : nbrs) {
            const double w = quadform_edge_length_2d(
                qx, qy,
                g.x[static_cast<std::size_t>(v)], g.y[static_cast<std::size_t>(v)],
                index_k,
                g.c1,
                g.c2
            );
            adj[static_cast<std::size_t>(qid)].push_back(v);
            weight[static_cast<std::size_t>(qid)].push_back(w);
            adj[static_cast<std::size_t>(v)].push_back(qid);
            weight[static_cast<std::size_t>(v)].push_back(w);
        }
    }

    Rcpp::NumericMatrix D(n_query, n_query);
    Rcpp::NumericMatrix oracle_D(n_query, n_query);
    Rcpp::IntegerMatrix oracle_n_points(n_query, n_query);
    Rcpp::IntegerMatrix oracle_status_code(n_query, n_query);
    Rcpp::List oracle_paths;
    if (with_oracle && return_oracle_paths) {
        oracle_paths = Rcpp::List(n_query * n_query);
    }
    std::vector<int> targets;
    targets.reserve(static_cast<std::size_t>(n_query));
    for (int j = 0; j < n_query; ++j) {
        targets.push_back(g.n_vertices + j);
    }
    for (int i = 0; i < n_query; ++i) {
        if (with_oracle) {
            DijkstraResult d = dijkstra_with_predecessor(
                adj, weight, g.n_vertices + i, targets
            );
            for (int j = 0; j < n_query; ++j) {
                const int target = g.n_vertices + j;
                D(i, j) = d.dist[static_cast<std::size_t>(target)];
                if (i == j) {
                    oracle_D(i, j) = 0.0;
                    oracle_n_points(i, j) = 1;
                    oracle_status_code(i, j) = 0; // same
                    if (return_oracle_paths) {
                        oracle_paths[i * n_query + j] = Rcpp::IntegerVector::create(i + 1);
                    }
                    continue;
                }
                std::vector<int> path = reconstruct_path(
                    g.n_vertices + i, target, d.pred
                );
                OraclePathResult oracle = sample_path_oracle_distance(
                    i, j, path, g, X, index_k, oracle_tube_radius
                );
                oracle_D(i, j) = oracle.distance;
                oracle_n_points(i, j) = oracle.n_points;
                oracle_status_code(i, j) = oracle.status_code;
                if (return_oracle_paths) {
                    Rcpp::IntegerVector sample_path(oracle.sample_indices.size());
                    for (std::size_t r = 0; r < oracle.sample_indices.size(); ++r) {
                        sample_path[static_cast<R_xlen_t>(r)] = oracle.sample_indices[r] + 1;
                    }
                    oracle_paths[i * n_query + j] = sample_path;
                }
            }
        } else {
            std::vector<double> d = dijkstra_to_targets(
                adj, weight, g.n_vertices + i, targets
            );
            for (int j = 0; j < n_query; ++j) {
                D(i, j) = d[static_cast<std::size_t>(j)];
            }
        }
    }

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("distances") = D,
        Rcpp::Named("n_sample_vertices") = n_query,
        Rcpp::Named("n_reference_vertices") = g.n_vertices,
        Rcpp::Named("n_edges") = g.n_edges + n_query * k,
        Rcpp::Named("grid_size") = grid_size,
        Rcpp::Named("sample_connection_k") = k,
        Rcpp::Named("index_k") = index_k,
        Rcpp::Named("coefficients") = coefficients,
        Rcpp::Named("domain_shape") = domain_shape,
        Rcpp::Named("domain_radius") = domain_radius
    );
    if (with_oracle) {
        out["oracle_distances"] = oracle_D;
        out["oracle_n_points"] = oracle_n_points;
        out["oracle_status_code"] = oracle_status_code;
        out["oracle_tube_radius"] = oracle_tube_radius;
        out["oracle_tube_k"] = oracle_tube_k;
        out["oracle_method"] = "sample.path";
        if (return_oracle_paths) {
            out["oracle_paths"] = oracle_paths;
        }
    }
    return out;
}
