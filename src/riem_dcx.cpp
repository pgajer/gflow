#include "riem_dcx.hpp"
#include "kNN_r.h"            // for S_kNN()

#include <unordered_set>
#include <numeric>    // For std::accumulate

#include <ANN/ANN.h>  // ANN library header

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void riem_dcx_t::build_nerve_from_knn(
    const spmat_t& X,
    index_t k,
    index_t max_p,
    bool use_counting_measure,
    bool directed_knn,
    double density_normalization
) {
    // Extract dimensions
    const Eigen::Index n_points = X.rows();
    const Eigen::Index n_features = X.cols();

    if (n_points <= 0 || n_features <= 0) {
        Rf_error("Invalid input dimensions: n_points=%ld, n_features=%ld",
                 (long)n_points, (long)n_features);
    }
    if (k < 2) Rf_error("k must be at least 2");
    if (k >= static_cast<index_t>(n_points)) {
        Rf_error("k=%zu must be less than n_points=%ld", k, (long)n_points);
    }

    // ==================== Phase 1: kNN Computation ====================

    // Convert sparse Eigen matrix to R matrix for ANN library
    SEXP s_X = PROTECT(Rf_allocMatrix(REALSXP, n_points, n_features));
    double* p_X = REAL(s_X);

    for (Eigen::Index j = 0; j < n_features; ++j) {
        for (Eigen::Index i = 0; i < n_points; ++i) {
            p_X[i + j * n_points] = X.coeff(i, j);
        }
    }

    SEXP s_k = PROTECT(Rf_ScalarInteger(k));

    // Call S_kNN to compute neighborhoods and distances
    SEXP knn_res = PROTECT(S_kNN(s_X, s_k));
    int* indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double* distances = REAL(VECTOR_ELT(knn_res, 1));

    // Build neighborhood structure
    // knn_neighbors[i] contains the k nearest neighbors of vertex i
    std::vector<std::vector<index_t>> knn_neighbors(n_points);
    std::vector<std::vector<double>> knn_distances(n_points);

    for (Eigen::Index i = 0; i < n_points; ++i) {
        knn_neighbors[i].resize(k);
        knn_distances[i].resize(k);
        for (index_t j = 0; j < k; ++j) {
            knn_neighbors[i][j] = static_cast<index_t>(indices[i + n_points * j]);
            knn_distances[i][j] = distances[i + n_points * j];
        }
    }

    UNPROTECT(3); // knn_res, s_k, s_X

    // Enforce mutual kNN if requested
    std::vector<std::unordered_set<index_t>> neighbor_sets(n_points);
    if (!directed_knn) {
        // Build mutual kNN: keep only edges where both vertices are neighbors
        for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
            neighbor_sets[i].clear();
            for (index_t j : knn_neighbors[i]) {
                // Check if i is also in j's neighborhood
                bool is_mutual = (i == j); // vertex is always in its own closed neighborhood
                if (!is_mutual) {
                    for (index_t neighbor_of_j : knn_neighbors[j]) {
                        if (neighbor_of_j == i) {
                            is_mutual = true;
                            break;
                        }
                    }
                }
                if (is_mutual) {
                    neighbor_sets[i].insert(j);
                }
            }
            // Always include the vertex itself in its closed neighborhood
            neighbor_sets[i].insert(i);
        }
    } else {
        // Use directed kNN: include all neighbors plus self
        for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
            neighbor_sets[i].insert(i); // closed neighborhood includes self
            for (index_t j : knn_neighbors[i]) {
                neighbor_sets[i].insert(j);
            }
        }
    }

    // ==================== Phase 2: Density Computation ====================

    std::vector<double> vertex_weights(n_points);

    if (use_counting_measure) {
        // Counting measure: each vertex has unit weight
        std::fill(vertex_weights.begin(), vertex_weights.end(), 1.0);
    } else {
        // Density-based measure using d_k distances
        // w(x) = (ε + d_k(x))^{-α}
        const double epsilon = 1e-10;
        const double alpha = 1.5; // exponent in [1, 2]

        // Compute d_k for each vertex (distance to kth nearest neighbor)
        std::vector<double> d_k_values(n_points);
        for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
            d_k_values[i] = knn_distances[i][k - 1]; // last neighbor
        }

        // Compute weights
        for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
            vertex_weights[i] = std::pow(epsilon + d_k_values[i], -alpha);
        }
    }

    // Normalize weights
    double total_weight = std::accumulate(vertex_weights.begin(),
                                          vertex_weights.end(), 0.0);
    double target_sum = (density_normalization > 0.0) ?
                        density_normalization : static_cast<double>(n_points);

    if (total_weight > 1e-15) {
        double scale = target_sum / total_weight;
        for (double& w : vertex_weights) {
            w *= scale;
        }
    }

    // ==================== Phase 3: Initialize riem_dcx_t Structure ====================

    // Determine simplex counts for initialization
    std::vector<index_t> n_by_dim(max_p + 1);
    n_by_dim[0] = static_cast<index_t>(n_points);

    // We'll count edges and higher simplices as we build them
    // For now, initialize with vertices only
    init_dims(max_p, n_by_dim);

    // ==================== Phase 4: Build 0-Simplices (Vertices) ====================

    S[0].simplex_verts.resize(n_points);
    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        S[0].simplex_verts[i] = {i};
        S[0].id_of[{i}] = i;
    }

    // Initialize vertex masses: m_i = μ(N̂_k(x_i))
    // For directed kNN: this is the weighted sum over the neighborhood
    // For mutual kNN: weighted sum over mutual neighbors
    vec_t vertex_masses = vec_t::Zero(n_points);
    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        double mass = 0.0;
        for (index_t j : neighbor_sets[i]) {
            mass += vertex_weights[j];
        }
        vertex_masses[i] = mass;
    }

    g.M[0] = spmat_t(n_points, n_points);
    g.M[0].reserve(Eigen::VectorXi::Constant(n_points, 1));
    for (Eigen::Index i = 0; i < n_points; ++i) {
        g.M[0].insert(i, i) = std::max(vertex_masses[i], 1e-15);
    }
    g.M[0].makeCompressed();

    // ==================== Phase 5: Build 1-Simplices (Edges) ====================

    // Helper function to compute measure of intersection
    auto compute_intersection_measure = [&](
        const std::unordered_set<index_t>& set_i,
        const std::unordered_set<index_t>& set_j
    ) -> double {
        double measure = 0.0;
        for (index_t v : set_i) {
            if (set_j.find(v) != set_j.end()) {
                measure += vertex_weights[v];
            }
        }
        return measure;
    };

    // Helper to compute intersection as a set
    auto compute_intersection_set = [](
        const std::unordered_set<index_t>& set_i,
        const std::unordered_set<index_t>& set_j
    ) -> std::unordered_set<index_t> {
        std::unordered_set<index_t> result;
        for (index_t v : set_i) {
            if (set_j.find(v) != set_j.end()) {
                result.insert(v);
            }
        }
        return result;
    };

    // Build edges: edge (i,j) exists iff N̂_k(x_i) ∩ N̂_k(x_j) ≠ ∅
    std::vector<std::array<index_t, 2>> edge_list;
    std::vector<double> edge_weights_vec;

    for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
        for (index_t j = i + 1; j < static_cast<index_t>(n_points); ++j) {
            double intersection_measure = compute_intersection_measure(
                neighbor_sets[i], neighbor_sets[j]
            );

            if (intersection_measure > 1e-15) {
                edge_list.push_back({i, j});
                edge_weights_vec.push_back(intersection_measure);
            }
        }
    }

    const index_t n_edges = edge_list.size();

    // Update dimension counts and reinitialize
    n_by_dim[1] = n_edges;
	extend_by_one_dim((index_t)n_edges);

    // Populate edge table
    S[1].simplex_verts.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        std::vector<index_t> edge_verts = {edge_list[e][0], edge_list[e][1]};
        S[1].simplex_verts[e] = edge_verts;
        S[1].id_of[edge_verts] = e;
    }

    // Build edge masses from intersection measures
    g.M[1] = spmat_t(n_edges, n_edges);
    g.M[1].reserve(Eigen::VectorXi::Constant(n_edges, 1));
    for (index_t e = 0; e < n_edges; ++e) {
        g.M[1].insert(e, e) = std::max(edge_weights_vec[e], 1e-15);
    }
    g.M[1].makeCompressed();

    // Build boundary operator B[1]: edges → vertices
    build_incidence_from_edges();

    // Store conductances for potential L0_sym construction
    L.c1 = vec_t(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        L.c1[e] = edge_weights_vec[e];
    }

    // ==================== Phase 6: Build Higher-Dimensional Simplices ====================

    if (max_p >= 2) {
        // Build 2-simplices (triangles) and higher if requested

        for (int p = 2; p <= static_cast<int>(max_p); ++p) {
            std::vector<std::vector<index_t>> p_simplices;
            std::vector<double> p_weights;

            if (p == 2) {
                // Triangles: check each triple of vertices
                for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
                    for (index_t j = i + 1; j < static_cast<index_t>(n_points); ++j) {
                        // Skip if no edge (i,j)
                        std::vector<index_t> edge_ij = {i, j};
                        if (S[1].id_of.find(edge_ij) == S[1].id_of.end()) continue;

                        for (index_t ell = j + 1; ell < static_cast<index_t>(n_points); ++ell) {
                            // Check if all three edges exist
                            std::vector<index_t> edge_i_ell = {i, ell};
                            std::vector<index_t> edge_j_ell = {j, ell};

                            if (S[1].id_of.find(edge_i_ell) == S[1].id_of.end()) continue;
                            if (S[1].id_of.find(edge_j_ell) == S[1].id_of.end()) continue;

                            // All three edges exist, compute triple intersection
                            std::unordered_set<index_t> intersection_ij =
                                compute_intersection_set(neighbor_sets[i], neighbor_sets[j]);

                            double triple_measure = 0.0;
                            for (index_t v : intersection_ij) {
                                if (neighbor_sets[ell].find(v) != neighbor_sets[ell].end()) {
                                    triple_measure += vertex_weights[v];
                                }
                            }

                            if (triple_measure > 1e-15) {
                                p_simplices.push_back({i, j, ell});
                                p_weights.push_back(triple_measure);
                            }
                        }
                    }
                }
            } else {
                // For p > 2: build from (p-1)-simplices
                // A p-simplex exists if all its (p-1)-faces exist and neighborhoods intersect

                const auto& prev_simplices = S[p-1].simplex_verts;

                for (index_t s1 = 0; s1 < prev_simplices.size(); ++s1) {
                    for (index_t s2 = s1 + 1; s2 < prev_simplices.size(); ++s2) {
                        const auto& simplex1 = prev_simplices[s1];
                        const auto& simplex2 = prev_simplices[s2];

                        // Check if they share exactly p-1 vertices
                        std::vector<index_t> common;
                        std::set_intersection(
                            simplex1.begin(), simplex1.end(),
                            simplex2.begin(), simplex2.end(),
                            std::back_inserter(common)
                        );

                        if (common.size() == static_cast<size_t>(p - 1)) {
                            // Merge to form p-simplex
                            std::vector<index_t> merged;
                            std::set_union(
                                simplex1.begin(), simplex1.end(),
                                simplex2.begin(), simplex2.end(),
                                std::back_inserter(merged)
                            );

                            // Compute measure of (p+1)-fold intersection
                            std::unordered_set<index_t> intersection = neighbor_sets[merged[0]];
                            for (size_t idx = 1; idx < merged.size(); ++idx) {
                                std::unordered_set<index_t> temp;
                                for (index_t v : intersection) {
                                    if (neighbor_sets[merged[idx]].find(v) !=
                                        neighbor_sets[merged[idx]].end()) {
                                        temp.insert(v);
                                    }
                                }
                                intersection = std::move(temp);
                            }

                            double measure = 0.0;
                            for (index_t v : intersection) {
                                measure += vertex_weights[v];
                            }

                            if (measure > 1e-15) {
                                p_simplices.push_back(merged);
                                p_weights.push_back(measure);
                            }
                        }
                    }
                }
            }

            if (p_simplices.empty()) {
                // No p-simplices exist; stop construction
				// Adjust pmax downward - no artificial limit
                pmax = p - 1;
                break;
            }

			// Update dimension count and populate structures
            n_by_dim[p] = p_simplices.size();

            // Reinitialize with updated dimensions
			extend_by_one_dim((index_t)p_simplices.size());

            // Add p-simplices
            S[p].simplex_verts = p_simplices;
            for (index_t s = 0; s < static_cast<index_t>(p_simplices.size()); ++s) {
                S[p].id_of[p_simplices[s]] = s;
            }

            g.M[p] = spmat_t(p_simplices.size(), p_simplices.size());
            g.M[p].reserve(Eigen::VectorXi::Constant(p_simplices.size(), 1));
            for (index_t s = 0; s < static_cast<index_t>(p_simplices.size()); ++s) {
                g.M[p].insert(s, s) = std::max(p_weights[s], 1e-15);
            }
            g.M[p].makeCompressed();
        }
    }

    // ==================== Phase 7: Build Star Tables ====================

    // Populate star tables for gradient computations
    for (int p = 0; p < pmax; ++p) {
        stars[p].resize(S[p].size());

        // For each (p+1)-simplex, register it in the stars of its p-faces
        for (index_t sp1 = 0; sp1 < S[p+1].size(); ++sp1) {
            const auto& vertices = S[p+1].simplex_verts[sp1];

            // Generate all p-faces
            if (p == 0) {
                // Faces of edges are vertices
                for (index_t v : vertices) {
                    stars[0].star_over[v].push_back(sp1);
                }
            } else {
                // Generate combinations of size p+1 from vertices
                std::vector<index_t> face(p + 1);
                std::function<void(size_t, size_t)> generate_faces;
                generate_faces = [&](size_t start, size_t depth) {
                    if (static_cast<int>(depth) == p + 1) {
                        std::vector<index_t> sorted_face = face;
                        std::sort(sorted_face.begin(), sorted_face.end());
                        index_t face_id = S[p].get_id(sorted_face);
                        stars[p].star_over[face_id].push_back(sp1);
                        return;
                    }
                    for (size_t i = start; i < vertices.size(); ++i) {
                        face[depth] = vertices[i];
                        generate_faces(i + 1, depth + 1);
                    }
                };
                generate_faces(0, 0);
            }
        }
    }

	// ==================== Phase 7.5: Build Boundary Operators for p >= 2 ====================

	for (int p = 2; p <= pmax; ++p) {
		const index_t n_p_simplices = S[p].size();
		const index_t n_p_minus_1_simplices = S[p-1].size();

		// Build boundary operator B[p]: C_p -> C_{p-1}
		std::vector<Eigen::Triplet<double>> triplets;
		triplets.reserve(n_p_simplices * (p + 1));  // Each p-simplex has p+1 faces

		for (index_t sp = 0; sp < n_p_simplices; ++sp) {
			const auto& vertices = S[p].simplex_verts[sp];

			// Generate all (p-1)-faces of this p-simplex
			for (size_t omit_idx = 0; omit_idx < vertices.size(); ++omit_idx) {
				// Create face by omitting vertex at position omit_idx
				std::vector<index_t> face;
				face.reserve(vertices.size() - 1);
				for (size_t j = 0; j < vertices.size(); ++j) {
					if (j != omit_idx) {
						face.push_back(vertices[j]);
					}
				}

				// Face should already be sorted since vertices is sorted
				// Look up this face in S[p-1]
				auto it = S[p-1].id_of.find(face);
				if (it != S[p-1].id_of.end()) {
					index_t face_id = it->second;
					double sign = (omit_idx % 2 == 0) ? 1.0 : -1.0;
					triplets.emplace_back(static_cast<Eigen::Index>(face_id),
										  static_cast<Eigen::Index>(sp),
										  sign);
				} else {
					// This should never happen if complex construction is correct
					Rf_error("Face not found in S[%d] during boundary operator construction", p-1);
				}
			}
		}

		// Assemble the sparse matrix
		L.B[p] = spmat_t(static_cast<Eigen::Index>(n_p_minus_1_simplices),
						 static_cast<Eigen::Index>(n_p_simplices));
		L.B[p].setFromTriplets(triplets.begin(), triplets.end());
		L.B[p].makeCompressed();
	}

    // ==================== Phase 8: Assemble Operators ====================

    g.normalize();
    assemble_operators();
}

void riem_dcx_t::build_knn_riem_dcx(
    const spmat_t& X,
    const vec_t& y,
    index_t k,
    index_t max_p,
    bool use_counting_measure,
    bool directed_knn,
    double density_normalization
) {
    // ==================== Input Validation ====================

    const Eigen::Index n_points = X.rows();
    const Eigen::Index n_features = X.cols();

    // Validate matrix dimensions
    if (n_points <= 0) {
        Rf_error("Feature matrix X must have at least one row (point)");
    }
    if (n_features <= 0) {
        Rf_error("Feature matrix X must have at least one column (feature)");
    }

    // Validate k parameter
    if (k < 2) {
        Rf_error("Parameter k must be at least 2 for meaningful neighborhoods");
    }
    if (k >= static_cast<index_t>(n_points)) {
        Rf_error("Parameter k=%zu cannot exceed n_points-1=%ld",
                 k, (long)(n_points - 1));
    }

    // Validate max_p
    if (max_p < 1) {
        Rf_error("Maximum simplex dimension max_p must be at least 1");
    }
	if (max_p > 10) {
        Rf_warning("Large max_p=%zu may result in excessive memory usage and computation time",
                   max_p);
    }

	// max_p cannot exceed n_points - 1 (simplex with n vertices has dimension n-1)
    if (max_p >= static_cast<index_t>(n_points)) {
        Rf_error("Maximum dimension max_p=%zu cannot exceed n_points-1=%ld",
                 max_p, (long)(n_points - 1));
    }

    // Validate response vector if provided
    if (y.size() > 0 && y.size() != n_points) {
        Rf_error("Response vector y has length %ld but X has %ld points",
                 (long)y.size(), (long)n_points);
    }

    // Validate density_normalization
    if (density_normalization < 0.0) {
        Rf_error("Density normalization must be non-negative, got %f",
                 density_normalization);
    }

    // ==================== Construction ====================

    // Build the nerve complex and Riemannian structure
	// init_dims() will automatically ensure
    // that S.size() == stars.size() == pmax + 1
    build_nerve_from_knn(X, k, max_p, use_counting_measure,
                        directed_knn, density_normalization);

    // ==================== Signal Initialization ====================

    // Store response vector if provided
    if (y.size() == n_points) {
        sig.y = y;

        // Initialize fitted value history (empty initially)
        sig.y_hat_hist.clear();
    } else {
        // No response provided; clear signal state
        sig.y = vec_t::Zero(0);
        sig.y_hat_hist.clear();
    }

    // ==================== Density Initialization ====================

    // Initialize density at vertices to vertex masses (normalized)
    if (!g.M.empty() && g.M[0].rows() == n_points) {
        rho.rho[0] = vec_t::Zero(n_points);
        for (Eigen::Index i = 0; i < n_points; ++i) {
            rho.rho[0][i] = g.M[0].coeff(i, i);
        }

        // Normalize density to sum to n_points
        real_t density_sum = rho.rho[0].sum();
        if (density_sum > 1e-15) {
            rho.rho[0] *= static_cast<real_t>(n_points) / density_sum;
        }
    }

    // Initialize higher-dimensional densities to zero
    for (int p = 1; p <= pmax; ++p) {
        if (rho.rho.size() > static_cast<size_t>(p)) {
            rho.rho[p].setZero();
        }
    }
}
