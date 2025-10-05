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
    double density_normalization
	) {

#define DEBUG_BUILD_NERVE_FROM_KNN 0

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
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 1 ...\n");
	#endif

    // Convert sparse Eigen matrix to R matrix for ANN library
    SEXP s_X = PROTECT(Rf_allocMatrix(REALSXP, n_points, n_features));
    double* p_X = REAL(s_X);

    for (Eigen::Index j = 0; j < n_features; ++j) {
        for (Eigen::Index i = 0; i < n_points; ++i) {
            p_X[i + j * n_points] = X.coeff(i, j);
        }
    }

    SEXP s_k = PROTECT(Rf_ScalarInteger(k));

	#if 0
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
	for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
		// neighbor_sets[i].insert(i); // closed neighborhood includes self - this is unnecessary as ANN includes self in the set of kNN's
		for (index_t j : knn_neighbors[i]) {
			neighbor_sets[i].insert(j);
		}
	}
    #endif

	// ---------
	struct iknn_vertex_t {
		size_t index; // index of the nearest neighbor
		size_t isize;  // the number of elements in the intersection of the set, N(x), of kNN of the vertext and the set, N(x_j), of kNN of the neighbor x_j
		double dist;  // distance between the vertex and the neighbor computed as the minimumum over k of d(x,x_k) + d(x_k,x_j), where x is the vertex, x_j is the j-th NN of x, x_k runs over elements of N(x) \cap N(x_j) and d(x,y is the distance between x and y as computed by kNN library ANN
	};

	SEXP knn_res = PROTECT(S_kNN(s_X, s_k));
	int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));
    UNPROTECT(3); // knn_res, s_k, s_X

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);

	std::vector<std::vector<iknn_vertex_t>> adjacency_list(n_points);
	std::vector<std::unordered_set<index_t>> neighbor_sets(n_points);

    // Perform k-NN search for each point
    size_t n_points_minus_one = n_points - 1;
    std::vector<int> intersection;
    for (size_t pt_i = 0; pt_i < n_points_minus_one; pt_i++) {
        // Copying indices of kNN of the pt_i point to nn_i
        for (size_t j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + n_points * j];
            sorted_nn_i[j] = nn_i[j];
			neighbor_sets[pt_i].insert(nn_i[j]);
        }
        std::sort(sorted_nn_i.begin(), sorted_nn_i.end()); // Ensure sorted for set intersection

        for (size_t pt_j = pt_i + 1; pt_j < static_cast<index_t>(n_points); pt_j++) {
            // Copying indices of kNN of the pt_j point to nn_j
            for (size_t j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + n_points * j];
                sorted_nn_j[j] = nn_j[j];
            }
            std::sort(sorted_nn_j.begin(), sorted_nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(), sorted_nn_j.begin(), sorted_nn_j.end(), std::back_inserter(intersection));

            size_t common_count = intersection.size();

            if (common_count > 0) {
                // Computing the minimum of d(x,x_k) + d(x_k,x_j)
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    size_t idx_i = std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin(); // The std::find function returns an iterator, not an index. We need to subtract the beginning iterator to get the index.
                    size_t idx_j = std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = distances[pt_i + n_points * idx_i];
                    double dist_j_k = distances[pt_j + n_points * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);

                }

                // Add edge from pt_i to pt_j and from pt_j to pt_i
                adjacency_list[pt_i].emplace_back(iknn_vertex_t{pt_j, common_count, min_dist});
                adjacency_list[pt_j].emplace_back(iknn_vertex_t{pt_i, common_count, min_dist});
            }
        }
    }

    // ==================== Phase 2: Density Computation ====================
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 2 ...\n");
    #endif

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
            //d_k_values[i] = knn_distances[i][k - 1]; // last neighbor
			d_k_values[i] = adjacency_list[i][k - 1].dist;
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
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 3 ...\n");
	#endif

	// Initialize with ONLY vertices (dimension 0)
	std::vector<index_t> n_by_dim = {static_cast<index_t>(n_points)};
	init_dims(0, n_by_dim);  // <-- Initialize pmax = 0, not max_p

	// ==================== Phase 4: Build 0-Simplices (Vertices) ====================
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 4 ...\n");
	#endif

	S[0].simplex_verts.resize(n_points);
	for (index_t i = 0; i < static_cast<index_t>(n_points); ++i) {
		S[0].simplex_verts[i] = {i};
		S[0].id_of[{i}] = i;
	}

	// Initialize vertex masses: m_i = μ(N̂_k(x_i))
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
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Starting edge construction...\n");
    #endif

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

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Building edge list...\n");
	#endif
	// Build edges: edge (i,j) exists iff N̂_k(x_i) ∩ N̂_k(x_j) ≠ ∅

	// neighbor_sets: std::vector<std::unordered_set<index_t>>   // CLOSED: contains i
	// compute_intersection_measure(const std::unordered_set<index_t>& A,
	//                              const std::unordered_set<index_t>& B) -> double
	// n_points: size_t (or index_t-castable)
	// eps threshold:
	const double eps_w = 1e-15;

	const index_t n = static_cast<index_t>(n_points);

	// Build OPEN inverse index: Pre(ell) = { i : ell in N_k(i) } (self excluded)
	std::vector<std::vector<index_t>> inv_neighbors;
	inv_neighbors.assign(n, {});
	for (index_t i = 0; i < n; ++i) {
		inv_neighbors[i].reserve(neighbor_sets[i].size()); // rough
	}
	for (index_t i = 0; i < n; ++i) {
		for (index_t neighbor : neighbor_sets[i]) {
			if (neighbor != i) inv_neighbors[neighbor].push_back(i);
		}
	}

	// Edge candidate generation via witnesses (closed on the forward side)
	auto pack_pair = [](index_t a, index_t b) -> uint64_t {
		if (a > b) std::swap(a, b);
		return (uint64_t(a) << 32) ^ uint64_t(b);
	};

	std::unordered_set<uint64_t> seen;
	seen.reserve(static_cast<size_t>(n) * 3u * (inv_neighbors.empty() ? 1u : inv_neighbors[0].capacity()));

	std::vector<std::array<index_t,2>> edge_list;
	std::vector<double> edge_weights_vec;
	// Rough pre-reserve: ~ 2k edges per vertex (tune as you like)
	edge_list.reserve(static_cast<size_t>(n) * 2u * (inv_neighbors.empty() ? 1u : inv_neighbors[0].capacity()));
	edge_weights_vec.reserve(edge_list.capacity());

	for (index_t i = 0; i < n; ++i) {
		// iterate witnesses in CLOSED forward set (ell may equal i; that’s fine for closed-cover nerve)
		for (index_t ell : neighbor_sets[i]) {
			// pull all j that also have ell as a neighbor (OPEN inverse index)
			const auto& pre = inv_neighbors[ell];
			for (index_t j : pre) {
				if (j == i) continue;

				const uint64_t key = pack_pair(i, j);
				if (!seen.insert(key).second) continue; // already processed this unordered pair

				// tiny constant-factor win: pass smaller closed set first
				const bool i_smaller = neighbor_sets[i].size() <= neighbor_sets[j].size();
				const auto& A = i_smaller ? neighbor_sets[i] : neighbor_sets[j];
				const auto& B = i_smaller ? neighbor_sets[j] : neighbor_sets[i];
				const double w = compute_intersection_measure(A, B);

				if (w > eps_w) {
					const index_t a = std::min(i, j), b = std::max(i, j);
					edge_list.push_back({a, b});
					edge_weights_vec.push_back(w);
				}
			}
		}
	}

	const index_t n_edges = edge_list.size();

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Found %zu edges\n", n_edges);
	#endif

	if (n_edges == 0) {
		Rf_error("No edges found in k-NN graph");
	}

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Calling extend_by_one_dim(%zu)\n", n_edges);
	#endif
	extend_by_one_dim(n_edges);

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Populating edge table...\n");
	#endif
	// Populate edge table
	S[1].simplex_verts.resize(n_edges);
	for (index_t e = 0; e < n_edges; ++e) {
		std::vector<index_t> edge_verts = {edge_list[e][0], edge_list[e][1]};
		S[1].simplex_verts[e] = edge_verts;
		S[1].id_of[edge_verts] = e;
	}

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Building edge metric...\n");
	#endif
	// Build edge masses from intersection measures
	g.M[1] = spmat_t(n_edges, n_edges);
	g.M[1].reserve(Eigen::VectorXi::Constant(n_edges, 1));
	for (index_t e = 0; e < n_edges; ++e) {
		g.M[1].insert(e, e) = std::max(edge_weights_vec[e], 1e-15);
	}
	g.M[1].makeCompressed();

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Building boundary operator...\n");
	#endif
	// Build boundary operator B[1]: edges → vertices
	build_incidence_from_edges();

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Storing conductances...\n");
	#endif
	// Store conductances for potential L0_sym construction
	L.c1 = vec_t(n_edges);
	for (index_t e = 0; e < n_edges; ++e) {
		L.c1[e] = edge_weights_vec[e];
	}

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 5: Complete\n");
	#endif

    // ==================== Phase 6: Build Higher-Dimensional Simplices ====================
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 6 max_p: %d\n", (int)max_p);
    #endif

	if (max_p >= 2) {
		// Build simplices dimension by dimension
		for (int p = 2; p <= static_cast<int>(max_p); ++p) {
			// Use unordered_map to automatically deduplicate and aggregate weights
			std::unordered_map<std::vector<index_t>, double, vec_hash_t> simplex_weight_map;

			// Consider all pairs of (p-1)-simplices
			const index_t n_prev = S[p-1].size();
			for (index_t i = 0; i < n_prev; ++i) {
				for (index_t j = i + 1; j < n_prev; ++j) {
					const auto& s1 = S[p-1].simplex_verts[i];
					const auto& s2 = S[p-1].simplex_verts[j];

					// Compute intersection size
					std::unordered_set<index_t> intersection;
					for (index_t v : s1) {
						if (std::find(s2.begin(), s2.end(), v) != s2.end()) {
							intersection.insert(v);
						}
					}

					// If they share exactly p vertices, merge to form p-simplex
					if (intersection.size() == static_cast<size_t>(p-1)) {
						// Merge the two simplices (union of vertices)
						std::vector<index_t> merged;
						merged.reserve(p + 1);
						std::set_union(s1.begin(), s1.end(),
									   s2.begin(), s2.end(),
									   std::back_inserter(merged));

						// Verify we got a p-simplex (p+1 vertices)
						if (merged.size() != static_cast<size_t>(p + 1)) {
							continue;
						}

						// Compute measure of the p-simplex (intersection of all neighborhoods)
						std::unordered_set<index_t> common_neighbors = neighbor_sets[merged[0]];
						for (size_t v_idx = 1; v_idx < merged.size(); ++v_idx) {
							std::unordered_set<index_t> temp;
							for (index_t v : common_neighbors) {
								if (neighbor_sets[merged[v_idx]].find(v) !=
									neighbor_sets[merged[v_idx]].end()) {
									temp.insert(v);
								}
							}
							common_neighbors = std::move(temp);
						}

						double measure = 0.0;
						for (index_t v : common_neighbors) {
							measure += vertex_weights[v];
						}

						if (measure > 1e-15) {
							// Aggregate weights: use max to handle duplicates
							auto it = simplex_weight_map.find(merged);
							if (it != simplex_weight_map.end()) {
								// Simplex already discovered via another pair
								it->second = std::max(it->second, measure);
							} else {
								// First time discovering this simplex
								simplex_weight_map[merged] = measure;
							}
						}
					}
				}
			}

			#if DEBUG_BUILD_NERVE_FROM_KNN
			Rprintf("Phase 6: Dimension p=%d, found %zu candidate simplices before deduplication\n",
					p, simplex_weight_map.size());
			#endif

			if (simplex_weight_map.empty()) {
				#if DEBUG_BUILD_NERVE_FROM_KNN
				Rprintf("Phase 6: No %d-simplices found, stopping at dimension %d\n",
						p, p - 1);
				#endif
				pmax = p - 1;
				break;
			}

			// Extract deduplicated simplices and weights
			std::vector<std::vector<index_t>> p_simplices;
			std::vector<double> p_weights;
			p_simplices.reserve(simplex_weight_map.size());
			p_weights.reserve(simplex_weight_map.size());

			for (const auto& [simplex, weight] : simplex_weight_map) {
				p_simplices.push_back(simplex);
				p_weights.push_back(weight);
			}

			// Extend complex to dimension p
			extend_by_one_dim(static_cast<index_t>(p_simplices.size()));

			// Add p-simplices to the complex
			S[p].simplex_verts = p_simplices;
			for (index_t s = 0; s < static_cast<index_t>(p_simplices.size()); ++s) {
				S[p].id_of[p_simplices[s]] = s;
			}

			// Build metric matrix for p-simplices
			g.M[p] = spmat_t(p_simplices.size(), p_simplices.size());
			g.M[p].reserve(Eigen::VectorXi::Constant(p_simplices.size(), 1));
			for (index_t s = 0; s < static_cast<index_t>(p_simplices.size()); ++s) {
				g.M[p].insert(s, s) = std::max(p_weights[s], 1e-15);
			}
			g.M[p].makeCompressed();
		}
	}

    // ==================== Phase 7: Build Star Tables ====================
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 7 ...\n");
    #endif

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
	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 7.5 ...\n");
	#endif

	for (int p = 2; p <= pmax; ++p) {
		const index_t n_p_simplices = S[p].size();
		const index_t n_p_minus_1_simplices = S[p-1].size();

		if (n_p_simplices == 0) {
			// No simplices at this dimension, skip
			continue;
		}

		// Build boundary operator B[p]: C_p -> C_{p-1}
		std::vector<Eigen::Triplet<double>> triplets;
		triplets.reserve(n_p_simplices * (p + 1));  // Each p-simplex has p+1 faces

		for (index_t sp = 0; sp < n_p_simplices; ++sp) {
			const auto& vertices = S[p].simplex_verts[sp];

			// Safety check: verify vertex count
			if (vertices.size() != static_cast<size_t>(p + 1)) {
				Rf_error("Simplex %zu in dimension %d has %zu vertices, expected %d",
						 sp, p, vertices.size(), p + 1);
			}

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

					// Safety check: verify face_id is in bounds
					if (face_id >= n_p_minus_1_simplices) {
						Rf_error("Face ID %zu out of bounds (max %zu) for dimension %d",
								 face_id, n_p_minus_1_simplices - 1, p - 1);
					}

					double sign = (omit_idx % 2 == 0) ? 1.0 : -1.0;
					triplets.emplace_back(static_cast<Eigen::Index>(face_id),
										  static_cast<Eigen::Index>(sp),
										  sign);
				} else {
					// Face not found - this indicates a serious error in complex construction
					Rf_error("Face {%s} of simplex %zu in dimension %d not found in S[%d]",
							 [&]() {
								 std::string s;
								 for (size_t i = 0; i < face.size(); ++i) {
									 if (i > 0) s += ",";
									 s += std::to_string(face[i]);
								 }
								 return s.c_str();
							 }(),
							 sp, p, p - 1);
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

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 8: Normalizing metrics...\n");
	#endif
	g.normalize();

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 8: Assembling operators (pmax=%d)...\n", pmax);
	#endif
	assemble_operators();

	#if DEBUG_BUILD_NERVE_FROM_KNN
	Rprintf("Phase 8: Complete\n");
	#endif
}

void riem_dcx_t::build_knn_riem_dcx(
    const spmat_t& X,
    const vec_t& y,
    index_t k,
    index_t max_p,
    bool use_counting_measure,
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
                        density_normalization);

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
/**
   Constructs the reference measure μ that assigns base weights to vertices before any geometric evolution.

   The reference measure provides the starting point for density estimation. Two
   approaches serve different data characteristics: counting measure treats all
   points equally (appropriate for uniform sampling), while distance-based
   measure down-weights points in dense regions (appropriate when sampling
   density varies). The formula $\mu(\{x\}) = d_\ell(x)^{-\alpha}$ with $\alpha
   \in [1,2]$ provides robust local density surrogate without bandwidth
   selection, addressing kernel density estimation instability in
   moderate-to-high dimensions.

 */
void riem_dcx_t::initialize_reference_measure(
    const std::vector<std::vector<index_t>>& knn_neighbors,
    const std::vector<std::vector<double>>& knn_distances,
    bool use_counting_measure,
    double density_normalization
) {
    const size_t n = S[0].size();
    std::vector<double> vertex_weights(n);

    if (use_counting_measure) {
        std::fill(vertex_weights.begin(), vertex_weights.end(), 1.0);
    } else {
        const double epsilon = 1e-10;
        const double alpha = 1.5;
        for (size_t i = 0; i < n; ++i) {
            double d_k = knn_distances[i].back();
            vertex_weights[i] = std::pow(epsilon + d_k, -alpha);
        }
    }

    // Normalize to target sum
    double total = std::accumulate(vertex_weights.begin(),
                                   vertex_weights.end(), 0.0);
    double target = (density_normalization > 0.0) ?
                    density_normalization : static_cast<double>(n);
    if (total > 1e-15) {
        double scale = target / total;
        for (double& w : vertex_weights) w *= scale;
    }

    // Store in member variable for use by other modules
    this->reference_measure = std::move(vertex_weights);
}
