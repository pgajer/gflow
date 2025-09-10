#include <Eigen/Dense>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include <vector>
#include <set>
#include <cmath>
#include <algorithm>

#include "set_wgraph.hpp"
#include "error_utils.h" // for REPORT_ERROR()

/**
 * Parameterize a circular graph structure using spectral methods.
 *
 * @param use_edge_lengths Whether to use edge lengths as weights in the graph Laplacian
 * @return A struct containing angles and eigenvectors for the circular parameterization
 */
circular_param_result_t set_wgraph_t::parameterize_circular_graph(bool use_edge_lengths) const {

    //Eigen::setNbThreads(12);

    const int n_vertices = static_cast<int>(adjacency_list.size());
    if (n_vertices < 3) {
        REPORT_ERROR("Graph must have at least 3 vertices to parameterize a circle.");
    }

    // Build symmetric weighted adjacency A
    Eigen::SparseMatrix<double> A(n_vertices, n_vertices);
    std::vector<Eigen::Triplet<double>> triples;
    triples.reserve(2 * n_vertices * 4); // heuristic
    for (int i = 0; i < n_vertices; ++i) {
        for (const auto& edge : adjacency_list[i]) {
            const int j = edge.vertex;
            const double w = use_edge_lengths ? edge.weight : 1.0;
            if (i < j) { // add each undirected edge once, both directions
                triples.emplace_back(i, j, w);
                triples.emplace_back(j, i, w);
            }
        }
    }
    A.setFromTriplets(triples.begin(), triples.end());

    // Degree matrix D and Laplacian L = D - A
    Eigen::SparseMatrix<double> D(n_vertices, n_vertices);
    for (int k = 0; k < A.outerSize(); ++k) {
        double sum = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) sum += it.value();
        if (sum > 0.0) D.insert(k, k) = sum;
    }
    Eigen::SparseMatrix<double> L = D - A;

    // Small diagonal regularization for numerical stability
    for (int i = 0; i < n_vertices; ++i) L.coeffRef(i, i) += 1e-8;

    // Spectra operator
    Spectra::SparseSymMatProd<double> op(L);

    // We want up to 6 smallest eigenpairs (beyond the trivial), but are limited by n.
    // Spectra requires: nev < ncv <= n, and typically ncv >= nev+1.
    const int want = std::min(6, std::max(1, n_vertices - 1)); // at most n-1 nontrivial
    int nev = want;
    int ncv = std::min(std::max(2 * nev + 1, nev + 1), n_vertices); // generous subspace, capped at n
    if (nev >= ncv) nev = ncv - 1;
    if (nev < 1)    nev = 1;

    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, nev, ncv);
    eigs.init();
    const int maxit = 1000;
    const double tol = 1e-16;
    eigs.compute(Spectra::SortRule::SmallestAlge, maxit, tol);

    if (eigs.info() != Spectra::CompInfo::Successful) {
        REPORT_ERROR("Eigenvalue estimation with Spectra failed.");
    }

    // Retrieve what was actually computed
    const Eigen::VectorXd eigenvalues = eigs.eigenvalues();
    const Eigen::MatrixXd evectors   = eigs.eigenvectors();
    const int m = static_cast<int>(eigenvalues.size());  // number of columns in evectors

    // Safe debug prints
    for (int i = 0; i < m; ++i) {
        Rprintf("Eigenvalue %d: %.8f\n", i + 1, eigenvalues(i));
    }

    // Need at least 3 eigenvectors (0th ~ constant, then 2nd & 3rd for embedding)
    // We access columns 1 and 2 (0-based): ensure they exist.
    if (m < 3) {
        REPORT_ERROR("Not enough eigenpairs to compute circular coordinates (need at least 3).");
    }

    // Extract available eigenvectors—only index what exists
    const Eigen::VectorXd eig_vec2_e = evectors.col(1);
    const Eigen::VectorXd eig_vec3_e = evectors.col(2);
    Eigen::VectorXd eig_vec4_e, eig_vec5_e, eig_vec6_e;
    const bool has4 = (m > 3);
    const bool has5 = (m > 4);
    const bool has6 = (m > 5);
    if (has4) eig_vec4_e = evectors.col(3);
    if (has5) eig_vec5_e = evectors.col(4);
    if (has6) eig_vec6_e = evectors.col(5);

    // Compute angles from (eig_vec2, eig_vec3)
    std::vector<double> angles(n_vertices);
    Eigen::VectorXd v2 = eig_vec2_e;
    Eigen::VectorXd v3 = eig_vec3_e;
    const double n2 = v2.norm();
    const double n3 = v3.norm();
    if (n2 == 0.0 || n3 == 0.0) {
        REPORT_ERROR("Encountered zero-norm eigenvector; cannot compute angles.");
    }
    v2 /= n2; v3 /= n3;

    for (int i = 0; i < n_vertices; ++i) {
        double ang = std::atan2(v3(i), v2(i));
        if (ang < 0.0) ang += 2.0 * M_PI;
        angles[i] = ang;
    }

    // Marshal to result struct
    circular_param_result_t out;
    out.angles = std::move(angles);

    out.eig_vec2.assign(v2.data(), v2.data() + v2.size());
    out.eig_vec3.assign(v3.data(), v3.data() + v3.size());
    if (has4) out.eig_vec4.assign(eig_vec4_e.data(), eig_vec4_e.data() + eig_vec4_e.size());
    if (has5) out.eig_vec5.assign(eig_vec5_e.data(), eig_vec5_e.data() + eig_vec5_e.size());
    if (has6) out.eig_vec6.assign(eig_vec6_e.data(), eig_vec6_e.data() + eig_vec6_e.size());
    // If not present, those vectors remain empty; the R wrapper returns NULL for them.

    return out;
}


/**
 * Alternative version that allows specifying a reference vertex to set as angle 0
 *
 * @param reference_vertex The vertex to use as reference (angle 0)
 * @param use_edge_lengths Whether to use edge lengths as weights in the graph Laplacian
 * @return A struct containing angles and eigenvectors for the circular parameterization
 */
circular_param_result_t set_wgraph_t::parameterize_circular_graph_with_reference(
	size_t reference_vertex,
	bool use_edge_lengths) const {

	// Get the results using spectral method
	circular_param_result_t result = parameterize_circular_graph(use_edge_lengths);

	if (reference_vertex < adjacency_list.size()) {
		// Get the angle of the reference vertex
		double ref_angle = result.angles[reference_vertex];

		// Shift all angles so that reference_vertex has angle 0
		for (size_t i = 0; i < result.angles.size(); i++) {
			result.angles[i] = fmod(result.angles[i] - ref_angle + 2 * M_PI, 2 * M_PI);
		}
	}

	return result;
}
