/*!
  Creating Intersetion Weighted kNN graph (IWkNN graph)

  Two points are connected by an edge in this graph if and only if they share at
  least one common point within their respective k-nearest neighbor sets. In
  other words, there's an intersection between their k-nearest neighborhoods.

  \textbf{Characteristics:}

  * \textbf{Reciprocity:} This connection is inherently reciprocal. If point A is within the k-nearest neighbors of point B, then point B will also be within the k-nearest neighbors of point A. This generally results in an undirected graph.
  * \textbf{Sensitivity to k:} The choice of the parameter k significantly influences the structure of this graph. Larger values of k tend to create denser graphs with more connections.
  * \textbf{Focus on Shared Neighbors:} This variant emphasizes points that are "locally similar," meaning they share commonalities within their local neighborhoods.

  \textbf{Potential Use Cases:}

  * \textbf{Clustering and Community Detection:} This type of graph could be valuable in identifying densely connected clusters of points that share neighbors, indicating groups exhibiting similar characteristics.
  * \textbf{Outlier Detection:}  Points with abnormally few connections in this graph might indicate potential outliers, as they don't consistently appear in the k-nearest neighborhoods of other points.
  * \textbf{Collaborative Filtering:} This construct might be applicable in recommendation systems where you want to emphasize items that have been "liked" by users with overlapping tastes.

  \textbf{Caveats:}

  * \textbf{Computational Cost:} Constructing this graph can be computationally more demanding than a standard k-NN graph because you need to compare k-nearest neighbor sets.
  * \textbf{Data Density:}  In very sparse datasets, this construct might lead to highly disconnected graphs.

  \textbf{Relationship with Other Graphs:}

  It's interesting to note that this graph would likely be denser than a Mutual k-Nearest Neighbor Graph (Mk-NN).  In an Mk-NN, two points need to be *each other's* k-nearest neighbors, while in this variant they just need to share at least one common neighbor in their sets.

*/

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <ANN/ANN.h>  // ANN library header
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm> // Required for std::set_intersection and std::sort
#include <memory>
#include <utility>
#include <limits>

#include "msr2.h"
#include "msr2_I_kNN_graph.h"

extern "C" {
    SEXP S_IW_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion, SEXP Rprune_version);
    SEXP S_IW_d1ecdf_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho);

    SEXP S_IWD_kNN_graph(SEXP RX, SEXP Rk);

    SEXP S_I_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_I_d1_kNN_graph(SEXP RX, SEXP Rk);
    SEXP S_I_d1ecdf_kNN_graph(SEXP RX, SEXP Rk);
}

std::unique_ptr<std::vector<double>> ecdf(const std::vector<double> &x);
std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<int>>>& adj_vect);
std::unique_ptr<std::vector<int>> union_find(std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& adj_vect);
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> prune_long_edges(const std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>>& wgraph, int version);
std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph);

/*!
  Creates an intersection weighted kNN graph

  This function creates a weighted graph that connects two points of X with an
  edge if the sets of kNN of these points intersect. Thus, the graph is the
  1-skeleton of the nerve simplicial complex of the kNN covering of X. The weight
  of an edge is the number of points in common between two sets of kNN's.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a vector of vector of pairs, adj_kNN, such that adj_kNN[i] is
  a vector of pairs: (i1, w1), (i2, w2), ... , (iN, wN), with each pair being
  the index of the neighbor and the corresponding edge weight.
*/
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX1 = nrX - 1;
    int common_count = 0;
    std::vector<int> intersection;
    for (int pt_i = 0; pt_i < nrX1; pt_i++) {
        for (int j = 0; j < k; j++)  // Copying indices of kNN of the pt_i point to nn_i
            nn_i[j] = indices[pt_i + nrX * j];

        std::sort(nn_i.begin(), nn_i.end()); // Ensure sorted for set intersection

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {

            for (int j = 0; j < k; j++) // Copying indices of kNN of the pt_j point to nn_j
                nn_j[j] = indices[pt_j +  nrX * j];

            std::sort(nn_j.begin(), nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(nn_i.begin(), nn_i.end(), nn_j.begin(), nn_j.end(), std::back_inserter(intersection));

            common_count = intersection.size();
            if (common_count > 0) {
                // Add edge from pt_i to pt_j and from pt_j to pt_i with the weight
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

//
// Version 2 utilizing std::unordered_set<int> set_i/set_j
//
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_v2(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX1 = nrX - 1;
    int common_count = 0;

    std::unordered_set<int> set_i;
    std::unordered_set<int> set_j;

    for (int pt_i = 0; pt_i < nrX1; pt_i++) {
        set_i.clear();
        // Populate set_i for pt_i
        for (int j = 0; j < k; j++) {
            set_i.insert(indices[pt_i + nrX * j]);
        }

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            set_j.clear();
            for (int j = 0; j < k; j++) {
                set_j.insert(indices[pt_j + nrX * j]);
            }

            common_count = 0;
            // If set_i is smaller, iterate over it for efficiency
            if (set_i.size() < set_j.size()) {
                for (const int& elem : set_i) {
                    if (set_j.count(elem) > 0) {
                        common_count++;
                    }
                }
            } else {
                for (const int& elem : set_j) {
                    if (set_i.count(elem) > 0) {
                        common_count++;
                    }
                }
            }

            if (common_count > 0) {
                // Add edge from pt_i to pt_j and from pt_j to pt_i with the weight
                // Assuming res is a properly initialized structure to hold these values
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

// ------------------------------------------------------------------------------------------
//
// IWD_kNN_graph
//
// ------------------------------------------------------------------------------------------

/**
 * Transforms an Intersection-Weighted-Distance k-Nearest Neighbor graph to an Intersection-Weighted k-Nearest Neighbor graph.
 *
 * This function converts a graph representation where each edge contains a distance measurement alongside the intersection size
 * to a representation that only considers the intersection size, discarding the distance information.
 *
 * @param graph The original graph represented as a vector of vectors of IkNN_vertex_t, where each IkNN_vertex_t contains
 *              the nearest neighbor index, intersection size, and distance.
 * @return A unique pointer to a vector of vectors of pairs, where each pair consists of a nearest neighbor index and an intersection size.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IWD_to_IW_kNN_graph(const std::vector<std::vector<IkNN_vertex_t>>& graph) {
    auto iw_graph = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(graph.size());

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto& vertex : graph[i]) {
            (*iw_graph)[i].emplace_back(vertex.nn_index, vertex.I_index);
        }
    }

    return iw_graph;
}


/**
 * @brief Constructs an Intersection-Weighted-Distance k-Nearest Neighbors (kNN) graph.
 *
 * This function builds an Intersection-Weighted-Distance kNN graph from a given data matrix.
 * For each point, it finds its k nearest neighbors, then constructs edges between points
 * based on the number of common nearest neighbors and computes the distance as the
 * minimum sum of distances through common neighbors.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @return A unique pointer to a 2D vector of IkNN_vertex_t structures representing the kNN graph.
 *
 * Each IkNN_vertex_t structure contains:
 * - `nn_index`: Index of the nearest neighbor.
 * - `I_index`: Number of common elements in the kNN sets.
 * - `dist`: Distance metric computed as the minimum sum of distances through common neighbors.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 */
std::unique_ptr<std::vector<std::vector<IkNN_vertex_t>>> IWD_kNN_graph(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *distances = REAL(VECTOR_ELT(knn_res, 1));

    std::vector<int> nn_i(k);
    std::vector<int> nn_j(k);
    std::vector<int> sorted_nn_i(k);
    std::vector<int> sorted_nn_j(k);

    auto res = std::make_unique<std::vector<std::vector<IkNN_vertex_t>>>(nrX);

    // Perform k-NN search for each point
    int nrX1 = nrX - 1;
    std::vector<int> intersection;
    for (int pt_i = 0; pt_i < nrX1; pt_i++) {
        // Copying indices of kNN of the pt_i point to nn_i
        for (int j = 0; j < k; j++) {
            nn_i[j] = indices[pt_i + nrX * j];
            sorted_nn_i[j] = nn_i[j];
        }

        std::sort(sorted_nn_i.begin(), sorted_nn_i.end()); // Ensure sorted for set intersection

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            // Copying indices of kNN of the pt_j point to nn_j
            for (int j = 0; j < k; j++) {
                nn_j[j] = indices[pt_j + nrX * j];
                sorted_nn_j[j] = nn_j[j];
            }

            std::sort(sorted_nn_j.begin(), sorted_nn_j.end()); // Ensure sorted for set intersection

            intersection.clear(); // Clear the intersection vector before reusing it
            std::set_intersection(sorted_nn_i.begin(), sorted_nn_i.end(), sorted_nn_j.begin(), sorted_nn_j.end(), std::back_inserter(intersection));

            int common_count = intersection.size();
            if (common_count > 0) {
                // Computing the minimum of d(x,x_k) + d(x_k,x_j)
                double min_dist = std::numeric_limits<double>::max();
                for (int x_k : intersection) {
                    auto idx_i = std::find(nn_i.begin(), nn_i.end(), x_k) - nn_i.begin(); // The std::find function returns an iterator, not an index. We need to subtract the beginning iterator to get the index.
                    auto idx_j = std::find(nn_j.begin(), nn_j.end(), x_k) - nn_j.begin();
                    double dist_i_k = distances[pt_i + nrX * idx_i];
                    double dist_j_k = distances[pt_j + nrX * idx_j];
                    min_dist = std::min(min_dist, dist_i_k + dist_j_k);
                }

                // Add edge from pt_i to pt_j and from pt_j to pt_i
                (*res)[pt_i].emplace_back(IkNN_vertex_t{pt_j, common_count, min_dist});
                (*res)[pt_j].emplace_back(IkNN_vertex_t{pt_i, common_count, min_dist});
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/**
 * S_IWD_kNN_graph
 *
 * This function computes a weighted k-nearest neighbors graph based on the input data matrix.
 * It returns a list containing the adjacency list of the graph, intersection sizes, distances,
 * connected components, pruned adjacency list, and sizes of cycles within the pruned graph.
 *
 * @param RX A numeric matrix (in column-major order) representing the data points.
 * @param Rk An integer vector of length 1 representing the number of nearest neighbors (k).
 *
 * @return A list with the following elements:
 *   - adj_list: A list of integer vectors where each vector represents the indices of the k-nearest neighbors for each point.
 *   - isize_list: A list of integer vectors where each vector represents the intersection sizes for the k-nearest neighbors.
 *   - dist_list: A list of numeric vectors where each vector represents the distances to the k-nearest neighbors.
 *   - conn_comps: An integer vector representing the connected components of the graph.
 *   - pruned_adj_list: A list of integer vectors where each vector represents the indices of the k-nearest neighbors for each point in the pruned graph.
 *
 * @note The function uses the base R C++ API and does not rely on Rcpp.
 *
 * @error Throws an error if the input data matrix cannot be coerced to a numeric matrix.
 * @error Throws an error if the k-nearest neighbors graph function returns a null pointer.
 *
 * Example:
 * \code
 * SEXP result = S_IWD_kNN_graph(RX, Rk);
 * \endcode
 */
SEXP S_IWD_kNN_graph(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Creating a kNN graph
    std::unique_ptr<std::vector<std::vector<IkNN_vertex_t>>> IWD_kNN_graph_res = IWD_kNN_graph(RX, Rk);
    if (!IWD_kNN_graph_res) Rf_error("IWD_kNN_graph function returned an invalid result (null pointer).");

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto dist_vect = std::make_unique<std::vector<std::vector<double>>>(nrX);
    int IWD_kNN_graph_n_vertices = IWD_kNN_graph_res->size();
    for (int i = 0; i < IWD_kNN_graph_n_vertices; i++) {
        for (auto nn_vertex : IWD_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_vertex.nn_index);
            (*intersection_size_vect)[i].push_back(nn_vertex.I_index);
            (*dist_vect)[i].push_back(nn_vertex.dist);
        }
    }

    // Preparing adjacency list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Preparing intersection size list
    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);
        for (auto size : intersection_size_vect->at(i))
            *W++ = size;
        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    // Preparing distance list
    SEXP dist_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RD = PROTECT(allocVector(REALSXP, dist_vect->at(i).size())); nprot++;
        double* D = REAL(RD);
        for (auto dist : dist_vect->at(i))
            *D++ = dist;
        SET_VECTOR_ELT(dist_list, i, RD);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the IWD graph by deriving a IW kNN graph and then pruning that graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph = IWD_to_IW_kNN_graph(*IWD_kNN_graph_res);

    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph, 1);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_vertex : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_vertex.first);
        }
    }

    // Preparing pruned adjacency list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    // Preparing the result list
    SEXP res = PROTECT(allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, dist_list);
    SET_VECTOR_ELT(res, 3, conn_comps);
    SET_VECTOR_ELT(res, 4, pruned_adj_list);

    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("dist_list"));
    SET_STRING_ELT(names, 3, mkChar("conn_comps"));
    SET_STRING_ELT(names, 4, mkChar("pruned_adj_list"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}



/*!
  Creates an intersection kNN graph without weights

  This function creates an undirected unweighted graph that connects two points
  of X with an edge if the sets of kNN of these points intersect. Thus, the
  graph is the 1-skeleton of the Cech simplicial complex of the kNN covering of
  X.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a pointer to a vector of vectors, adj_kNN, such that
  adj_kNN[i] is the set of neighbors of the i-th vertex.
*/
std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_I_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    SEXP knn_res = PROTECT(S_kNN(RX, Rk)); nprot++;
    int *nn_indices = INTEGER(VECTOR_ELT(knn_res, 0));

#if DEBUG_I_kNN_graph
    // printing the content of nn_indices
    Rprintf("\nnn_indices\n");
    for (int point_i = 0; point_i < nrX; point_i++) {
        Rprintf("i: %d: ", point_i + 1);
        for (int j = 0; j < k; j++)
            Rprintf(" %d", nn_indices[point_i + nrX * j] + 1);
        Rprintf("\n");
    }
#endif

    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX);

#if DEBUG_I_kNN_graph
    Rprintf("In I_kNN_graph before the main loop\n");
#endif

    int n_points_minus_one = nrX - 1;
    std::set<int> knn_set_i;
    std::set<int> knn_set_j;
    std::vector<int> intersection; // For storing the intersection result
    for (int point_i = 0; point_i < n_points_minus_one; point_i++) {
        knn_set_i.clear();
        for (int j = 0; j < k; j++)
            knn_set_i.insert(nn_indices[point_i + nrX * j]);

#if DEBUG_I_kNN_graph
        Rprintf("\n---> pt_i: %d\n", point_i);
        Rprintf("knn_set_i: ");
        for (auto el : knn_set_i)
            Rprintf(" %d,", el);
        Rprintf("\n");
#endif

        for (int point_j = point_i + 1; point_j < nrX; point_j++) {
            knn_set_j.clear();
            for (int j = 0; j < k; j++)  // j = 0 is the point whose NN we are computing, so we skip it
                knn_set_j.insert(nn_indices[point_j + nrX * j]);

#if DEBUG_I_kNN_graph
            Rprintf("\npt_j: %d\n", point_j);
            Rprintf("knn_set_j: ");
            for (auto el : knn_set_j)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            intersection.clear();
            std::set_intersection(knn_set_i.begin(), knn_set_i.end(),
                                  knn_set_j.begin(), knn_set_j.end(),
                                  std::back_inserter(intersection));

#if DEBUG_I_kNN_graph
            Rprintf("intersection: ");
            for (auto el : intersection)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            if (!intersection.empty()) {
                (*res)[point_i].push_back(point_j);
                (*res)[point_j].push_back(point_i);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}


/*!
  Creates an intersection kNN graph accounting for the density of the data in
  the construct of kNN's using the distance to the nearest neighbor.

  This function creates an undirected graph that connects two points of X with
  an edge if the sets of kNN of these points intersect. The distance to nearest
  neighbors (NN) is modified by changing the distance, d_j, to the j-th nearest
  neighbor of the given point to be, d_j * d1[j], where d1[j] is the distance to
  the first nearest neighbor for the j-th NN of the given point. An idea behind
  this modification is like this: Suppose d_j = d_{j+1} but the density of X at
  the j-th NN is higher than at the (j+1)-st NN. Since the only distance that
  makes sense in a complex data is the geodesic distance, we need to modify d_j
  and d_{j+1} in a way consistent with what the geodesic distance between the
  query point and the corresponding NN's suppose to be. If X was a conneted
  Riemannian manifold or more generally, a path connected metric space, a
  geodesic between two points, x0 and x1, would be the shortest path, within X,
  that connects x0 with x1. The question is how 'within X' statement needs to be
  modified in the context of a discrete dataset X embedded in some R^d. We posit
  that a geodesic in this case is a path in R^d connecting x0 and x1 that is
  minimizing the energy functional that is the intergal of the inverse of the
  density of X. Such a functional makes the geodesics of X traverse through the
  densest possible regions of X to connect any two points. In this function I
  use 1/d1 as a surrogate of the density of X at the given point, thus the
  inverse of it multiplies the distances by d1.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a pointer to a vector of vectors, adj_kNN, such that
  adj_kNN[i] is the set of neighbors of the i-th vertex.
*/
std::unique_ptr<std::vector<std::vector<int>>> I_d1_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_I_d1_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];
    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    // using the double of k
    SEXP Rtwo_k = PROTECT(allocVector(INTSXP, 1)); nprot++;
    int two_k = 2 * k;
    if (two_k >= nrX) {
        two_k = nrX - 1;
    }

    INTEGER(Rtwo_k)[0] = two_k;

    SEXP knn_res = PROTECT(S_kNN(RX, Rtwo_k)); nprot++;
    int *two_k_indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *two_k_distances = REAL(VECTOR_ELT(knn_res, 1));

    // Creating an array of distances to the nearest neighbor
    double *dist_to_NN = two_k_distances + nrX;

    // Indices of kNN's computed with d1 distance corrections
    std::vector<int> k_indices(k * nrX);

    std::vector<std::pair<int,double>> nn_idist(two_k - 1); // ANN includes the query point as the first of the kNN's, we are sorting kNN's excluding the query point
    int nn_j;
    for (int i = 0; i < nrX; i++) {
        nn_idist.clear();  // Reset for the next point
        for (int j = 1; j < two_k; j++) {
            nn_j = two_k_indices[i + nrX * j];
            nn_idist.emplace_back(nn_j, two_k_distances[i + nrX * j] * dist_to_NN[nn_j]);
        }

        // Sorting nn_idist in the ascending order with respect to the second component of each pair
        std::sort(nn_idist.begin(), nn_idist.end(), [&](std::pair<int,double> a, std::pair<int,double> b) {
            return a.second < b.second;
        });

        // Populating k_indices with d1-corrected kNN's indices
        k_indices[i] = i;
        for (int j = 1; j < k; j++)
            k_indices[i + nrX * j] = nn_idist[j - 1].first;
    }
    
    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX);

#if DEBUG_I_d1_kNN_graph
    Rprintf("In I_d1_kNN_graph before the main loop\n");
#endif

    int n_points_minus_one = nrX - 1;
    std::set<int> knn_set_i;
    std::set<int> knn_set_j;
    std::vector<int> intersection; // For storing the intersection result
    for (int point_i = 0; point_i < n_points_minus_one; point_i++) {
        knn_set_i.clear();
        for (int j = 0; j < k; j++) // j = 0 is the point whose NN we are computing, so we skip it
            knn_set_i.insert(k_indices[point_i + nrX * j]);

#if DEBUG_I_d1_kNN_graph
        Rprintf("\n---> pt_i: %d\n", point_i);
        Rprintf("knn_set_i: ");
        for (auto el : knn_set_i)
            Rprintf(" %d,", el);
        Rprintf("\n");
#endif

        for (int point_j = point_i + 1; point_j < nrX; point_j++) {
            knn_set_j.clear();
            for (int j = 0; j < k; j++)  // j = 0 is the point whose NN we are computing, so we skip it
                knn_set_j.insert(k_indices[point_j + nrX * j]);

#if DEBUG_I_d1_kNN_graph
            Rprintf("\npt_j: %d\n", point_j);
            Rprintf("knn_set_j: ");
            for (auto el : knn_set_j)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            intersection.clear();
            std::set_intersection(knn_set_i.begin(), knn_set_i.end(),
                                  knn_set_j.begin(), knn_set_j.end(),
                                  std::back_inserter(intersection));

#if DEBUG_I_d1_kNN_graph
            Rprintf("intersection: ");
            for (auto el : intersection)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            if (!intersection.empty()) {
                (*res)[point_i].push_back(point_j);
                (*res)[point_j].push_back(point_i);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates an intersection kNN graph accounting for the density of the data in
  the construct of kNN's using the distance to the nearest neighbor.

  This function creates an undirected graph that connects two points of X with
  an edge if the sets of kNN of these points intersect. The distance to nearest
  neighbors (NN) is modified by changing the distance, d_j, to the j-th nearest
  neighbor of the given point to be, d_j * ecdfd1_j, where ecdfd1_j is the
  empirical cummulative density function of 1 / d1_j, where d1_j is the distance
  to the first nearest neighbor for the j-th NN of the given point.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a pointer to a vector of vectors, adj_kNN, such that
  adj_kNN[i] is the set of neighbors of the i-th vertex.

 */
std::unique_ptr<std::vector<std::vector<int>>> I_d1ecdf_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_I_d1ecdf_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];
    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X
    // using the double of k
    SEXP Rtwo_k = PROTECT(allocVector(INTSXP, 1)); nprot++;
    int two_k = 2 * k;
    if (two_k >= nrX) {
        two_k = nrX - 1;
    }

    INTEGER(Rtwo_k)[0] = two_k;

    // Finding kNN for K = 2*k
    SEXP knn_res = PROTECT(S_kNN(RX, Rtwo_k)); nprot++;
    int *two_k_indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *two_k_distances = REAL(VECTOR_ELT(knn_res, 1));

    // Creating an array of distances to the nearest neighbor
    double *dist_to_NN = two_k_distances + nrX;

    // Copying the values  1 / dist_to_NN[j] to d1_inv
    std::vector<double> d1_inv(nrX);
    for (int i = 0; i < nrX; i++)
        d1_inv[i] = 1.0 / dist_to_NN[i];

    // Call the C++ ecdf function
    std::unique_ptr<std::vector<double>> ecdf_d1_inv = ecdf(d1_inv);

    // Indices of kNN's computed with d1 distance corrections
    std::vector<int> k_indices(k * nrX);

    std::vector<std::pair<int,double>> nn_idist(two_k - 1); // ANN includes the query point as the first of the kNN's, we are sorting kNN's excluding the query point
    int nn_j;
    for (int i = 0; i < nrX; i++) {
        nn_idist.clear();  // Reset for the next point
        for (int j = 1; j < two_k; j++) {
            nn_j = two_k_indices[i + nrX * j];
            nn_idist.emplace_back(nn_j, two_k_distances[i + nrX * j] / (*ecdf_d1_inv)[nn_j]);
        }

        // Sorting nn_idist in the ascending order with respect to the second component of each pair
        std::sort(nn_idist.begin(), nn_idist.end(), [&](std::pair<int,double> a, std::pair<int,double> b) {
            return a.second < b.second;
        });

        // Populating k_indices with d1-corrected kNN's indices
        k_indices[i] = i;
        for (int j = 1; j < k; j++)
            k_indices[i + nrX * j] = nn_idist[j - 1].first;
    }

    auto res = std::make_unique<std::vector<std::vector<int>>>(nrX);

#if DEBUG_I_d1ecdf_kNN_graph
    Rprintf("In I_d1ecdf_kNN_graph before the main loop\n");
#endif

    int n_points_minus_one = nrX - 1;
    std::set<int> knn_set_i;
    std::set<int> knn_set_j;
    std::vector<int> intersection; // For storing the intersection result
    for (int point_i = 0; point_i < n_points_minus_one; point_i++) {
        knn_set_i.clear();
        for (int j = 0; j < k; j++) // j = 0 is the point whose NN we are computing, so we skip it
            knn_set_i.insert(k_indices[point_i + nrX * j]);

#if DEBUG_I_d1ecdf_kNN_graph
        Rprintf("\n---> pt_i: %d\n", point_i);
        Rprintf("knn_set_i: ");
        for (auto el : knn_set_i)
            Rprintf(" %d,", el);
        Rprintf("\n");
#endif

        for (int point_j = point_i + 1; point_j < nrX; point_j++) {
            knn_set_j.clear();
            for (int j = 0; j < k; j++)  // j = 0 is the point whose NN we are computing, so we skip it
                knn_set_j.insert(k_indices[point_j + nrX * j]);

#if DEBUG_I_d1ecdf_kNN_graph
            Rprintf("\npt_j: %d\n", point_j);
            Rprintf("knn_set_j: ");
            for (auto el : knn_set_j)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            intersection.clear();
            std::set_intersection(knn_set_i.begin(), knn_set_i.end(),
                                  knn_set_j.begin(), knn_set_j.end(),
                                  std::back_inserter(intersection));

#if DEBUG_I_d1ecdf_kNN_graph
            Rprintf("intersection: ");
            for (auto el : intersection)
                Rprintf(" %d,", el);
            Rprintf("\n");
#endif

            if (!intersection.empty()) {
                (*res)[point_i].push_back(point_j);
                (*res)[point_j].push_back(point_i);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}


/*! Creates a weighted intersection kNN graph accounting for the density of the
  data in the construct of kNN's using the empirical cummulative density
  function of the distance to the nearest neighbor.

  \param X A matrix of points.
  \param k The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns a pointer to a vector of vector pairs, adj_kNN, such that
  adj_kNN[i] is the vector of neighbors and the corresponding intersection sizes
  of the i-th vertex.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_d1ecdf_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG_IW_d1ecdf_kNN_graph 0

#if DEBUG_IW_d1ecdf_kNN_graph
    Rprintf("In IW_d1ecdf_kNN_graph\n");
#endif

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];

    int k = INTEGER(Rk)[0];

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X with K = 2 * k
    SEXP Rtwo_k = PROTECT(allocVector(INTSXP, 1)); nprot++;
    int two_k = 2 * k;
    if (two_k >= nrX) {
        two_k = nrX - 1;
    }
    INTEGER(Rtwo_k)[0] = two_k;

    // Finding kNN for K = 2*k
    SEXP knn_res = PROTECT(S_kNN(RX, Rtwo_k)); nprot++;
    int *two_k_indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *two_k_distances = REAL(VECTOR_ELT(knn_res, 1));

    // Creating an array of distances to the nearest neighbor
    double *dist_to_NN = two_k_distances + nrX;

    // Copying the values  1 / dist_to_NN[j] to d1_inv
    std::vector<double> d1_inv(nrX);
    for (int i = 0; i < nrX; i++)
        d1_inv[i] = 1.0 / dist_to_NN[i];

    // Call the C++ ecdf function
    std::unique_ptr<std::vector<double>> ecdf_d1_inv = ecdf(d1_inv);

    // Indices of kNN's computed with d1 distance corrections
    std::vector<int> k_indices(k * nrX);
    std::vector<std::pair<int,double>> nn_idist(two_k); // ANN includes the query point as the first of the kNN's, we are sorting kNN's excluding the query point
    int nn_j;

    for (int i = 0; i < nrX; i++) {

#if DEBUG_IW_d1ecdf_kNN_graph
        Rprintf("\ni: %d\n", i);
#endif

        nn_idist.clear();  // Reset for the next point
        for (int j = 1; j < two_k; j++) {
            nn_j = two_k_indices[i + nrX * j];
            nn_idist.emplace_back(nn_j, two_k_distances[i + nrX * j] / (*ecdf_d1_inv)[nn_j]);
        }
#if DEBUG_IW_d1ecdf_kNN_graph
        Rprintf("nn_idist before sorting\n");
        for (auto p :  nn_idist)
            Rprintf("%d\t%.6f\n", p.first, p.second);
        Rprintf("\n");
#endif
        // Sorting nn_idist in the ascending order with respect to the second component of each pair
        std::sort(nn_idist.begin(), nn_idist.end(),
                  [&](std::pair<int,double> a, std::pair<int,double> b) {
            return a.second < b.second;
        });

#if DEBUG_IW_d1ecdf_kNN_graph
        Rprintf("nn_idist after sorting\n");
        for (auto p :  nn_idist)
            Rprintf("%d\t%.6f\n", p.first, p.second);
        Rprintf("\n");
#endif

        // Populating k_indices with d1-corrected kNN's indices
        k_indices[i] = i;
        for (int j = 1; j < k; j++)
            k_indices[i + nrX * j] = nn_idist[j - 1].first;

#if DEBUG_IW_d1ecdf_kNN_graph
        Rprintf("k_indices[i]\n");
        for (int j = 0; j < k; j++)
            Rprintf("%d\n", k_indices[i + nrX * j]);
        Rprintf("\n");
#endif

    }

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    int common_count = 0;

    std::unordered_set<int> set_i;
    std::unordered_set<int> set_j;

#if DEBUG_IW_d1ecdf_kNN_graph
    Rprintf("In IW_d1ecdf_kNN_graph before pt_i loop\n");
#endif

    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        set_i.clear();
        // Populate set_i for pt_i
        for (int j = 0; j < k; j++)
            set_i.insert(k_indices[pt_i + nrX * j]);

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            set_j.clear();
            for (int j = 0; j < k; j++)
                set_j.insert(k_indices[pt_j + nrX * j]);

            common_count = 0;
            // If set_i is smaller, iterate over it for efficiency
            if (set_i.size() < set_j.size()) {
                for (const int& elem : set_i) {
                    if (set_j.count(elem) > 0) {
                        common_count++;
                    }
                }
            } else {
                for (const int& elem : set_j) {
                    if (set_i.count(elem) > 0) {
                        common_count++;
                    }
                }
            }

            if (common_count > 0) {
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}





/*!
  Creates an Intersection kNN graph (with no weights)

  This function calls an internal C++ function I_kNN_graph() to create an interesection kNN graph.

  \param RX          A matrix of points.
  \param Rk          The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns an adjacency list of point indices adjacent to each vertex.
*/
SEXP S_I_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) { // This is here only in case the R function that calls this one is changed in a way that removes the conversion of X to a matrix
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    // Creating an intersection kNN graph
    std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph_res = I_kNN_graph(RX, Rk);

    // Preparing adj_list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, I_kNN_graph_res->at(i).size())); nprot++; // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : I_kNN_graph_res->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> I_kNN_graph_conn_comp = union_find(I_kNN_graph_res);

    SEXP conn_comps = PROTECT(allocVector(INTSXP, I_kNN_graph_conn_comp->size())); nprot++; // Graph connected components
    std::copy(I_kNN_graph_conn_comp->begin(), I_kNN_graph_conn_comp->end(), INTEGER(conn_comps));

    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, conn_comps);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("conn_comps"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates a d1-density corrected Intersection kNN graph

  This function calls an internal C++ function I_d1_kNN_graph() to create an density corrected interesection kNN graph.

  \param RX          A matrix of points.
  \param Rk          The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns an adjacency list of point indices adjacent to each vertex.
*/
SEXP S_I_d1_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) { // This is here only in case the R function that calls this one is changed in a way that removes the conversion of X to a matrix
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph_res = I_d1_kNN_graph(RX, Rk);

    // Preparing return list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, I_kNN_graph_res->at(i).size())); nprot++; // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : I_kNN_graph_res->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> I_kNN_graph_conn_comp = union_find(I_kNN_graph_res);

    SEXP conn_comps = PROTECT(allocVector(INTSXP, I_kNN_graph_conn_comp->size())); nprot++; // Graph connected components
    std::copy(I_kNN_graph_conn_comp->begin(), I_kNN_graph_conn_comp->end(), INTEGER(conn_comps));

    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, conn_comps);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("conn_comps"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates a d1ecdf-density corrected Intersection kNN graph

  This function calls an internal C++ function I_d1ecdf_kNN_graph() to create an density corrected interesection kNN graph.

  \param RX          A matrix of points.
  \param Rk          The number of nearest neighbors to use for the construction of the kNN graph.

  \return Returns an adjacency list of point indices adjacent to each vertex.
*/
SEXP S_I_d1ecdf_kNN_graph(SEXP RX, SEXP Rk) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) { // This is here only in case the R function that calls this one is changed in a way that removes the conversion of X to a matrix
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    std::unique_ptr<std::vector<std::vector<int>>> I_kNN_graph_res = I_d1ecdf_kNN_graph(RX, Rk);

    // Preparing return list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++; // List with nrX components; each component hold indices of points adjacent to the given point (index of the list)
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, I_kNN_graph_res->at(i).size())); nprot++; // Adjacency list of the i-th vertex
        int* A = INTEGER(RA);

        for (auto neighbor : I_kNN_graph_res->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> I_kNN_graph_conn_comp = union_find(I_kNN_graph_res);

    SEXP conn_comps = PROTECT(allocVector(INTSXP, I_kNN_graph_conn_comp->size())); nprot++; // Graph connected components
    std::copy(I_kNN_graph_conn_comp->begin(), I_kNN_graph_conn_comp->end(), INTEGER(conn_comps));

    SEXP res = PROTECT(allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, conn_comps);

    SEXP names = PROTECT(allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("conn_comps"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

/*!
  Creates a weighted d1ecdf-density corrected Intersection kNN graph.

  This function performs the following steps:
  1. Calls an internal C++ function IW_d1ecdf_kNN_graph() to create a density-corrected weighted intersection kNN graph.
  2. Computes the connected components of the graph using the union_find algorithm.
  3. Prunes the long edges within the graph using the prune_long_edges function, producing a copy of the graph with pruned long edges.

  \param RX A matrix of points.
  \param Rk The number of nearest neighbors to use for the construction of the kNN graph.

  \return A list with three components:
  1. adj_list: An adjacency list of point indices adjacent to each vertex.
  2. conn_comps: A vector of graph connected component indicators.
  3. pruned_adj_list: The pruned graph with 'long' edges removed.
*/
SEXP S_IW_d1ecdf_kNN_graph(SEXP RX, SEXP Rk) {

#define DEBUG__S_IW_d1ecdf_kNN_graph 0

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

#if DEBUG__S_IW_d1ecdf_kNN_graph
    Rprintf("In S_IW_d1ecdf_kNN_graph\n");
#endif

    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_res = IW_d1ecdf_kNN_graph(RX, Rk);

#if DEBUG__S_IW_d1ecdf_kNN_graph
    Rprintf("\nIW_kNN_graph_res:\n");
    for (int i = 0; i < IW_kNN_graph_res->size(); i++) {
        Rprintf("i: %d\n", i + 1);
        for (auto neighbor_pair : (*IW_kNN_graph_res)[i])
            Rprintf("%d\t%d\n", neighbor_pair.first + 1, neighbor_pair.second);
        Rprintf("\n");
    }
#endif

#if DEBUG__S_IW_d1ecdf_kNN_graph
    Rprintf("In S_IW_d1ecdf_kNN_graph after IW_d1ecdf_kNN_graph(RX, Rk)\n");
#endif

    // Preparing I_kNN_graph adj list
    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, IW_kNN_graph_res->at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor_pair : IW_kNN_graph_res->at(i))
            *A++ = neighbor_pair.first + 1;
        SET_VECTOR_ELT(adj_list, i, RA);
    }

    // Identifying connected components of the graph
    std::unique_ptr<std::vector<int>> I_kNN_graph_conn_comp = union_find(IW_kNN_graph_res);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, I_kNN_graph_conn_comp->size())); nprot++;
    std::copy(I_kNN_graph_conn_comp->begin(), I_kNN_graph_conn_comp->end(), INTEGER(conn_comps));

#if DEBUG__S_IW_d1ecdf_kNN_graph
    Rprintf("In S_IW_d1ecdf_kNN_graph before prune_long_edges(IW_kNN_graph_res)\n");
#endif

    //
    // Pruning the graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph_res, 1);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_pair : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_pair.first);
        }
    }

    // Preparing pruned_I_kNN_graph adj list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    // Preparing the return list
    SEXP res = PROTECT(allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, conn_comps);
    SET_VECTOR_ELT(res, 2, pruned_adj_list);

    SEXP names = PROTECT(allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("conn_comps"));
    SET_STRING_ELT(names, 2, mkChar("pruned_adj_list"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}


/*!
  Creates an intersection weighted kNN graph

  This function calls an internal C++ function kNN_graph() to create a kNN graph.

  \param RX               A matrix of points.
  \param Rk               The number of nearest neighbors to use for the construction of the kNN graph.
  \param Rversion         A version of the IW_kNN_graph() function.
  \param Rprune_version   A version of the prune_long_edges() function to use.

  \return A list with four components:
  1. adj_list: An adjacency list of point indices adjacent to each vertex.
  2. isize_list: A list of the sizes of the intersections of kNN sets.
  3. conn_comps: A vector of graph connected component indicators.
  4. pruned_adj_list: The pruned graph with 'long' edges removed.

*/
SEXP S_IW_kNN_graph(SEXP RX, SEXP Rk, SEXP Rversion, SEXP Rprune_version) {

#define DEBUG__S_IW_kNN_graph 0
#define DEBUG__S_IW_kNN_graph_L1 0

    int version = INTEGER(Rversion)[0];
    if (version < 1 || version > 2) {
        Rf_error("Invalid 'Rversion' value. Must be 1 or 2");
    }

    int prune_version = INTEGER(Rprune_version)[0];
    if (prune_version < 1 || prune_version > 2) {
        Rf_error("Invalid 'Rprune_version' value. Must be 1 or 2");
    }

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("\nIn S_IW_kNN_graph()\n");
    Rprintf("Before call to IW_kNN_graph\n");
#endif

    //
    // Creating a kNN graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_res;
    if (version == 1) {
        IW_kNN_graph_res = IW_kNN_graph(RX, Rk);
        if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph function returned an invalid result (null pointer).");
    } else if (version == 2) {
        IW_kNN_graph_res = IW_kNN_graph_v2(RX, Rk);
        if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph_v2 function returned an invalid result (null pointer).");
    }

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("After call to IW_kNN_graph\n");
#endif

#if DEBUG__S_IW_kNN_graph
    Rprintf("\nIn S_IW_kNN_graph()\n");
    Rprintf("\nIW_kNN_graph:\n");
    for (int i = 0; i < IW_kNN_graph_res->size(); i++) {
        Rprintf("i: %d\n", i);
        for (auto neighbor_pair : (*IW_kNN_graph_res)[i])
            Rprintf("%d\t%d\n", neighbor_pair.first, neighbor_pair.second);
        Rprintf("\n");
    }
#endif

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    int IW_kNN_graph_n_vertices = IW_kNN_graph_res->size();
    for (int i = 0; i < IW_kNN_graph_n_vertices; i++) {
        for (auto nn_pair : IW_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_pair.first);
            (*intersection_size_vect)[i].push_back(nn_pair.second);
        }
    }

    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);

        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);

        for (auto dist : intersection_size_vect->at(i))
            *W++ = dist;

        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    //
    // Identifying connected components of the graph
    //
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("Before pruning\n");
#endif

    //
    // Pruning the graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph_res, prune_version);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_pair : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_pair.first);
        }
    }

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("After pruning\n");
#endif


#if DEBUG__S_IW_kNN_graph
    Rprintf("\nS_IW_kNN_graph(): pruned_graph_adj_vect:\n");
    for (int i = 0; i < pruned_graph_adj_vect.size(); i++) {
        Rprintf("i: %d\n", i);
        for (auto neighbor : pruned_graph_adj_vect[i])
            Rprintf("%d\t%d\n", neighbor);
        Rprintf("\n");
    }
#endif

    // Preparing pruned_I_kNN_graph adj list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("Before finding cycles\n");
#endif

    //
    // Finding the sizes of cycles within the pruned graph
    //
    #if 0
    std::unique_ptr<std::vector<int>> pruned_graph_cycle_lengths_ptr = cycle_sizes(pruned_graph_adj_vect);
    SEXP pruned_graph_cycle_lengths_Rvector = PROTECT(allocVector(INTSXP, pruned_graph_cycle_lengths_ptr->size())); nprot++;
    std::copy(pruned_graph_cycle_lengths_ptr->begin(), pruned_graph_cycle_lengths_ptr->end(), INTEGER(pruned_graph_cycle_lengths_Rvector));
    #endif

#if DEBUG__S_IW_kNN_graph_L1
    Rprintf("After finding cycles\n");
#endif

    SEXP res = PROTECT(allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, conn_comps);
    SET_VECTOR_ELT(res, 3, pruned_adj_list);
    //SET_VECTOR_ELT(res, 4, pruned_graph_cycle_lengths_Rvector);

    SEXP names = PROTECT(allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("conn_comps"));
    SET_STRING_ELT(names, 3, mkChar("pruned_adj_list"));
    //SET_STRING_ELT(names, 4, mkChar("pruned_cycles_sizes"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}

// ------------------------------------------------------------------------------------------
//
// IW_rho_kNN_graph functions
//
// ------------------------------------------------------------------------------------------

/**
 * @brief Constructs a rho-corrected intersection-weighted k-Nearest Neighbors (kNN) graph.
 *
 * This function builds a rho-corrected intersection-weighted kNN graph from a given data matrix.
 * For each point, it finds its k nearest neighbors with a correction factor rho, then constructs
 * edges between points based on the number of common nearest neighbors.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @param Rrho A numeric vector (SEXP) representing the correction factors for each point.
 * @return A unique pointer to a 2D vector of pairs representing the kNN graph. Each pair contains:
 *         - The index of the nearest neighbor.
 *         - The number of common elements in the kNN sets.
 *
 * Each point in the graph is connected to its nearest neighbors, with weights determined by the
 * intersection of their kNN sets. The distance between points is adjusted by a correction factor
 * rho, which normalizes the distances.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 *
 * @note The kNN search is performed with 2*k neighbors initially to allow for rho correction.
 * This function efficiently handles the case where the number of points is less than 2*k.
 */
std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = coerceVector(Rk, INTSXP)); nprot++;
    PROTECT(Rrho = coerceVector(Rrho, REALSXP)); nprot++;

    double *X = REAL(RX);
    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];
    int ncX = dimX[1];
    int k = INTEGER(Rk)[0];
    double *rho = REAL(Rrho);

    // Convert RX matrix to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nrX, ncX);
    for (int i = 0; i < nrX; i++) {
        for (int j = 0; j < ncX; j++) {
            data_pts[i][j] = X[i + nrX * j]; // Column-major to row-major
        }
    }

    // Finding kNN's for all points of X with K = 2 * k
    SEXP Rtwo_k = PROTECT(allocVector(INTSXP, 1)); nprot++;
    int two_k = 2 * k;
    if (two_k >= nrX) {
        two_k = nrX - 1;
    }
    INTEGER(Rtwo_k)[0] = two_k;

    // Finding kNN for K = 2*k
    SEXP knn_res = PROTECT(S_kNN(RX, Rtwo_k)); nprot++;
    int *two_k_indices = INTEGER(VECTOR_ELT(knn_res, 0));
    double *two_k_distances = REAL(VECTOR_ELT(knn_res, 1));

    // Indices of kNN's computed with d1 distance corrections
    std::vector<int> k_indices(k * nrX);
    std::vector<std::pair<int,double>> nn_idist(two_k); // ANN includes the query point as the first of the kNN's, we are sorting kNN's excluding the query point
    int nn_j;

    for (int i = 0; i < nrX; i++) {
        nn_idist.clear();  // Reset for the next point
        for (int j = 1; j < two_k; j++) {
            nn_j = two_k_indices[i + nrX * j];
            nn_idist.emplace_back(nn_j, two_k_distances[i + nrX * j] / rho[nn_j]);
        }
        // Sorting nn_idist in the ascending order with respect to the second component of each pair
        std::sort(nn_idist.begin(), nn_idist.end(),
                  [&](std::pair<int,double> a, std::pair<int,double> b) {
            return a.second < b.second;
        });

        // Populating k_indices with d1-corrected kNN's indices
        k_indices[i] = i;
        for (int j = 1; j < k; j++)
            k_indices[i + nrX * j] = nn_idist[j - 1].first;
    }

    auto res = std::make_unique<std::vector<std::vector<std::pair<int, int>>>>(nrX); // Initialize res with nrX empty vectors of int-int pairs

    // Perform k-NN search for each point
    int nrX_minus_one = nrX - 1;
    int common_count = 0;

    std::unordered_set<int> set_i;
    std::unordered_set<int> set_j;

    for (int pt_i = 0; pt_i < nrX_minus_one; pt_i++) {
        set_i.clear();
        // Populate set_i for pt_i
        for (int j = 0; j < k; j++)
            set_i.insert(k_indices[pt_i + nrX * j]);

        for (int pt_j = pt_i + 1; pt_j < nrX; pt_j++) {
            set_j.clear();
            for (int j = 0; j < k; j++)
                set_j.insert(k_indices[pt_j + nrX * j]);

            common_count = 0;
            // If set_i is smaller, iterate over it for efficiency
            if (set_i.size() < set_j.size()) {
                for (const int& elem : set_i) {
                    if (set_j.count(elem) > 0) {
                        common_count++;
                    }
                }
            } else {
                for (const int& elem : set_j) {
                    if (set_i.count(elem) > 0) {
                        common_count++;
                    }
                }
            }

            if (common_count > 0) {
                (*res)[pt_i].emplace_back(pt_j, common_count);
                (*res)[pt_j].emplace_back(pt_i, common_count);
            }
        }
    }

    // Clean up
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    UNPROTECT(nprot);

    return res;
}

/**
 * @brief Constructs a rho-corrected intersection-weighted k-Nearest Neighbors (kNN) graph and processes it.
 *
 * This function builds a rho-corrected intersection-weighted kNN graph from a given data matrix and processes
 * the graph to produce several outputs: adjacency list, intersection sizes, connected components, pruned
 * adjacency list, and cycle sizes in the pruned graph.
 *
 * @param RX A numeric matrix (SEXP) representing the data points (rows are points, columns are features).
 *           It is coerced to a numeric matrix if not already in the correct format.
 * @param Rk An integer (SEXP) representing the number of nearest neighbors to consider.
 * @param Rrho A numeric vector (SEXP) representing the correction factors for each point.
 * @return A list (SEXP) with the following components:
 *         - `adj_list`: A list of integer vectors, each representing the adjacency list for a point.
 *         - `isize_list`: A list of integer vectors, each representing the intersection size with neighbors.
 *         - `conn_comps`: An integer vector representing the connected components of the graph.
 *         - `pruned_adj_list`: A list of integer vectors, each representing the adjacency list of the pruned graph.
 *         - `pruned_cycles_sizes`: An integer vector representing the sizes of cycles within the pruned graph.
 *
 * The function performs the following steps:
 * 1. Constructs a rho-corrected intersection-weighted kNN graph using `IW_rho_kNN_graph`.
 * 2. Extracts the adjacency lists and intersection sizes from the graph.
 * 3. Identifies connected components of the graph using a union-find algorithm.
 * 4. Prunes long edges from the graph using `prune_long_edges`.
 * 5. Finds the sizes of cycles within the pruned graph using `cycle_sizes`.
 *
 * @note The function uses the ANN library for kNN computation and R's SEXP structures.
 * Ensure the ANN library is properly linked and included in your project.
 *
 * @note The kNN search is performed with 2*k neighbors initially to allow for rho correction.
 * This function efficiently handles the case where the number of points is less than 2*k.
 *
 * @throws Rf_error if RX cannot be coerced to a numeric matrix or if `IW_rho_kNN_graph` returns a null pointer.
 */
SEXP S_IW_rho_kNN_graph(SEXP RX, SEXP Rk, SEXP Rrho) {

    int nprot = 0;
    PROTECT(RX = coerceVector(RX, REALSXP)); nprot++;
    if (TYPEOF(RX) != REALSXP) {
        Rf_error("X could not be coerced to a numeric matrix. X has to be a numeric matrix (cannot be a data frame).");
    }

    int *dimX = INTEGER(getAttrib(RX, R_DimSymbol));
    int nrX = dimX[0];

    //
    // Creating a kNN graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> IW_kNN_graph_res = IW_rho_kNN_graph(RX, Rk, Rrho);
    if (!IW_kNN_graph_res) Rf_error("IW_kNN_graph function returned an invalid result (null pointer).");

    auto adj_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    auto intersection_size_vect = std::make_unique<std::vector<std::vector<int>>>(nrX);
    int IW_kNN_graph_n_vertices = IW_kNN_graph_res->size();
    for (int i = 0; i < IW_kNN_graph_n_vertices; i++) {
        for (auto nn_pair : IW_kNN_graph_res->at(i)) {
            (*adj_vect)[i].push_back(nn_pair.first);
            (*intersection_size_vect)[i].push_back(nn_pair.second);
        }
    }

    SEXP adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, adj_vect->at(i).size())); nprot++;
        int* A = INTEGER(RA);

        for (auto neighbor : adj_vect->at(i))
            *A++ = neighbor + 1;

        SET_VECTOR_ELT(adj_list, i, RA);
    }

    SEXP intersection_size_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RW = PROTECT(allocVector(INTSXP, intersection_size_vect->at(i).size())); nprot++;
        int* W = INTEGER(RW);

        for (auto dist : intersection_size_vect->at(i))
            *W++ = dist;

        SET_VECTOR_ELT(intersection_size_list, i, RW);
    }

    //
    // Identifying connected components of the graph
    //
    std::unique_ptr<std::vector<int>> conn_comps_ptr = union_find(adj_vect);
    SEXP conn_comps = PROTECT(allocVector(INTSXP, conn_comps_ptr->size())); nprot++;
    std::copy(conn_comps_ptr->begin(), conn_comps_ptr->end(), INTEGER(conn_comps));

    //
    // Pruning the graph
    //
    std::unique_ptr<std::vector<std::vector<std::pair<int, int>>>> pruned_graph = prune_long_edges(IW_kNN_graph_res, 1);
    auto pruned_graph_adj_vect = std::vector<std::vector<int>>(nrX);
    int pruned_graph_n_vertices = pruned_graph->size();
    for (int i = 0; i < pruned_graph_n_vertices; i++) {
        for (auto nn_pair : pruned_graph->at(i)) {
            pruned_graph_adj_vect[i].push_back(nn_pair.first);
        }
    }

    // Preparing pruned_I_kNN_graph adj list
    SEXP pruned_adj_list = PROTECT(allocVector(VECSXP, nrX)); nprot++;
    for (int i = 0; i < nrX; i++) {
        SEXP RA = PROTECT(allocVector(INTSXP, pruned_graph_adj_vect.at(i).size())); nprot++;
        int* A = INTEGER(RA);
        for (auto neighbor : pruned_graph_adj_vect.at(i))
            *A++ = neighbor + 1;
        SET_VECTOR_ELT(pruned_adj_list, i, RA);
    }

    //
    // Finding the sizes of cycles within the pruned graph
    //
    std::unique_ptr<std::vector<int>> pruned_graph_cycle_lengths_ptr = cycle_sizes(pruned_graph_adj_vect);
    SEXP pruned_graph_cycle_lengths_Rvector = PROTECT(allocVector(INTSXP, pruned_graph_cycle_lengths_ptr->size())); nprot++;
    std::copy(pruned_graph_cycle_lengths_ptr->begin(), pruned_graph_cycle_lengths_ptr->end(), INTEGER(pruned_graph_cycle_lengths_Rvector));


    SEXP res = PROTECT(allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(res, 0, adj_list);
    SET_VECTOR_ELT(res, 1, intersection_size_list);
    SET_VECTOR_ELT(res, 2, conn_comps);
    SET_VECTOR_ELT(res, 3, pruned_adj_list);
    SET_VECTOR_ELT(res, 4, pruned_graph_cycle_lengths_Rvector);

    SEXP names = PROTECT(allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(names, 0, mkChar("adj_list"));
    SET_STRING_ELT(names, 1, mkChar("isize_list"));
    SET_STRING_ELT(names, 2, mkChar("conn_comps"));
    SET_STRING_ELT(names, 3, mkChar("pruned_adj_list"));
    SET_STRING_ELT(names, 4, mkChar("pruned_cycles_sizes"));

    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return res;
}
