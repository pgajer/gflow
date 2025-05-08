#include <R.h>
#include <Rinternals.h>
// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <unordered_set>
#include <limits>
#include <ANN/ANN.h>

#include "cpp_utils.hpp"

extern "C" {
    SEXP S_mstree(SEXP X);
}

struct edge_t {
    int from;
    int to;
    double dist;
    edge_t(int f, int t, double d) : from(f), to(t), dist(d) {}
};

struct compare_edge_t {
    bool operator()(const edge_t& e1, const edge_t& e2) {
        return e1.dist > e2.dist;
    }
};

/**
 * @brief Prints the contents of a priority queue.
 *
 * This function creates a copy of the input priority queue and prints its elements
 * in the order they would be popped (highest priority first). The original queue
 * remains unmodified.
 *
 * @tparam T The type of elements in the priority queue.
 * @tparam Container The container type used by the priority queue (default: std::vector<T>).
 * @tparam Compare The comparison function type used by the priority queue.
 *
 * @param pq The priority queue to print. Passed by const reference to prevent modification.
 * @param name An optional name to be printed before the queue contents (default: "").
 * @param n The maximum number of elements to print. If 0 or greater than queue size,
 *          prints all elements (default: 0).
 *
 * @note This function creates a copy of the entire queue, which may be memory-intensive
 *       for very large queues.
 *
 * @example
 *     std::priority_queue<int> my_pq({3, 1, 4, 1, 5, 9});
 *     print_pq(my_pq, "My Queue", 3);
 *     // Output: My Queue: 9 5 4
 */

// Add this operator overload
std::ostream& operator<<(std::ostream& os, const edge_t& edge) {
    os << "(" << edge.from << "->" << edge.to << ": " << edge.dist << ")";
    return os;
}

void print_edge(const edge_t& edge) {
    std::cout << "(" << edge.from << "->" << edge.to << ": " << edge.dist << ")";
}

template <typename T, typename Container, typename Compare>
void print_pq(const std::priority_queue<T, Container, Compare>& pq, const std::string& name = "", size_t n = 0) {
    std::priority_queue<T, Container, Compare> pq_copy = pq;
    if (n == 0) {
        n = pq.size();
    }
    if (!name.empty()) {
        std::cout << name << ": ";
    }
    size_t count = 0;
    while (!pq_copy.empty() && count < n) {
        print_edge(pq_copy.top());  // Use the helper function here
        std::cout << " ";
        pq_copy.pop();
        ++count;
    }
    std::cout << std::endl;
}


/**
 * @brief Computes the Minimum Spanning Tree (MST) using Prim's algorithm with ANN library for efficient nearest neighbor searches.
 *
 * This function implements Prim's algorithm to find the Minimum Spanning Tree of a given set of points in a multi-dimensional space.
 * It uses the ANN (Approximate Nearest Neighbor) library to perform efficient nearest neighbor searches, which is crucial for
 * the performance of the algorithm, especially in high-dimensional spaces.
 *
 * @details The algorithm works as follows:
 * 1. Convert the input data to ANN's point array format and build a kd-tree for efficient searching.
 * 2. Initialize the MST with the first point (index 0).
 * 3. Repeatedly find the shortest edge connecting a point in the MST to a point not in the MST:
 *    a. For each point in the MST, find its nearest neighbor that is not in the MST.
 *    b. Among these edges, choose the shortest one and add it to the MST.
 * 4. Repeat step 3 until all points are in the MST.
 *
 * The implementation uses a priority queue to efficiently select the next edge to add to the MST.
 * It also uses an unordered set to keep track of the vertices currently in the MST, allowing for
 * efficient checking of whether a point is in the MST or not.
 *
 * @param X A vector of doubles representing the flat data matrix. Each nc_X consecutive elements represent one point.
 * @param nr_X The number of points (rows) in the data matrix.
 * @param nc_X The number of dimensions (columns) for each point in the data matrix.
 *
 * @return A vector of edge_t structures representing the edges in the Minimum Spanning Tree.
 *         Each edge_t contains the indices of the two points it connects and the distance between them.
 *
 * @note The function uses Euclidean distance as the distance metric. If a different metric is needed,
 *       modifications to the distance calculation would be required.
 *
 * @warning This function assumes that the input data is valid and contains at least one point.
 *          It does not perform extensive error checking on the input parameters.
 *
 * @see edge_t, compare_edge_t
 *
 * Time Complexity: O(n^2 log n), where n is the number of points.
 * Space Complexity: O(n), where n is the number of points.
 *
 * @todo Consider implementing a more efficient version that doesn't recompute all nearest neighbors
 *       after each point is added to the MST.
 */
std::vector<edge_t> data_mstree(const std::vector<double>& X, int nr_X, int nc_X) {

    // Convert flat X to ANNpointArray
    ANNpointArray data_pts = annAllocPts(nr_X, nc_X);
    for (int i = 0; i < nr_X; i++) {
        for (int j = 0; j < nc_X; j++) {
            data_pts[i][j] = X[i +  nr_X * j];
        }
    }

    // Build kd-tree
    ANNkd_tree* kd_tree = new ANNkd_tree(data_pts, nr_X, nc_X);

    // Initialize MST algorithm
    std::vector<bool> inMST(nr_X, false);
    std::unordered_set<int> tree_vertices;
    std::vector<edge_t> result;
    std::priority_queue<edge_t, std::vector<edge_t>, compare_edge_t> pq;

    // Variables for nearest neighbor search
    ANNidxArray nn_idx = new ANNidx[nr_X];  // Allocate for maximum possible k
    ANNdistArray nn_dist = new ANNdist[nr_X];

    // Start from the first point
    inMST[0] = true;
    tree_vertices.insert(0);

    // Function to add edges from points in MST to points not in MST
    auto add_edges = [&]() {
        for (const auto &v : tree_vertices) {
            kd_tree->annkSearch(data_pts[v], nr_X, nn_idx, nn_dist, 0);
            for (int i = 1; i < nr_X; i++) {  // Start from 1 to skip the query point itself
                if (!inMST[nn_idx[i]]) {
                    pq.push(edge_t(v, nn_idx[i], std::sqrt(nn_dist[i])));
                    break;  // We only need the closest non-MST point
                }
            }
        }
    };

    // Add edges from the first point
    add_edges();

    while (!pq.empty() && result.size() < nr_X - 1) {
        edge_t e = pq.top();
        pq.pop();

        if (inMST[e.to]) continue;

        // Add this edge to MST
        result.push_back(e);
        inMST[e.to] = true;
        tree_vertices.insert(e.to);

        // Add edges from the newly updated MST
        add_edges();
    }

    // Clean up
    delete[] nn_idx;
    delete[] nn_dist;
    delete kd_tree;
    annDeallocPts(data_pts);
    annClose(); // Close ANN

    return result;
}

/**
 * @brief R interface for computing Minimum Spanning Tree using Prim's algorithm.
 *
 * This function serves as an interface between R and the C++ implementation of Prim's algorithm
 * for computing the Minimum Spanning Tree (MST). It takes an R matrix as input, computes the MST,
 * and returns the result as an R matrix.
 *
 * @param X An R matrix (SEXP) where each column represents a point in the dataset.
 *          The matrix should be numeric (double precision).
 *
 * @return An R matrix (SEXP) with three columns:
 *         1. start: The index of the starting vertex of each edge (1-based for R compatibility)
 *         2. end: The index of the ending vertex of each edge (1-based for R compatibility)
 *         3. length: The length (distance) of each edge
 *         The number of rows in the output matrix is equal to the number of edges in the MST,
 *         which is always one less than the number of input points.
 *
 * @note This function uses the data_mstree C++ function to compute the MST.
 *       It handles the conversion between R and C++ data structures.
 *
 * @warning The input matrix X must be a numeric matrix. The function will raise an R error
 *          if this condition is not met.
 *
 * @see data_mstree
 *
 * @example
 * In R:
 * ```R
 * # Assuming the shared object is loaded
 * X <- matrix(runif(200), ncol=2)  # 100 random 2D points
 * result <- .Call("S_mstree", X)
 * ```
 *
 * @note The resulting edge indices are 1-based to conform with R's indexing convention.
 *
 * @note This function uses R's C API and should be compiled into a shared object
 *       before being called from R.
 *
 * @todo Consider adding more robust error checking and handling for production use.
 *
 * @todo If performance is critical, consider using Rcpp for easier and potentially
 *       more efficient R and C++ integration.
 */
SEXP S_mstree(SEXP X) {
    // Check if X is a numeric matrix
    if (!isMatrix(X) || !isReal(X)) {
        error("Input must be a numeric matrix");
    }

    // Get dimensions of X
    SEXP dim = getAttrib(X, R_DimSymbol);
    int nr_X = INTEGER(dim)[0];
    int nc_X = INTEGER(dim)[1];

    // Convert X to std::vector<double>
    std::vector<double> x_vec(REAL(X), REAL(X) + nr_X * nc_X);

    // Call data_mstree
    std::vector<edge_t> mst = data_mstree(x_vec, nr_X, nc_X);

    // Create result matrix
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, mst.size(), 3));
    double *result_ptr = REAL(result);

    // Fill result matrix
    for (size_t i = 0; i < mst.size(); ++i) {
        result_ptr[i] = mst[i].from + 1; // R uses 1-based indexing
        result_ptr[i + mst.size()] = mst[i].to + 1;
        result_ptr[i + 2 * mst.size()] = mst[i].dist;
    }

    // Set column names
    SEXP colnames;
    PROTECT(colnames = allocVector(STRSXP, 3));
    SET_STRING_ELT(colnames, 0, mkChar("start"));
    SET_STRING_ELT(colnames, 1, mkChar("end"));
    SET_STRING_ELT(colnames, 2, mkChar("length"));
    setAttrib(result, R_DimNamesSymbol,
              PROTECT(list2(R_NilValue, colnames)));

    UNPROTECT(3);
    return result;
}

