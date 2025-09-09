// Initial version of diffusion smoother

#include <R.h>
#include <Rinternals.h>

// Undefine conflicting macros after including R headers
#undef length

#include <vector>
#include <queue>
#include <memory>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <stack>
#include <numeric>

#include "msr2.h"
#include "msr2_cpp_utils.h"
#include "msr2_Cpp_to_R_utils.h"


double calculate_smallest_difference(const std::vector<double>& y);
std::unique_ptr<std::pair<std::vector<int>, std::vector<int>>> find_local_extrema(const std::vector<std::vector<int>>& graph,
                                                                                  const std::vector<double>& y);
std::unique_ptr<std::vector<int>> loc_const_vertices(const std::vector<std::vector<int>>& graph,
                                                     const std::vector<double>& y,
                                                     double prec);
bool has_constant_vertex(const std::vector<std::vector<int>>& graph,
                         const std::vector<double>& y,
                         double prec);
bool has_constant_vertex(const std::vector<std::vector<int>>& graph,
                         const std::vector<double>& y,
                         double prec,
                         const std::vector<int>& indices_to_check);

extern "C" {

    SEXP S_heat_smoother(SEXP Rgraph,
                         SEXP Ry,
                         SEXP Rexclude,
                         SEXP Rstep_factor,
                         SEXP Rprec,
                         SEXP Rn_itrs,
                         SEXP Rstop_criterion);

    SEXP S_weighted_heat_smoother(SEXP Rgraph,
                                  SEXP Ry,
                                  SEXP Rweights,
                                  SEXP Rstep_factor,
                                  SEXP Rprec,
                                  SEXP Rn_itrs,
                                  SEXP Rstop_criterion);

    SEXP S_weighted_heat_smoother_trajectory(SEXP Rgraph,
                                             SEXP Ry,
                                             SEXP Rweights,
                                             SEXP Rstep_factor,
                                             SEXP Rprec,
                                             SEXP Rn_itrs,
                                             SEXP Rstop_criterion,
                                             SEXP Rmean_adjust);
}


// ------------------------------------------------------------------------------------------
//
// Heat equation smoothing
//
// ------------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------------
//
//  Heat Smoother
//
// ------------------------------------------------------------------------------------------

/**
 * Smooths a graph function using the heat diffusion equation.
 *
 * This function smooths the input function `y` defined on the vertices of a graph using the heat
 * diffusion equation. The smoothing process iteratively adjusts the function values at each vertex
 * based on the average value of its neighboring vertices. The process can be stopped based on two
 * criteria: either when the smoothed function is not locally constant at any vertex (determined by
 * the precision threshold `prec`), or when a specified number of iterations (`n_itrs`) is reached.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph.
 * @param exclude Indices of vertices to exclude from the smoothing process. The function values at these
 *                vertices will remain unchanged. This can be used to preserve local extrema of the function
 *                and maintain the range of the function during smoothing.
 * @param step_factor A positive number used to calculate the step size. The step size is the product of
 *                    `step_factor` and the smallest difference between the unique values of `y`.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference between
 *             the value at a vertex and any of its neighbors is greater than this threshold, the vertex is
 *             not considered locally constant. This parameter is used only when `stop_criterion` is set to 1.
 *             Default value is 1e-8.
 * @param n_itrs A positive integer indicating the maximum number of iterations for the smoothing process.
 *               This parameter is used only when `stop_criterion` is set to 2. Default value is 10.
 * @param stop_criterion An integer indicating the type of stopping criterion for the smoothing process.
 *                       If set to 1 (default), the process stops when the smoothed function is not locally
 *                       constant at any vertex, determined by the `prec` parameter. If set to 2, the process
 *                       stops after `n_itrs` iterations.
 * @return A unique pointer to the smoothed version of the input function.
 */
std::unique_ptr<std::vector<double>> heat_smoother(const std::vector<std::vector<int>>& graph,
                                                   const std::vector<double>& y,
                                                   const std::vector<int>& exclude,
                                                   double step_factor,
                                                   double prec = 1e-8,
                                                   int n_itrs = 10,
                                                   int stop_criterion = 1) {
    int n_vertices = y.size();

    if (n_vertices != (int)graph.size())
        error("The lengths of graph and y must be the same.");

    if (prec < 0)
        error("prec has to be a positive number.");

    if (n_itrs < 1)
        error("n_itrs has to be greater than or equal to 1.");

    if (stop_criterion < 0 || stop_criterion > 2)
        error("stop_criterion can only take values 1 or 2.");

    int n_vertices_minus_one = n_vertices - 1;
    std::vector<int> all_indices(n_vertices);
    std::iota(all_indices.begin(), all_indices.end(), 0);
    std::vector<int> indices_to_process;

    if (!exclude.empty()) {
        // Checking that all values of exclude (if non-empty) are within the range between 0 and (n_vertices - 1)
        for (const auto& val : exclude)
            if (val < 0 || val > n_vertices_minus_one) {
                std::cerr << "Found value " << val << std::endl;
                std::cerr << "exclude values need to be between 0 and " << n_vertices_minus_one << std::endl;
                error("Invalid exclude value.");
            }

        // Setting the set of indices to process to be the set difference of all_indices and exclude.
        std::set_difference(all_indices.begin(), all_indices.end(),
                            exclude.begin(), exclude.end(),
                            std::back_inserter(indices_to_process));
    } else {
        std::copy(all_indices.begin(), all_indices.end(), std::back_inserter(indices_to_process));
    }

    double step_size = calculate_smallest_difference(y);
    step_size *= step_factor;

    // computing the mean of y
    double Ey = 0;
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        Ey += y[vertex];
    Ey /= n_vertices;

    auto y_t = std::make_unique<std::vector<double>>(n_vertices); // y at time t
    std::copy(y.begin(), y.end(), y_t->begin()); // Initializing y_t with the starting values

    auto n_neighbors = std::vector<double>(n_vertices);
    for (const auto& vertex : indices_to_process)
        n_neighbors[vertex] = graph[vertex].size();

    double average_neighbor_value = 0;
    double Ey_t = 0;
    double Delta;
    if (stop_criterion == 1) {

        // Identifying constant vertices
        std::unique_ptr<std::vector<int>> loc_const_vertices_vect = loc_const_vertices(graph, y, prec);

        std::vector<int> sorted_vertices(n_vertices);

        if ( (*loc_const_vertices_vect).empty() ) {
            return y_t;
        } else {
            // Sorting indices so that constant vertices are at the beginning of the sorted vector
            auto const_vertices = *loc_const_vertices_vect;

            std::vector<int> all_indices(n_vertices);
            std::iota(all_indices.begin(), all_indices.end(), 0);
            std::vector<int> non_const_indices;

            std::set_difference(all_indices.begin(), all_indices.end(),
                                const_vertices.begin(), const_vertices.end(),
                                std::back_inserter(non_const_indices));

            // Concatenate the vectors, maintaining the order where constant vertices come first
            std::copy(const_vertices.begin(),  const_vertices.end(), sorted_vertices.begin());
            std::copy(non_const_indices.begin(), non_const_indices.end(), sorted_vertices.begin() + const_vertices.size());
        }

        while (has_constant_vertex(graph, *y_t, prec, sorted_vertices)) {
            for (const auto& vertex : indices_to_process) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += (*y_t)[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                (*y_t)[vertex] += step_size * (average_neighbor_value - (*y_t)[vertex]);
            }

            // Adjusting y_t so that its mean = mean(y)
            Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += (*y_t)[vertex];
            Ey_t /= n_vertices;

            Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                (*y_t)[vertex] += Delta;
        }
    } else if (stop_criterion == 2) {
        int itr = 0;
        while (itr < n_itrs) {
            for (const auto& vertex : indices_to_process) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += (*y_t)[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                (*y_t)[vertex] += step_size * (average_neighbor_value - (*y_t)[vertex]);
            }

            // Adjusting y_t so that its mean = mean(y)
            Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += (*y_t)[vertex];
            Ey_t /= n_vertices;

            Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                (*y_t)[vertex] += Delta;

            itr++;
        }
    }

    return y_t;
}

/**
 * An R interface to the heat_smoother function.
 *
 * This function takes an R list representing a graph adjacency list, a numeric vector
 * of function values at each vertex, a vector of indices to exclude from smoothing,
 * and various parameters for the smoothing process. It calls the C++ heat_smoother
 * function to smooth the function values using the heat diffusion equation and returns
 * the smoothed function values as an R numeric vector.
 *
 * @param Rgraph An R list representing the graph adjacency list. Each element of the list
 *               should be an integer vector specifying the neighboring vertex indices for
 *               a given vertex.
 * @param Ry A numeric vector representing the function values at each vertex of the graph.
 * @param Rexclude An integer vector specifying the indices of vertices to exclude from the
 *                 smoothing process. The function values at these vertices will remain unchanged.
 * @param Rstep_factor A numeric scalar specifying the step factor used to calculate the step size
 *                     in the smoothing process. The step size is the product of the step factor
 *                     and the smallest difference between the unique values of the function.
 * @param Rprec A numeric scalar specifying the precision threshold used to determine local
 *              constancy. If the absolute difference between the value at a vertex and any
 *              of its neighbors is greater than this threshold, the vertex is not considered
 *              locally constant.
 * @param Rn_itrs An integer scalar specifying the maximum number of iterations for the smoothing
 *                process. This parameter is used only when the stopping criterion is set to 2.
 * @param Rstop_criterion An integer scalar specifying the stopping criterion for the smoothing process.
 *                        A value of 1 indicates that the process should stop when there are no locally
 *                        constant vertices remaining. A value of 2 indicates that the process should stop
 *                        after the specified number of iterations (Rn_itrs).
 *
 * @return An R numeric vector containing the smoothed function values at each vertex of the graph.
 */
SEXP S_heat_smoother(SEXP Rgraph,
                     SEXP Ry,
                     SEXP Rexclude,
                     SEXP Rstep_factor,
                     SEXP Rprec,
                     SEXP Rn_itrs,
                     SEXP Rstop_criterion) {

    std::vector<std::vector<int>> graph_adj_list = std::move(*Rlist_to_Cpp(Rgraph));

    int n_vertices = LENGTH(Ry);

    int nprot = 0;
    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    double *y = REAL(Ry);

    PROTECT(Rexclude = coerceVector(Rexclude, INTSXP)); nprot++;
    int *exclude = INTEGER(Rexclude);

    PROTECT(Rstep_factor = coerceVector(Rstep_factor, REALSXP)); nprot++;
    double step_factor = REAL(Rstep_factor)[0];

    PROTECT(Rprec = coerceVector(Rprec, REALSXP)); nprot++;
    double prec = REAL(Rprec)[0];

    PROTECT(Rn_itrs = coerceVector(Rn_itrs, INTSXP)); nprot++;
    int n_itrs = INTEGER(Rn_itrs)[0];

    PROTECT(Rstop_criterion = coerceVector(Rstop_criterion, INTSXP)); nprot++;
    int stop_criterion = INTEGER(Rstop_criterion)[0];

    std::unique_ptr<std::vector<double>> y_smoothed = heat_smoother(graph_adj_list,
                                                                    std::vector<double>(y, y + n_vertices),
                                                                    std::vector<int>(exclude, exclude + LENGTH(Rexclude)),
                                                                    step_factor,
                                                                    prec,
                                                                    n_itrs,
                                                                    stop_criterion);

    SEXP result = PROTECT(allocVector(REALSXP, n_vertices)); nprot++;
    std::copy(y_smoothed->begin(), y_smoothed->end(), REAL(result));

    UNPROTECT(nprot);

    return result;
}

// ------------------------------------------------------------------------------------------
//
// Weighted Heat Smoother
//
// ------------------------------------------------------------------------------------------

/**
 * Smooths a graph function using the weighted heat diffusion equation.
 *
 * This function smooths the input function `y` defined on the vertices of a graph using the weighted
 * heat diffusion equation. The smoothing process iteratively adjusts the function values at each vertex
 * based on the weighted average of its neighboring vertices. The process can be stopped based on two
 * criteria: either when the smoothed function is not locally constant at any vertex (determined by
 * the precision threshold `prec`), or when a specified number of iterations (`n_itrs`) is reached.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the function values at each vertex of the graph.
 * @param weights A vector of vertex weights used in the weighted smoothing process. The length of `weights`
 *                must be equal to the number of vertices in the graph.
 * @param step_factor A positive number used to calculate the step size. The step size is the product of
 *                    `step_factor` and the smallest difference between the unique values of `y`.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference between
 *             the value at a vertex and any of its neighbors is greater than this threshold, the vertex is
 *             not considered locally constant. This parameter is used only when `stop_criterion` is set to 1.
 *             Default value is 1e-8.
 * @param n_itrs A positive integer indicating the maximum number of iterations for the smoothing process.
 *               This parameter is used only when `stop_criterion` is set to 2. Default value is 10.
 * @param stop_criterion An integer indicating the type of stopping criterion for the smoothing process.
 *                       If set to 1 (default), the process stops when the smoothed function is not locally
 *                       constant at any vertex, determined by the `prec` parameter. If set to 2, the process
 *                       stops after `n_itrs` iterations.
 * @return A unique pointer to the smoothed version of the input function.
 */
std::unique_ptr<std::vector<double>> weighted_heat_smoother(const std::vector<std::vector<int>>& graph,
                                                            const std::vector<double>& y,
                                                            const std::vector<double>& weights,
                                                            double step_factor,
                                                            double prec = 1e-8,
                                                            int n_itrs = 10,
                                                            int stop_criterion = 1) {
    int n_vertices = y.size();

    if (n_vertices != (int)graph.size())
        error("The lengths of graph and y must be the same.");

    if (n_vertices != (int)weights.size())
        error("The lengths of weights and y must be the same.");

    if (prec < 0)
        error("prec has to be a positive number.");

    if (n_itrs < 1)
        error("n_itrs has to be greater than or equal to 1.");

    if (stop_criterion < 0 || stop_criterion > 2)
        error("stop_criterion can only take values 1 or 2.");

    double step_size = calculate_smallest_difference(y);
    step_size *= step_factor;

    // computing the mean of y
    double Ey = 0;
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        Ey += y[vertex];
    Ey /= n_vertices;

    auto y_t = std::make_unique<std::vector<double>>(n_vertices); // y at time t
    std::copy(y.begin(), y.end(), y_t->begin()); // Initializing y_t to y

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    double average_neighbor_value = 0;
    double Ey_t = 0;
    double Delta;
    if (stop_criterion == 1) {

        // Identifying constant vertices
        std::unique_ptr<std::vector<int>> loc_const_vertices_vect = loc_const_vertices(graph, y, prec);

        std::vector<int> sorted_vertices(n_vertices);

        if ( (*loc_const_vertices_vect).empty() ) {
            return y_t;
        } else {
            // Sorting indices so that constant vertices are at the beginning of the sorted vector
            auto const_vertices = *loc_const_vertices_vect;

            std::vector<int> all_indices(n_vertices);
            std::iota(all_indices.begin(), all_indices.end(), 0);
            std::vector<int> non_const_indices;

            std::set_difference(all_indices.begin(), all_indices.end(),
                                const_vertices.begin(), const_vertices.end(),
                                std::back_inserter(non_const_indices));

            // Concatenate the vectors, maintaining the order where constant vertices come first
            std::copy(const_vertices.begin(),  const_vertices.end(), sorted_vertices.begin());
            std::copy(non_const_indices.begin(), non_const_indices.end(), sorted_vertices.begin() + const_vertices.size());
        }

        while (has_constant_vertex(graph, *y_t, prec, sorted_vertices)) {
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += (*y_t)[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                (*y_t)[vertex] += weights[vertex] * step_size * (average_neighbor_value - (*y_t)[vertex]);
            }

            // Adjusting y_t so that its mean = mean(y)
            Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += (*y_t)[vertex];
            Ey_t /= n_vertices;

            Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                (*y_t)[vertex] += Delta;
        }
    } else if (stop_criterion == 2) {
        int itr = 0;
        while (itr < n_itrs) {
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += (*y_t)[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                (*y_t)[vertex] += weights[vertex] * step_size * (average_neighbor_value - (*y_t)[vertex]);
            }

            // Adjusting y_t so that its mean = mean(y)
            Ey_t = 0;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                Ey_t += (*y_t)[vertex];
            Ey_t /= n_vertices;

            Delta = Ey - Ey_t;
            for (int vertex = 0; vertex < n_vertices; ++vertex)
                (*y_t)[vertex] += Delta;

            itr++;
        }
    }

    return y_t;
}

/**
 * An R interface to the weighted_heat_smoother function.
 *
 * This function takes an R list representing a graph adjacency list, a numeric vector
 * of function values at each vertex, a numeric vector of vertex weights, and various
 * parameters for the smoothing process. It calls the C++ weighted_heat_smoother function
 * to smooth the function values using the weighted heat diffusion equation and returns
 * the smoothed function values as an R numeric vector.
 *
 * @param Rgraph An R list representing the graph adjacency list. Each element of the list
 *               should be an integer vector specifying the neighboring vertex indices for
 *               a given vertex.
 * @param Ry A numeric vector representing the function values at each vertex of the graph.
 * @param Rweights A numeric vector representing the weights of each vertex in the graph.
 *                 The length of Rweights must be equal to the number of vertices.
 * @param Rstep_factor A numeric scalar specifying the step factor used to calculate the step size
 *                     in the smoothing process. The step size is the product of the step factor
 *                     and the smallest difference between the unique values of the function.
 * @param Rprec A numeric scalar specifying the precision threshold used to determine local
 *              constancy. If the absolute difference between the value at a vertex and any
 *              of its neighbors is greater than this threshold, the vertex is not considered
 *              locally constant.
 * @param Rn_itrs An integer scalar specifying the maximum number of iterations for the smoothing
 *                process. This parameter is used only when the stopping criterion is set to 2.
 * @param Rstop_criterion An integer scalar specifying the stopping criterion for the smoothing process.
 *                        A value of 1 indicates that the process should stop when there are no locally
 *                        constant vertices remaining. A value of 2 indicates that the process should stop
 *                        after the specified number of iterations (Rn_itrs).
 *
 * @return An R numeric vector containing the smoothed function values at each vertex of the graph.
 */
SEXP S_weighted_heat_smoother(SEXP Rgraph,
                              SEXP Ry,
                              SEXP Rweights,
                              SEXP Rstep_factor,
                              SEXP Rprec,
                              SEXP Rn_itrs,
                              SEXP Rstop_criterion) {
    std::vector<std::vector<int>> graph_adj_list = std::move(*Rlist_to_Cpp(Rgraph));
    int n_vertices = LENGTH(Ry);

    int nprot = 0;

    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    double *y = REAL(Ry);

    PROTECT(Rweights = coerceVector(Rweights, REALSXP)); nprot++;
    double *weights = REAL(Rweights);

    PROTECT(Rstep_factor = coerceVector(Rstep_factor, REALSXP)); nprot++;
    double step_factor = REAL(Rstep_factor)[0];

    PROTECT(Rprec = coerceVector(Rprec, REALSXP)); nprot++;
    double prec = REAL(Rprec)[0];

    PROTECT(Rn_itrs = coerceVector(Rn_itrs, INTSXP)); nprot++;
    int n_itrs = INTEGER(Rn_itrs)[0];

    PROTECT(Rstop_criterion = coerceVector(Rstop_criterion, INTSXP)); nprot++;
    int stop_criterion = INTEGER(Rstop_criterion)[0];

    std::unique_ptr<std::vector<double>> y_smoothed = weighted_heat_smoother(graph_adj_list,
                                                                             std::vector<double>(y, y + n_vertices),
                                                                             std::vector<double>(weights, weights + n_vertices),
                                                                             step_factor,
                                                                             prec,
                                                                             n_itrs,
                                                                             stop_criterion);

    SEXP result = PROTECT(allocVector(REALSXP, n_vertices)); nprot++;
    std::copy(y_smoothed->begin(), y_smoothed->end(), REAL(result));

    UNPROTECT(nprot);
    return result;
}



// ------------------------------------------------------------------------------------------
//
// Weighted Heat Smoother Trajectory
//
// ------------------------------------------------------------------------------------------

/**
 * Smooths a graph function using the weighted heat diffusion equation and returns the trajectory of the function over time.
 *
 * This function smooths the input function `y` defined on the vertices of a graph using the weighted
 * heat diffusion equation. The smoothing process iteratively adjusts the function values at each vertex
 * based on the weighted average of its neighboring vertices. The process can be stopped based on two
 * criteria: either when the smoothed function is not locally constant at any vertex (determined by
 * the precision threshold `prec`), or when a specified number of iterations (`n_itrs`) is reached.
 *
 * The function returns the entire trajectory of the smoothed function over time, capturing the intermediate
 * states of the function at each iteration of the smoothing process. The trajectory is represented as a
 * vector of vectors, where each inner vector corresponds to the function values at a specific iteration.
 *
 * @param graph The graph represented as an adjacency list, where each element of the outer vector
 *              represents a vertex, and the inner vector contains the indices of its neighboring vertices.
 * @param y A vector representing the initial function values at each vertex of the graph.
 * @param weights A vector of vertex weights used in the weighted smoothing process. The length of `weights`
 *                must be equal to the number of vertices in the graph.
 * @param step_factor A positive number used to calculate the step size. The step size is the product of
 *                    `step_factor` and the smallest difference between the unique values of `y`.
 * @param prec The precision threshold used to determine local constancy. If the absolute difference between
 *             the value at a vertex and any of its neighbors is greater than this threshold, the vertex is
 *             not considered locally constant. This parameter is used only when `stop_criterion` is set to 1.
 *             Default value is 1e-8.
 * @param n_itrs A positive integer indicating the maximum number of iterations for the smoothing process.
 *               This parameter is used in both `stop_criterion`. Default value is 1000.
 * @param stop_criterion An integer indicating the type of stopping criterion for the smoothing process.
 *                       If set to 1 (default), the process stops when the smoothed function is not locally
 *                       constant at any vertex, determined by the `prec` parameter. If set to 2, the process
 *                       stops after `n_itrs` iterations.
 * @param mean_adjust A boolean scalar specifying whether to adjust the mean of the smoothed function
 *                     values at each iteration to match the mean of the initial function values.
 *                     If Rmean_adjust is true, the mean adjustment is performed; otherwise, it is skipped.
 *
 * @return A unique pointer to a vector of vectors representing the trajectory of the smoothed function over time.
 *         Each inner vector corresponds to the function values at a specific iteration of the smoothing process.
 */
std::unique_ptr<std::vector<std::vector<double>>> weighted_heat_smoother_trajectory(const std::vector<std::vector<int>>& graph,
                                                                                    const std::vector<double>& y,
                                                                                    const std::vector<double>& weights,
                                                                                    double step_factor,
                                                                                    double prec = 1e-8,
                                                                                    int n_itrs = 10,
                                                                                    int stop_criterion = 1,
                                                                                    bool mean_adjust = true) {
    int n_vertices = y.size();

    if (n_vertices != (int)graph.size())
        error("The lengths of graph and y must be the same.");

    if (n_vertices != (int)weights.size())
        error("The lengths of weights and y must be the same.");

    if (prec < 0)
        error("prec has to be a positive number.");

    if (n_itrs < 1)
        error("n_itrs has to be greater than or equal to 1.");

    if (stop_criterion < 0 || stop_criterion > 2)
        error("stop_criterion can only take values 1 or 2.");

    double step_size = calculate_smallest_difference(y);
    step_size *= step_factor;

    // computing the mean of y
    double Ey = 0;
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        Ey += y[vertex];
    Ey /= n_vertices;

    auto y_traj = std::make_unique<std::vector<std::vector<double>>>();
    auto y_t = y; // Initialize y_t to y

    auto n_neighbors = std::vector<double>(n_vertices);
    for (int vertex = 0; vertex < n_vertices; ++vertex)
        n_neighbors[vertex] = graph[vertex].size();

    double average_neighbor_value = 0;
    double Ey_t = 0;
    double Delta;
    if (stop_criterion == 1) {
        // Identifying constant vertices
        std::unique_ptr<std::vector<int>> loc_const_vertices_vect = loc_const_vertices(graph, y, prec);

        std::vector<int> sorted_vertices(n_vertices);

        if ((*loc_const_vertices_vect).empty()) {
            y_traj->push_back(y_t);
            return y_traj;
        } else {
            // Sorting indices so that constant vertices are at the beginning of the sorted vector
            auto const_vertices = *loc_const_vertices_vect;

            std::vector<int> all_indices(n_vertices);
            std::iota(all_indices.begin(), all_indices.end(), 0);
            std::vector<int> non_const_indices;

            std::set_difference(all_indices.begin(), all_indices.end(),
                                const_vertices.begin(), const_vertices.end(),
                                std::back_inserter(non_const_indices));

            // Concatenate the vectors, maintaining the order where constant vertices come first
            std::copy(const_vertices.begin(), const_vertices.end(), sorted_vertices.begin());
            std::copy(non_const_indices.begin(), non_const_indices.end(), sorted_vertices.begin() + const_vertices.size());
        }

        while (has_constant_vertex(graph, y_t, prec, sorted_vertices)) {
            y_traj->push_back(y_t);
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += y_t[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                y_t[vertex] += weights[vertex] * step_size * (average_neighbor_value - y_t[vertex]);
            }

            if (mean_adjust) {
                // Adjusting y_t so that its mean = mean(y)
                Ey_t = 0;
                for (int vertex = 0; vertex < n_vertices; ++vertex)
                    Ey_t += y_t[vertex];
                Ey_t /= n_vertices;

                Delta = Ey - Ey_t;
                for (int vertex = 0; vertex < n_vertices; ++vertex)
                    y_t[vertex] += Delta;
            }
        }
    } else if (stop_criterion == 2) {
        int itr = 0;
        while (itr < n_itrs) {
            y_traj->push_back(y_t);
            for (int vertex = 0; vertex < n_vertices; ++vertex) {
                average_neighbor_value = 0;
                for (int neighbor : graph[vertex])
                    average_neighbor_value += y_t[neighbor];
                average_neighbor_value /= n_neighbors[vertex];

                y_t[vertex] += weights[vertex] * step_size * (average_neighbor_value - y_t[vertex]);
            }

            if (mean_adjust) {
                // Adjusting y_t so that its mean = mean(y)
                Ey_t = 0;
                for (int vertex = 0; vertex < n_vertices; ++vertex)
                    Ey_t += y_t[vertex];
                Ey_t /= n_vertices;

                Delta = Ey - Ey_t;
                for (int vertex = 0; vertex < n_vertices; ++vertex)
                    y_t[vertex] += Delta;
            }

            itr++;
        }
    }

    y_traj->push_back(y_t);
    return y_traj;
}


/**
 * An R interface to the weighted_heat_smoother_trajectory function.
 *
 * This function takes an R list representing a graph adjacency list, a numeric vector
 * of function values at each vertex, a numeric vector of vertex weights, and various
 * parameters for the smoothing process. It calls the C++ weighted_heat_smoother_trajectory
 * function to smooth the function values using the weighted heat diffusion equation and
 * returns a list containing the trajectory of the smoothed function values over time,
 * along with the local minima and local maxima at each iteration.
 *
 * @param Rgraph An R list representing the graph adjacency list. Each element of the list
 *               should be an integer vector specifying the neighboring vertex indices for
 *               a given vertex.
 * @param Ry A numeric vector representing the initial function values at each vertex of the graph.
 * @param Rweights A numeric vector representing the weights of each vertex in the graph.
 *                 The length of Rweights must be equal to the number of vertices.
 * @param Rstep_factor A numeric scalar specifying the step factor used to calculate the step size
 *                     in the smoothing process. The step size is the product of the step factor
 *                     and the smallest difference between the unique values of the function.
 * @param Rprec A numeric scalar specifying the precision threshold used to determine local
 *              constancy. If the absolute difference between the value at a vertex and any
 *              of its neighbors is greater than this threshold, the vertex is not considered
 *              locally constant.
 * @param Rn_itrs An integer scalar specifying the maximum number of iterations for the smoothing
 *                process. This parameter is used only when the stopping criterion is set to 2.
 * @param Rstop_criterion An integer scalar specifying the stopping criterion for the smoothing process.
 *                        A value of 1 indicates that the process should stop when there are no locally
 *                        constant vertices remaining. A value of 2 indicates that the process should stop
 *                        after the specified number of iterations (Rn_itrs).
 * @param Rmean_adjust A boolean scalar specifying whether to adjust the mean of the smoothed function
 *                     values at each iteration to match the mean of the initial function values.
 *                     If Rmean_adjust is true, the mean adjustment is performed; otherwise, it is skipped.
 *
 * @return An R list containing three components:
 *         1. A matrix representing the trajectory of the smoothed function values over time, where each
 *            column corresponds to a specific iteration and each row represents the function values at
 *            a particular vertex.
 *         2. A list of integer vectors, where each vector represents the indices of the local minima
 *            at a specific iteration.
 *         3. A list of integer vectors, where each vector represents the indices of the local maxima
 *            at a specific iteration.
 */
SEXP S_weighted_heat_smoother_trajectory(SEXP Rgraph,
                                         SEXP Ry,
                                         SEXP Rweights,
                                         SEXP Rstep_factor,
                                         SEXP Rprec,
                                         SEXP Rn_itrs,
                                         SEXP Rstop_criterion,
                                         SEXP Rmean_adjust) {

    std::vector<std::vector<int>> graph = std::move(*Rlist_to_Cpp(Rgraph));
    int n_vertices = LENGTH(Ry);

    int nprot = 0;

    PROTECT(Ry = coerceVector(Ry, REALSXP)); nprot++;
    double *y = REAL(Ry);

    PROTECT(Rweights = coerceVector(Rweights, REALSXP)); nprot++;
    double *weights = REAL(Rweights);

    PROTECT(Rstep_factor = coerceVector(Rstep_factor, REALSXP)); nprot++;
    double step_factor = REAL(Rstep_factor)[0];

    PROTECT(Rprec = coerceVector(Rprec, REALSXP)); nprot++;
    double prec = REAL(Rprec)[0];

    PROTECT(Rn_itrs = coerceVector(Rn_itrs, INTSXP)); nprot++;
    int n_itrs = INTEGER(Rn_itrs)[0];

    PROTECT(Rstop_criterion = coerceVector(Rstop_criterion, INTSXP)); nprot++;
    int stop_criterion = INTEGER(Rstop_criterion)[0];

    PROTECT(Rmean_adjust = coerceVector(Rmean_adjust, LGLSXP)); nprot++;
    bool mean_adjust = LOGICAL(Rmean_adjust)[0];

    std::unique_ptr<std::vector<std::vector<double>>> y_traj = weighted_heat_smoother_trajectory(graph,
                                                                                                 std::vector<double>(y, y + n_vertices),
                                                                                                 std::vector<double>(weights, weights + n_vertices),
                                                                                                 step_factor,
                                                                                                 prec,
                                                                                                 n_itrs,
                                                                                                 stop_criterion,
                                                                                                 mean_adjust);

    int n_iterations = y_traj->size();

    SEXP y_traj_matrix = PROTECT(allocMatrix(REALSXP, n_vertices, n_iterations)); nprot++;
    double *y_traj_ptr = REAL(y_traj_matrix);

    for (int i = 0; i < n_iterations; i++) {
        std::copy((*y_traj)[i].begin(), (*y_traj)[i].end(), y_traj_ptr + i * n_vertices);
    }

    SEXP lmin_list = PROTECT(allocVector(VECSXP, n_iterations)); nprot++;
    SEXP lmax_list = PROTECT(allocVector(VECSXP, n_iterations)); nprot++;

    for (int i = 0; i < n_iterations; i++) {
        auto extrema = find_local_extrema(graph, (*y_traj)[i]);
        std::vector<int>& lmin = extrema->first;
        std::vector<int>& lmax = extrema->second;

        SEXP lmin_vec = PROTECT(allocVector(INTSXP, lmin.size()));
        int *lmin_ptr = INTEGER(lmin_vec);
        std::copy(lmin.begin(), lmin.end(), lmin_ptr);
        SET_VECTOR_ELT(lmin_list, i, lmin_vec);
        UNPROTECT(1);

        SEXP lmax_vec = PROTECT(allocVector(INTSXP, lmax.size()));
        int *lmax_ptr = INTEGER(lmax_vec);
        std::copy(lmax.begin(), lmax.end(), lmax_ptr);
        SET_VECTOR_ELT(lmax_list, i, lmax_vec);
        UNPROTECT(1);
    }

    SEXP result_list = PROTECT(allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(result_list, 0, y_traj_matrix);
    SET_VECTOR_ELT(result_list, 1, lmin_list);
    SET_VECTOR_ELT(result_list, 2, lmax_list);

    SEXP names = PROTECT(allocVector(STRSXP, 3)); ++nprot;
    SET_STRING_ELT(names, 0, mkChar("y_traj_matrix"));
    SET_STRING_ELT(names, 1, mkChar("lmin_list"));
    SET_STRING_ELT(names, 2, mkChar("lmax_list"));
    setAttrib(result_list, R_NamesSymbol, names);

    UNPROTECT(nprot);

    return result_list;
}
