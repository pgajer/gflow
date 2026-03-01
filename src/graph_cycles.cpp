#include "graph_cycles_r.h"
#include "cpp_utils.hpp"

#include <vector>
#include <queue>
#include <memory>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <stack>
#include <numeric>

#include <R.h>
#include <Rinternals.h>

// ------------------------------------------------------------------------------------------
//
// Cycle Detection DFS Recursive Algorithm
//
// ------------------------------------------------------------------------------------------

// (s <- sort(cycle.sizes(graph.adj.list.5)))
// + table(s)
//  [1] 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 5
// > s
//  3  4  5
// 10  5  1

// The number of 3-cycles in K_5 is C(5,3) = 4*5/2 = 10;
// The number of 4-cycles in K_5 is C(5,4) = 5

#if 0

/**
 * @brief Finds cycles in a graph using depth-first search (DFS).
 *
 * This function recursively explores the graph, starting from a given vertex.
 * When a path revisits its starting vertex (and the path length is sufficient),
 * a cycle is detected.
 *
 * @param vertex The current vertex being visited.
 * @param start_vertex The starting vertex of the current DFS path.
 * @param parent_vertex The parent of the current vertex in the DFS tree.
 * @param graph The graph represented as an adjacency list.
 * @param visited A vector indicating whether a vertex has been visited in the current DFS path.
 * @param path A vector storing the vertices of the current DFS path.
 * @param candidate_cycles A set to store unique cycles found during the search.
 *
 * @return None. Cycles are directly stored within the 'candidate_cycles' parameter.
 */
void dfs_find_cycles(int vertex,
                     int start_vertex,
                     int parent_vertex,
                     const std::vector<std::vector<int>>& graph,
                     std::vector<bool>& visited, std::vector<int>& path,
                     std::set<std::vector<int>>& candidate_cycles) {
    visited[vertex] = true;
    path.push_back(vertex);

    for (int neighbor : graph[vertex]) {
        if (neighbor == start_vertex && path.size() > 2) {
            // Found a cycle
            std::vector<int> cycle(path.begin(), path.end());
            // Sort the cycle to avoid duplicates
            std::sort(cycle.begin(), cycle.end());
            candidate_cycles.insert(cycle);
        } else if (!visited[neighbor] && neighbor != parent_vertex) {
            dfs_find_cycles(neighbor, start_vertex, vertex, graph, visited, path, candidate_cycles);
        }
    }

    visited[vertex] = false;
    path.pop_back();
}

/**
 * @brief Finds cycles within a graph and returns their lengths.
 *
 * This function iterates over vertices in the graph, initiating a DFS from each vertex.
 * Detected cycles are recorded, and their lengths are returned.
 *
 * @param graph The graph represented as an adjacency list.
 * @return A unique pointer to a vector containing the lengths of unique cycles in the graph.
 */
std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph) {
    int num_vertices = graph.size();
    std::vector<bool> visited(num_vertices, false);
    std::vector<int> path;
    auto cycle_lengths = std::make_unique<std::vector<int>>();
    std::set<std::vector<int>> candidate_cycles;

#if 0
    int min_deg = 1;
    // finding the number of vertices of degree equal or greater than min_deg
    int n_min_deg_vertices = 0;
    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        int deg = graph[vertex].size(); // degree of the vertex
        if ( deg >= min_deg ) n_min_deg_vertices++;
    }
    Rprintf("In cycle_sizes\n");
    Rprintf("n_min_deg_vertices: %d\n", n_min_deg_vertices);
    int counter = 0;
#endif

    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        // int deg = graph[vertex].size(); // degree of the vertex
        // if ( deg < min_deg ) continue;
        // Rprintf("%d / %d\n", counter, n_min_deg_vertices);
        // counter++;
        dfs_find_cycles(vertex, vertex, -1, graph, visited, path, candidate_cycles);
    }

    // Rprintf("candidate_cycles\n");
    // for(const auto& cycle : candidate_cycles)
    //     print_vector(cycle);

    // Copy cycles from set to vector
    std::vector<std::vector<int>> cycles_vector(candidate_cycles.begin(), candidate_cycles.end());

    // Record cycle lengths
    for (const auto& cycle : cycles_vector) {
        cycle_lengths->push_back(cycle.size());
    }

    return cycle_lengths;
}

#endif


#if 1
//
// [claude3 iterative version]
//
void dfs_find_cycles(int vertex, int start_vertex, const std::vector<std::vector<int>>& graph,
                     std::vector<bool>& visited, std::vector<int>& path,
                     std::set<std::vector<int>>& candidate_cycles) {
    std::stack<std::pair<int, int>> stack;
    stack.push({vertex, -1});
    visited[vertex] = true;

    while (!stack.empty()) {
        int current_vertex = stack.top().first;
        //int parent_vertex = stack.top().second;
        stack.pop();

        path.push_back(current_vertex);

        for (int neighbor : graph[current_vertex]) {
            if (neighbor == start_vertex && path.size() > 2) {
                // Found a cycle
                std::vector<int> cycle(path.begin(), path.end());
                // Sort the cycle to avoid duplicates
                std::sort(cycle.begin(), cycle.end());
                candidate_cycles.insert(cycle);
            }

            if (!visited[neighbor]) {
                visited[neighbor] = true;
                stack.push({neighbor, current_vertex});
            }
        }

        if (stack.empty() || stack.top().second != current_vertex) {
            // Backtrack
            visited[current_vertex] = false;
            path.pop_back();
        }
    }
}

std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph) {
    int num_vertices = graph.size();
    std::vector<bool> visited(num_vertices, false);
    std::vector<int> path;
    auto cycle_lengths = std::make_unique<std::vector<int>>();
    std::set<std::vector<int>> candidate_cycles;

    int min_deg = 3;
    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        int deg = graph[vertex].size(); // degree of the vertex
        if (deg < min_deg) continue;
        dfs_find_cycles(vertex, vertex, graph, visited, path, candidate_cycles);
    }

    // Copy cycles from set to vector
    std::vector<std::vector<int>> cycles_vector(candidate_cycles.begin(), candidate_cycles.end());

    // Record cycle lengths
    for (const auto& cycle : cycles_vector) {
        cycle_lengths->push_back(cycle.size());
    }

    return cycle_lengths;
}
#endif


#if 0
//
// [Gemini iterative version]
//
void dfs_find_cycles(const std::vector<std::vector<int>>& graph, int start_vertex,
                     std::set<std::vector<int>>& candidate_cycles) {
  std::vector<bool> visited(graph.size(), false);
  std::stack<std::pair<int, int>> stack;
  std::vector<int> path;

  stack.push({start_vertex, -1});

  while (!stack.empty()) {
    int vertex = stack.top().first;
    //int parent_vertex = stack.top().second;
    stack.pop();

    if (visited[vertex]) {
      // ... (cycle handling remains the same) ...
    } else {
      visited[vertex] = true;
      path.push_back(vertex);

      for (auto it = graph[vertex].rbegin(); it != graph[vertex].rend(); ++it) {
         // ... (neighbor logic remains the same) ...
      }

      // Selective unmarking upon backtracking
      path.pop_back();
      visited[vertex] = false;
    }
  }
}

std::unique_ptr<std::vector<int>> cycle_sizes(const std::vector<std::vector<int>>& graph) {
  auto cycle_lengths = std::make_unique<std::vector<int>>();
  std::set<std::vector<int>> candidate_cycles;

  int min_deg = 3;
  for (int vertex = 0; vertex < graph.size(); ++vertex) {
    int deg = graph[vertex].size(); // degree of the vertex
    if (deg < min_deg) continue;
    dfs_find_cycles(graph, vertex, candidate_cycles);
  }

  // Record cycle lengths
  for (const auto& cycle : candidate_cycles) {
    cycle_lengths->push_back(cycle.size());
  }

  return cycle_lengths;
}
#endif


/**
 * @brief R interface for the cycle_sizes function.
 *
 * This function takes an R graph adjacency list as input, converts it to a C++ representation,
 * calls the cycle_sizes function to find the lengths of the cycles in the graph, and returns
 * the result as an R integer vector.
 *
 * @param RA An R graph adjacency list representing the graph. It should be a list of integer
 *           vectors, where each vector contains the neighboring vertex indices for a given vertex.
 *           The vertex indices should be 1-based.
 *
 * @return An R integer vector containing the lengths of the cycles found in the graph.
 *
 * @throw An Rf_error is thrown if the input is not a valid R list or if any element of the list is
 *        not an integer vector.
 */
SEXP S_cycle_sizes(SEXP RA) {
    // Check if the input is a list
    if (!Rf_isNewList(RA))
        Rf_error("Input must be a list");

    // Get the length of the list
    int n = LENGTH(RA);

    // Convert R adjacency list to std::vector<std::vector<int>>
    std::vector<std::vector<int>> graph;
    for (int i = 0; i < n; ++i) {
        SEXP vec = VECTOR_ELT(RA, i);
        if (!Rf_isInteger(vec))
            Rf_error("Each element of the list must be an integer vector");

        int* ptr = INTEGER(vec);
        int len = LENGTH(vec);
        std::vector<int> neighbors(ptr, ptr + len);
        graph.push_back(neighbors);
    }

    // Call the cycle_sizes function
    std::unique_ptr<std::vector<int>> cycle_lengths = cycle_sizes(graph);

    // Convert the output to SEXP
    SEXP result = PROTECT(Rf_allocVector(INTSXP, cycle_lengths->size()));
    int* result_ptr = INTEGER(result);
    std::copy(cycle_lengths->begin(), cycle_lengths->end(), result_ptr);
    UNPROTECT(1);

    return result;
}
