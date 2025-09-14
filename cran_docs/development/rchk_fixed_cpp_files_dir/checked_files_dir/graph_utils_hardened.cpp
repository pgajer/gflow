// graph_utils_hardened.cpp — rchk-hardened drop-in for S_join_graphs
//
// Changes vs. provided graph_utils.cpp:
// - Replace PROTECT_WITH_INDEX coercions with Rf_asInteger() (no allocation to track)
// - Add strict argument validation (types, lengths, NA checks)
// - Add bounds checks for i1, i2 against graph sizes
// - Container-first allocation & single PROTECT for the only newly-allocated SEXP (result)
// - Literal UNPROTECT(1) — no variable counters
//
// This file is intended as a drop-in replacement for the S_join_graphs() entry point.

#include <vector>
#include <memory>
#include <stdexcept>
#include <R.h>
#include <Rinternals.h>

// Forward declarations for helper functions (assumed available in the build)
extern std::vector<std::vector<int>> convert_adj_list_from_R(SEXP);
extern SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>&);

// Core function declaration (assumed available)
std::unique_ptr<std::vector<std::vector<int>>> join_graphs(
    const std::vector<std::vector<int>>& graph1,
    const std::vector<std::vector<int>>& graph2,
    int i1, int i2);

extern "C" {

SEXP S_join_graphs(SEXP Rgraph1, SEXP Rgraph2, SEXP Ri1, SEXP Ri2) {
    // Basic type sanity checks (fast-fail, defensive)
    if (TYPEOF(Rgraph1) != VECSXP)
        Rf_error("S_join_graphs(): Rgraph1 must be a list (adjacency list).");
    if (TYPEOF(Rgraph2) != VECSXP)
        Rf_error("S_join_graphs(): Rgraph2 must be a list (adjacency list).");

    // Convert adjacency lists from R to C++ containers
    std::vector<std::vector<int>> graph1 = convert_adj_list_from_R(Rgraph1);
    std::vector<std::vector<int>> graph2 = convert_adj_list_from_R(Rgraph2);

    // Indices — use Rf_asInteger (no need to PROTECT, result is a plain int)
    int i1 = Rf_asInteger(Ri1);
    int i2 = Rf_asInteger(Ri2);
    if (i1 == NA_INTEGER || i2 == NA_INTEGER) {
        Rf_error("S_join_graphs(): i1 and i2 must be valid integers (not NA).");
    }
    if (i1 < 0 || i2 < 0) {
        Rf_error("S_join_graphs(): i1 and i2 must be non-negative.");
    }

    // Bounds checks vs. graph sizes
    const R_xlen_t n1 = static_cast<R_xlen_t>(graph1.size());
    const R_xlen_t n2 = static_cast<R_xlen_t>(graph2.size());
    if (n1 == 0 || n2 == 0) {
        Rf_error("S_join_graphs(): adjacency lists cannot be empty.");
    }
    if (static_cast<R_xlen_t>(i1) >= n1) {
        Rf_error("S_join_graphs(): i1=%d out of range [0, %lld).", i1, (long long)n1);
    }
    if (static_cast<R_xlen_t>(i2) >= n2) {
        Rf_error("S_join_graphs(): i2=%d out of range [0, %lld).", i2, (long long)n2);
    }

    // Core operation
    std::unique_ptr<std::vector<std::vector<int>>> joined_graph =
        join_graphs(graph1, graph2, i1, i2);
    if (!joined_graph) {
        Rf_error("S_join_graphs(): join_graphs() returned null.");
    }

    // Convert back to R. Only new allocation that needs protection.
    SEXP result = PROTECT(convert_vector_vector_int_to_R(*joined_graph));

    UNPROTECT(1);
    return result;
}

} // extern "C"
