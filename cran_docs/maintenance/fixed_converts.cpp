/**
 * @file fixed_converts.cpp
 * @brief Refactored conversion functions with proper PROTECT/UNPROTECT handling
 * 
 * This file contains refactored versions of conversion functions from
 * SEXP_cpp_conversion_utils.cpp to fix PROTECT/UNPROTECT issues identified by rchk.
 * 
 * Key principles applied:
 * A) One disciplined counter per function
 * B) Protect immediately; do not let a fresh SEXP cross an allocating call
 * C) Avoid "unsupported form of unprotect"
 * D) REPROTECT when replacing a protected object
 * E) Lists and attributes - allocate container first, then set attributes/elements
 * F) Strings/names - Rf_mkChar allocates but is safe when target is protected
 * G) Loops - UNPROTECT per iteration if PROTECT per iteration
 * H) Early returns and errors - ensure balanced PROTECT/UNPROTECT
 */

#include <vector>
#include <unordered_map>
#include <set>
#include <R.h>
#include <Rinternals.h>

// Assume set_wgraph_t is defined elsewhere with appropriate structure
// For compilation, you'll need to include the proper header

/**
 * @brief Converts a C++ vector of doubles to an R numeric vector
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_double_to_R(const std::vector<double>& vec) {
    int nprot = 0;
    int n = vec.size();
    
    SEXP Rvec = PROTECT(Rf_allocVector(REALSXP, n)); ++nprot;
    double* ptr = REAL(Rvec);
    
    for (int j = 0; j < n; ++j) {
        ptr[j] = vec[j];
    }
    
    UNPROTECT(nprot);
    return Rvec;
}

/**
 * @brief Converts a C++ vector of integers to an R integer vector
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_int_to_R(const std::vector<int>& vec) {
    int nprot = 0;
    int n = vec.size();
    
    SEXP Rvec = PROTECT(Rf_allocVector(INTSXP, n)); ++nprot;
    int* ptr = INTEGER(Rvec);
    
    for (int j = 0; j < n; ++j) {
        ptr[j] = vec[j];
    }
    
    UNPROTECT(nprot);
    return Rvec;
}

/**
 * @brief Converts a C++ vector of booleans to an R logical vector
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_bool_to_R(const std::vector<bool>& vec) {
    int nprot = 0;
    int n = vec.size();
    
    SEXP Rvec = PROTECT(Rf_allocVector(LGLSXP, n)); ++nprot;
    int* ptr = LOGICAL(Rvec);
    
    for (int j = 0; j < n; ++j) {
        ptr[j] = vec[j] ? TRUE : FALSE;
    }
    
    UNPROTECT(nprot);
    return Rvec;
}

/**
 * @brief Converts a C++ vector of vectors of doubles to an R list
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>& vec) {
    int nprot = 0;
    int n = vec.size();
    
    // Allocate and protect the container first
    SEXP Rlist = PROTECT(Rf_allocVector(VECSXP, n)); ++nprot;
    
    for (int i = 0; i < n; ++i) {
        int m = vec[i].size();
        
        // Allocate and protect the element
        SEXP Rvec = PROTECT(Rf_allocVector(REALSXP, m)); ++nprot;
        double* ptr = REAL(Rvec);
        
        // Fill the element
        for (int j = 0; j < m; ++j) {
            ptr[j] = vec[i][j];
        }
        
        // Set element in container (container now holds reference)
        SET_VECTOR_ELT(Rlist, i, Rvec);
        
        // Unprotect the element (safe because container holds it)
        UNPROTECT(1); --nprot;
    }
    
    UNPROTECT(nprot);
    return Rlist;
}

/**
 * @brief Converts a C++ vector of vectors of integers to an R list
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>& cpp_vec_vec) {
    int nprot = 0;
    int n = cpp_vec_vec.size();
    
    // Allocate and protect the container first
    SEXP Rlist = PROTECT(Rf_allocVector(VECSXP, n)); ++nprot;
    
    for (int i = 0; i < n; ++i) {
        int m = cpp_vec_vec[i].size();
        
        // Allocate and protect the element
        SEXP Rvec = PROTECT(Rf_allocVector(INTSXP, m)); ++nprot;
        int* ptr = INTEGER(Rvec);
        
        // Fill the element
        for (int j = 0; j < m; ++j) {
            ptr[j] = cpp_vec_vec[i][j];
        }
        
        // Set element in container
        SET_VECTOR_ELT(Rlist, i, Rvec);
        
        // Unprotect the element
        UNPROTECT(1); --nprot;
    }
    
    UNPROTECT(nprot);
    return Rlist;
}

/**
 * @brief Converts a C++ vector of vectors of booleans to an R list
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_vector_bool_to_R(const std::vector<std::vector<bool>>& vec) {
    int nprot = 0;
    int n = vec.size();
    
    // Allocate and protect the container first
    SEXP Rlist = PROTECT(Rf_allocVector(VECSXP, n)); ++nprot;
    
    for (int i = 0; i < n; ++i) {
        int m = vec[i].size();
        
        // Allocate and protect the element
        SEXP Rvec = PROTECT(Rf_allocVector(LGLSXP, m)); ++nprot;
        int* ptr = LOGICAL(Rvec);
        
        // Fill the element
        for (int j = 0; j < m; ++j) {
            ptr[j] = vec[i][j] ? TRUE : FALSE;
        }
        
        // Set element in container
        SET_VECTOR_ELT(Rlist, i, Rvec);
        
        // Unprotect the element
        UNPROTECT(1); --nprot;
    }
    
    UNPROTECT(nprot);
    return Rlist;
}

/**
 * @brief Converts a C++ vector of vectors to an R matrix
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_vector_vector_double_to_matrix(const std::vector<std::vector<double>>& data) {
    // Handle empty input
    if (data.empty()) {
        return R_NilValue;
    }
    
    int nprot = 0;
    int nrow = data.size();
    int ncol = data[0].size();
    
    // Allocate and protect the matrix
    SEXP matrix = PROTECT(Rf_allocMatrix(REALSXP, nrow, ncol)); ++nprot;
    double* matrix_ptr = REAL(matrix);
    
    // R matrices are column-major, so we need to transpose during copying
    for (int j = 0; j < ncol; ++j) {
        for (int i = 0; i < nrow; ++i) {
            matrix_ptr[i + j * nrow] = data[i][j];
        }
    }
    
    UNPROTECT(nprot);
    return matrix;
}

/**
 * @brief Converts a C++ unordered map to an R named list
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_map_int_vector_int_to_R(
    const std::unordered_map<int, std::vector<int>>& cpp_map_int_vect_int,
    const std::vector<int>& names) {
    
    int nprot = 0;
    
    // First pass: validate and find max key (no PROTECTs yet)
    if (cpp_map_int_vect_int.empty() && names.empty()) {
        // Return an empty list
        SEXP empty = PROTECT(Rf_allocVector(VECSXP, 0)); ++nprot;
        UNPROTECT(nprot);
        return empty;
    }
    
    int max_key = -1;
    for (const auto& pair : cpp_map_int_vect_int) {
        if (pair.first > max_key) {
            max_key = pair.first;
        }
    }
    
    int Rlist_len = max_key + 1;
    if (Rlist_len != (int)names.size()) {
        Rf_error("convert_map_int_vector_int_to_R: Rlist_len != names.size(): "
                 "Rlist_len: %d\tnames.size(): %d", Rlist_len, (int)names.size());
    }
    
    // Allocate and protect the container first
    SEXP Rlist = PROTECT(Rf_allocVector(VECSXP, Rlist_len)); ++nprot;
    
    // Allocate and protect the names vector
    SEXP Rnames = PROTECT(Rf_allocVector(STRSXP, names.size())); ++nprot;
    
    // Fill the names vector (safe because Rnames is protected)
    for (size_t i = 0; i < names.size(); ++i) {
        SET_STRING_ELT(Rnames, i, Rf_mkChar(std::to_string(names[i]).c_str()));
    }
    
    // Set the names attribute (safe because both are protected)
    Rf_setAttrib(Rlist, R_NamesSymbol, Rnames);
    
    // Now fill the list elements
    for (const auto& [key, vect] : cpp_map_int_vect_int) {
        int m = vect.size();
        
        // Allocate and protect the element
        SEXP Rvec = PROTECT(Rf_allocVector(INTSXP, m)); ++nprot;
        int* ptr = INTEGER(Rvec);
        
        // Fill the element
        for (int j = 0; j < m; ++j) {
            ptr[j] = vect[j];
        }
        
        // Set element in container
        SET_VECTOR_ELT(Rlist, key, Rvec);
        
        // Unprotect the element
        UNPROTECT(1); --nprot;
    }
    
    UNPROTECT(nprot);
    return Rlist;
}

// For convert_wgraph_to_R, we need the actual definition of set_wgraph_t
// Here's a stub that shows the proper PROTECT/UNPROTECT pattern

struct edge_t {
    size_t vertex;
    double weight;
};

struct set_wgraph_t {
    std::vector<std::set<edge_t>> adjacency_list;
    
    size_t num_vertices() const {
        return adjacency_list.size();
    }
};

/**
 * @brief Converts a weighted graph to R representation
 * 
 * Fixed version with proper PROTECT/UNPROTECT handling.
 * The returned SEXP is NOT protected - caller must protect if needed.
 */
SEXP convert_wgraph_to_R(const set_wgraph_t& graph) {
    int nprot = 0;
    size_t n_vertices = graph.num_vertices();
    
    // Allocate and protect containers first
    SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); ++nprot;
    SEXP weight_list = PROTECT(Rf_allocVector(VECSXP, n_vertices)); ++nprot;
    
    // Process each vertex
    for (size_t i = 0; i < n_vertices; i++) {
        size_t n_neighbors = graph.adjacency_list[i].size();
        
        // Allocate and protect elements for this vertex
        SEXP RA = PROTECT(Rf_allocVector(INTSXP, n_neighbors));
        SEXP RW = PROTECT(Rf_allocVector(REALSXP, n_neighbors));
        
        int* A = INTEGER(RA);
        double* W = REAL(RW);
        
        // Fill the elements
        size_t j = 0;
        for (const auto& edge : graph.adjacency_list[i]) {
            A[j] = static_cast<int>(edge.vertex + 1);  // Convert to 1-based indexing
            W[j] = edge.weight;
            j++;
        }
        
        // Set elements in containers (containers now hold references)
        SET_VECTOR_ELT(adj_list, i, RA);
        SET_VECTOR_ELT(weight_list, i, RW);
        
        // Unprotect the elements (safe because containers hold them)
        UNPROTECT(2); // RA, RW
    }
    
    // Create the result list
    SEXP r_list = PROTECT(Rf_allocVector(VECSXP, 2)); ++nprot;
    SET_VECTOR_ELT(r_list, 0, adj_list);
    SET_VECTOR_ELT(r_list, 1, weight_list);
    
    // Create and set names
    SEXP r_list_names = PROTECT(Rf_allocVector(STRSXP, 2)); ++nprot;
    SET_STRING_ELT(r_list_names, 0, Rf_mkChar("adj_list"));
    SET_STRING_ELT(r_list_names, 1, Rf_mkChar("weight_list"));
    Rf_setAttrib(r_list, R_NamesSymbol, r_list_names);
    
    UNPROTECT(nprot);
    return r_list;
}

/* Additional notes on the refactoring:
 * 
 * 1. Each function now has a single nprot counter that tracks protections
 * 2. All PROTECTs are followed immediately by ++nprot
 * 3. All UNPROTECTs use --nprot when in loops
 * 4. Final UNPROTECT(nprot) at single exit point
 * 5. Container allocated and protected before elements
 * 6. Elements unprotected after being placed in container
 * 7. Names/attributes set while container is protected
 * 8. No unprotected SEXPs cross allocating calls
 * 9. Early returns (like empty case) properly balance PROTECT/UNPROTECT
 * 
 * These patterns ensure:
 * - No stack imbalance
 * - No GC hazards
 * - No "unsupported form of unprotect" errors
 * - Proper cleanup on all paths
 */
