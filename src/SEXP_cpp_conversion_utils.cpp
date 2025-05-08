// The helper functions follow these patterns:

// Sets are converted to R vectors
// Maps are converted to named R lists
// Nested structures (vector of sets, etc.) become nested R lists
// All numeric values are properly converted to R integer or real types

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
#include <unordered_map>

#include "msr2.h"
#include "cpp_utils.hpp"
#include "MS_complex.h"

/**
 * @brief Converts an R matrix to a C++ vector of vectors.
 *
 * This function takes an R matrix (SEXP) as input and converts it to a C++
 * std::vector<std::vector<double>>. It handles the column-major storage of R matrices
 * and returns a row-major C++ representation.
 *
 * @param Rmatrix SEXP representing an R matrix. Must be a numeric (double) matrix.
 *
 * @return std::unique_ptr<std::vector<std::vector<double>>> A unique pointer to a vector of vectors
 *         where each inner vector represents a row of the input matrix.
 *
 * @throws Rf_error If the input is not a numeric matrix.
 *
 * @note This function assumes column-major storage for the input R matrix and converts it
 *       to a row-major format in the returned C++ structure.
 *
 * @warning The function does not perform any data type conversion. It assumes that the input
 *          matrix contains double precision floating point numbers.
 *
 * Example usage:
 * @code
 * SEXP R_matrix = ...; // Assume this is a valid R matrix
 * try {
 *     auto cpp_matrix = Rmatrix_to_cpp(R_matrix);
 *     // Use cpp_matrix...
 * } catch (const std::exception& e) {
 *     REprintf("Error: %s\n", e.what());
 * }
 * @endcode
 *
 * @see R_list_of_dvectors_to_cpp_vector_of_dvectors for converting R list of vectors to C++
 */
std::unique_ptr<std::vector<std::vector<double>>> Rmatrix_to_cpp(SEXP Rmatrix) {
    // Check if the input is a matrix
    if (!isMatrix(Rmatrix) || !isReal(Rmatrix)) {
        Rf_error("Rmatrix_to_cpp: Input must be a numeric matrix");
    }

    // Get matrix dimensions
    SEXP Rdim = getAttrib(Rmatrix, R_DimSymbol);
    if (Rdim == R_NilValue || Rf_length(Rdim) != 2) {
        Rf_error("Rmatrix_to_cpp: Input must be a matrix");
    }

    int nrows = INTEGER(Rdim)[0];
    int ncols = INTEGER(Rdim)[1];

    // Get pointer to matrix data
    double* matrix_data = REAL(Rmatrix);

    // Create and populate the C++ vector of vectors
    auto cpp_matrix = std::make_unique<std::vector<std::vector<double>>>(nrows);

    for (int i = 0; i < nrows; ++i) {
        (*cpp_matrix)[i].reserve(ncols);
        for (int j = 0; j < ncols; ++j) {
            // R matrices are column-major, so we need to calculate the correct index
            (*cpp_matrix)[i].push_back(matrix_data[i + j * nrows]);
        }
    }

    return cpp_matrix;
}

/**
 * Converts an R numeric vector to a C++ vector of doubles.
 *
 * This function takes an R numeric vector (`Ry`) and converts it to a C++ vector
 * (`std::vector<double>`). The function iterates over each element of the R vector,
 * extracts the double value, and assigns it to the corresponding element of the
 * C++ vector.
 *
 * Caller Responsibility: The caller should ensure that `Ry` is protected before
 * passing it to `Rvect_to_CppVect_double`.
 *
 * Function Assumption: The `Rvect_to_CppVect_double` function assumes that the passed
 * `Ry` is protected and does not need to protect it again.
 *
 * Protection Management: The caller should manage the protection lifecycle of `Ry`,
 * ensuring it remains protected for the duration of its use.
 *
 * @param Ry An R numeric vector.
 *
 * @return A unique pointer to a C++ vector of doubles (`std::unique_ptr<std::vector<double>>`)
 *         representing the values in the R numeric vector.
 */
std::unique_ptr<std::vector<double>> Rvect_to_CppVect_double(SEXP Ry) {
    std::vector<double> y(REAL(Ry), REAL(Ry) + LENGTH(Ry));
    return std::make_unique<std::vector<double>>(std::move(y));
}

#if 0
std::unique_ptr<std::vector<double>> Rvect_to_CppVect_double(SEXP Ry) {
    int y_length = LENGTH(Ry);
    auto y = std::make_unique<std::vector<double>>(y_length);

    for (int i = 0; i < y_length; ++i)
        (*y)[i] = REAL(Ry)[i];

    return y;
}
#endif


/**
 * @brief Converts an R list of numeric vectors to a C++ vector of vectors of doubles.
 *
 * This function takes an R object (SEXP) representing a list of double vectors and
 * converts it to a C++ std::vector<std::vector<double>>.
 *
 * @param Rvectvect An R object (SEXP) representing a list of double vectors.
 * @return A unique pointer to a vector of vectors of doubles.
 *
 * @note The function assumes that Rvectvect is a valid R list of double vectors.
 *       It does not perform extensive error checking on the input.
 */
std::unique_ptr<std::vector<std::vector<double>>> R_list_of_dvectors_to_cpp_vector_of_dvectors(SEXP Rvectvect) {
    int n_vertices = LENGTH(Rvectvect);
    auto cpp_vect_list = std::make_unique<std::vector<std::vector<double>>>(n_vertices);

    for (int i = 0; i < n_vertices; ++i) {
        SEXP R_vector = VECTOR_ELT(Rvectvect, i);
        int vector_length = LENGTH(R_vector);
        double* double_vector = REAL(R_vector);

        (*cpp_vect_list)[i].reserve(vector_length);
        (*cpp_vect_list)[i].assign(double_vector, double_vector + vector_length);
    }

    return cpp_vect_list;
}

/**
 * @brief Converts an R adjacency list to a C++ vector representation
 *
 * @details This function takes an R SEXP object representing an adjacency list
 * and converts it to a C++ nested vector representation. The function expects
 * the input to be a list of integer vectors, where each vector contains the
 * indices of adjacent vertices. The function performs input validation to ensure
 * the proper format.
 *
 * @param s_adj_list An R SEXP object representing an adjacency list
 * @return std::vector<std::vector<int>> A C++ nested vector representation of the adjacency list
 * @throws R error if the input is not a list or contains non-integer vectors
 *
 * @note The function assumes indices are already adjusted for C++ (0-based indexing)
 */
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list) {
    // Input validation
    if (!isNewList(s_adj_list)) {
        error("Expected a list for adjacency list");
    }

    int n_vertices = LENGTH(s_adj_list);
    std::vector<std::vector<int>> adj_list(n_vertices);

    for (int i = 0; i < n_vertices; i++) {
        SEXP vertex_neighbors = VECTOR_ELT(s_adj_list, i);

        if (!isInteger(vertex_neighbors)) {
            error("Expected integer vector for vertex neighbors");
        }

        int* neighbors = INTEGER(vertex_neighbors);
        int n_neighbors = LENGTH(vertex_neighbors);

        // Efficient bulk assignment - no index conversion
        adj_list[i].assign(neighbors, neighbors + n_neighbors);
    }
    return adj_list;
}

/**
 * @brief Converts an R edge weight list to a C++ vector representation
 *
 * @details This function takes an R SEXP object representing a list of edge weights
 * and converts it to a C++ nested vector representation. The function expects
 * the input to be a list of numeric vectors, where each vector contains the
 * weights of edges corresponding to the adjacency list. The function performs
 * input validation to ensure the proper format.
 *
 * @param s_weight_list An R SEXP object representing a list of edge weights
 * @return std::vector<std::vector<double>> A C++ nested vector representation of the weight list
 * @throws R error if the input is not a list or contains non-numeric vectors
 *
 * @note The function assumes the weight list structure matches the corresponding adjacency list
 */
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list) {
    // Input validation
    if (!isNewList(s_weight_list)) {
        error("Expected a list for weight list");
    }

    int n_vertices = LENGTH(s_weight_list);
    std::vector<std::vector<double>> weight_list(n_vertices);

    for (int i = 0; i < n_vertices; i++) {
        SEXP vertex_weights = VECTOR_ELT(s_weight_list, i);

        if (!isReal(vertex_weights)) {
            error("Expected numeric vector for vertex weights");
        }

        double* weights = REAL(vertex_weights);
        int n_weights = LENGTH(vertex_weights);

        // Efficient bulk assignment - no index conversion
        weight_list[i].assign(weights, weights + n_weights);
    }
    return weight_list;
}


/**
 * Converts a C++ vector of vectors of integers to an R list of integer vectors.
 *
 * This function takes a C++ std::vector<std::vector<int>> and converts it to an R list
 * of integer vectors. Each inner vector of the input is converted to an R integer vector,
 * and these vectors are combined into an R list.
 *
 * @param cpp_vec_vec The input C++ vector of vectors of integers to be converted.
 * @return An R list of integer vectors representing the converted input.
 */
SEXP convert_vector_vector_int_to_R(const std::vector<std::vector<int>>& cpp_vec_vec) {

    int n = cpp_vec_vec.size();

    SEXP Rlist = PROTECT(allocVector(VECSXP, n));

    for (int i = 0; i < n; ++i) {
        int m = cpp_vec_vec[i].size();

        SEXP Rvec = PROTECT(allocVector(INTSXP, m));
        int* ptr = INTEGER(Rvec);

        for (int j = 0; j < m; ++j) {
            ptr[j] = cpp_vec_vec[i][j];
        }

        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return Rlist;
}

/**
 * @brief Converts a C++ vector of vectors of doubles to an R list of numeric vectors
 *
 * @param vec Input vector of vectors of doubles to convert
 * @return SEXP A protected R list where each element is a numeric vector
 *             corresponding to the inner vectors of the input
 *
 * @details This function takes a C++ nested vector structure containing double values
 *          and converts it to an R list where each element is a numeric vector.
 *          The function handles all necessary R object protection and unprotection.
 *          The double precision values are preserved in the conversion process.
 *
 * @note The returned SEXP object is protected once and should be unprotected by the caller
 *       if necessary
 *
 * @warning The function assumes the input vector is valid and non-null. Empty vectors
 *          are handled correctly but not specially treated.
 *
 * @see REAL, PROTECT, UNPROTECT, allocVector, SET_VECTOR_ELT
 */
SEXP convert_vector_vector_double_to_R(const std::vector<std::vector<double>>& vec) {
    int n = vec.size();

    SEXP Rlist = PROTECT(allocVector(VECSXP, n));

    for (int i = 0; i < n; ++i) {
        int m = vec[i].size();

        SEXP Rvec = PROTECT(allocVector(REALSXP, m));
        double* ptr = REAL(Rvec);

        for (int j = 0; j < m; ++j) {
            ptr[j] = vec[i][j];
        }

        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return Rlist;
}


/**
 * @brief Converts a C++ vector of vectors to an R matrix
 *
 * @details
 * This function handles the conversion from C++ std::vector<std::vector<double>>
 * to an R matrix (SEXP), accounting for:
 * - R's column-major order storage
 * - Memory protection
 * - Empty input handling
 *
 * The function performs the following steps:
 * 1. Checks for empty input
 * 2. Allocates R matrix
 * 3. Copies data with transposition for column-major order
 * 4. Handles memory protection
 *
 * @param data Input vector of vectors where:
 *             - Outer vector represents rows
 *             - Inner vectors represent columns
 *             - All inner vectors must have the same size
 *             - data[i][j] represents element at row i, column j
 *
 * @return SEXP (REALSXP matrix) where:
 *         - Returns R_NilValue if input is empty
 *         - Returns nrow × ncol matrix otherwise
 *         - Matrix is column-major ordered as required by R
 *         - Matrix elements are copied from input with proper transposition
 *
 * @note Memory Management:
 *       - Uses PROTECT/UNPROTECT for R object safety
 *       - Caller does not need to UNPROTECT the result
 *       - Properly handles cleanup on error
 *
 * @warning
 *       - Assumes all inner vectors have the same length
 *       - Does not verify inner vector size consistency
 *       - May cause undefined behavior if inner vectors differ in size
 *
 * @see convert_vector_double_to_R For single vector conversion
 * @see convert_vector_int_to_R For integer vector conversion
 *
 * @example
 * ```cpp
 * std::vector<std::vector<double>> data = {{1.0, 2.0}, {3.0, 4.0}};
 * SEXP matrix = convert_vector_vector_double_to_matrix(data);
 * // Results in 2×2 R matrix:
 * // [,1] [,2]
 * // [1,]  1.0  2.0
 * // [2,]  3.0  4.0
 * ```
 */
SEXP convert_vector_vector_double_to_matrix(const std::vector<std::vector<double>>& data) {
    if (data.empty()) return R_NilValue;

    int nrow = data.size();
    int ncol = data[0].size();

    SEXP matrix = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    double* matrix_ptr = REAL(matrix);

    // R matrices are column-major, so we need to transpose during copying
    for (int j = 0; j < ncol; ++j) {
        for (int i = 0; i < nrow; ++i) {
            matrix_ptr[i + j * nrow] = data[i][j];
        }
    }

    UNPROTECT(1);
    return matrix;
}


/**
 * @brief Converts a C++ vector of vectors of booleans to an R list of logical vectors
 *
 * @param vec Input vector of vectors of booleans to convert
 * @return SEXP A protected R list where each element is a logical vector
 *             corresponding to the inner vectors of the input
 *
 * @details This function takes a C++ nested vector structure containing boolean values
 *          and converts it to an R list where each element is a logical vector.
 *          The function handles all necessary R object protection and unprotection.
 *          Each inner vector becomes a logical vector in R, preserving the boolean values.
 *
 * @note The returned SEXP object is protected once and should be unprotected by the caller
 *       if necessary
 *
 * @warning The function assumes the input vector is valid and non-null. Empty vectors
 *          are handled correctly but not specially treated.
 */
SEXP convert_vector_vector_bool_to_R(const std::vector<std::vector<bool>>& vec) {
    int n = vec.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    for (int i = 0; i < n; ++i) {
        int m = vec[i].size();
        SEXP Rvec = PROTECT(allocVector(LGLSXP, m));
        int* ptr = LOGICAL(Rvec);
        for (int j = 0; j < m; ++j) {
            ptr[j] = vec[i][j];
        }
        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return Rlist;
}



/**
 * Converts a C++ vector of doubles to an R numeric vector.
 *
 * This function takes a C++ std::vector<double> and converts it to an R numeric vector (SEXP).
 * The function allocates a new R numeric vector with the same length as the input C++ vector
 * and copies the values from the C++ vector to the R vector.
 *
 * @param cpp_vec The input C++ vector of doubles to be converted.
 * @return An R numeric vector (SEXP) containing the values from the input C++ vector.
 *
 * @note The returned R vector is protected from garbage collection using the PROTECT macro,
 *       and the protection is removed before returning the vector using the UNPROTECT macro.
 *
 * Example usage:
 * @code
 * std::vector<double> cpp_vec = {1.5, 2.3, 3.7, 4.2, 5.1};
 * SEXP Rvec = convert_vector_double_to_R(cpp_vec);
 * // Rvec now contains the values {1.5, 2.3, 3.7, 4.2, 5.1} as an R numeric vector
 * @endcode
 */
SEXP convert_vector_double_to_R(const std::vector<double>& vec) {
    int n = vec.size();

    SEXP Rvec = PROTECT(allocVector(REALSXP, n));
    double* ptr = REAL(Rvec);

    for (int j = 0; j < n; ++j)
        ptr[j] = vec[j];

    UNPROTECT(1);
    return Rvec;
}

/**
 * @brief Converts a C++ vector of integers to an R integer vector
 *
 * @param vec Input vector of integers to convert
 * @return SEXP A protected R integer vector containing the same values as the input
 *
 * @details This function takes a C++ vector of integers and converts it to an R integer
 *          vector, preserving all values. The function handles the necessary R object
 *          protection and unprotection.
 *
 * @note The returned SEXP object is protected once and should be unprotected by the caller
 *       if necessary
 *
 * @warning The function assumes the input vector is valid. No range checking is performed
 *          on the integer values, so values outside the range of R's integers may cause
 *          undefined behavior
 *
 * @see INTEGER, PROTECT, UNPROTECT, allocVector
 */
SEXP convert_vector_int_to_R(const std::vector<int>& vec) {
    int n = vec.size();

    SEXP Rvec = PROTECT(allocVector(INTSXP, n));
    int* ptr = INTEGER(Rvec);

    for (int j = 0; j < n; ++j)
        ptr[j] = vec[j];

    UNPROTECT(1);
    return Rvec;
}

/**
 * @brief Converts a C++ vector of booleans to an R logical vector
 *
 * @param vec Input vector of booleans to convert
 * @return SEXP A protected R logical vector containing the same values as the input
 *
 * @details This function takes a C++ vector of booleans and converts it to an R logical
 *          vector, preserving the boolean values. The function handles the necessary R
 *          object protection and unprotection. Note that R's logical values are stored
 *          as integers internally (0 for FALSE, 1 for TRUE).
 *
 * @note The returned SEXP object is protected once and should be unprotected by the caller
 *       if necessary
 *
 * @warning Special consideration should be given to std::vector<bool> which is a specialized
 *          template that packs booleans into bits for space efficiency
 *
 * @see LOGICAL, PROTECT, UNPROTECT, allocVector
 */
SEXP convert_vector_bool_to_R(const std::vector<bool>& vec) {
    int n = vec.size();
    SEXP Rvec = PROTECT(allocVector(LGLSXP, n));
    int* ptr = LOGICAL(Rvec);
    for (int j = 0; j < n; ++j)
        ptr[j] = vec[j];
    UNPROTECT(1);
    return Rvec;
}

/**
 * @brief Converts a C++ unordered map to an R list of integer vectors
 *
 * @param cpp_map_int_vect_int An unordered_map from integers to vectors of integers
 * @return SEXP A protected R list where each element is an integer vector
 *
 * @details The function creates an R list where:
 *          - Each map key determines the position in the list
 *          - Each map value becomes an integer vector
 *          - List length equals (maximum key + 1)
 *          - Missing keys result in NULL elements
 *
 * @example
 * Input: {{1, {1,2}}, {3, {4,5}}}
 * Result: List of length 4:
 *   [[1]] = c(1,2)
 *   [[2]] = NULL
 *   [[3]] = c(4,5)
 *   [[4]] = NULL
 */
SEXP convert_map_int_vector_int_to_R(const std::unordered_map<int, std::vector<int>>& cpp_map_int_vect_int) {

    // Identifying the largest key of the map
    int max_key = 0;
    for (const auto& pair : cpp_map_int_vect_int) {
        if (pair.first > max_key)
            max_key = pair.first;
    }

    int Rlist_len = max_key + 1;
    SEXP Rlist = PROTECT(allocVector(VECSXP, Rlist_len));

    for (const auto& [key, vect] : cpp_map_int_vect_int) {
        int m = vect.size();
        SEXP Rvec = PROTECT(allocVector(INTSXP, m));
        int* ptr = INTEGER(Rvec);

        for (int j = 0; j < m; ++j)
            ptr[j] = vect[j];

        SET_VECTOR_ELT(Rlist, key, Rvec);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return Rlist;
}

/**
 * @brief Converts a C++ unordered map to an R named list of integer vectors
 *
 * @param cpp_map_int_vect_int An unordered_map from integers to vectors of integers
 * @param names Vector of integers to use as component names
 * @return SEXP A protected R list where each element is an integer vector with named components
 *
 * @details The function creates an R list where:
 *          - Each map key determines the position in the list
 *          - Each map value becomes an integer vector
 *          - List length equals (maximum key + 1)
 *          - Missing keys result in NULL elements
 *          - Names vector must match the list length
 *
 * @throws R error if names vector length doesn't match list length
 *
 * @example
 * Input map: {{1, {1,2}}, {3, {4,5}}}
 * Names: {10, 20, 30, 40}
 * Result: List of length 4:
 *   [[1]] = c(1,2)      name: "10"
 *   [[2]] = NULL        name: "20"
 *   [[3]] = c(4,5)      name: "30"
 *   [[4]] = NULL        name: "40"
 */
SEXP convert_map_int_vector_int_to_R(const std::unordered_map<int, std::vector<int>>& cpp_map_int_vect_int,
                                     const std::vector<int>& names) {
    // Identifying the largest key of the map
    int max_key = 0;
    for (const auto& pair : cpp_map_int_vect_int) {
        if (pair.first > max_key)
            max_key = pair.first;
    }

    int Rlist_len = max_key + 1;
    if (Rlist_len != (int)names.size()) {
        error("Rlist_len != names.size(): Rlist_len: %d\tnames.size(): %d\n\n", Rlist_len, (int)names.size());
    }

    SEXP Rlist = PROTECT(allocVector(VECSXP, Rlist_len));

    for (const auto& [key, vect] : cpp_map_int_vect_int) {
        int m = vect.size();
        SEXP Rvec = PROTECT(allocVector(INTSXP, m));
        int* ptr = INTEGER(Rvec);

        for (int j = 0; j < m; ++j)
            ptr[j] = vect[j];

        SET_VECTOR_ELT(Rlist, key, Rvec);
        UNPROTECT(1);
    }

    // Convert the names vector to an R character vector
    SEXP Rnames = PROTECT(allocVector(STRSXP, names.size()));
    for (size_t i = 0; i < names.size(); ++i) {
        SET_STRING_ELT(Rnames, i, mkChar(std::to_string(names[i]).c_str()));
    }

    // Set the names of the list elements
    setAttrib(Rlist, R_NamesSymbol, Rnames);

    UNPROTECT(2);
    return Rlist;
}

/**
 * Converts a C++ unordered map of integers to integer sets into an R list of integer vectors.
 *
 * This function takes a C++ std::unordered_map<int, std::set<int>> and converts
 * it to an R list of integer vectors. Each inner set of the input is converted
 * to an R integer vector, and these vectors are combined into an R list.
 *
 * @param cpp_map_int_set_int The input C++ unordered map from integers to sets of integers to be converted.
 * @return An R list of integer vectors representing the converted input.
 *
 * @note The function assumes that the keys of the input map are non-negative integers.
 * The resulting R list will have a length equal to the maximum key value plus one,
 * and the elements corresponding to missing keys will be NULL.
 *
 * @note The function uses 1-based indexing for the R list, consistent with R's indexing convention.
 * @note The function protects the allocated R objects to prevent memory leaks.
 */
SEXP Cpp_map_int_set_int_to_Rlist(const std::unordered_map<int, std::set<int>>& cpp_map_int_set_int) {
    // Identifying the largest key of the map
    int max_key = 0;
    for (const auto& pair : cpp_map_int_set_int) {
        if (pair.first > max_key)
            max_key = pair.first;
    }

    int Rlist_len = max_key + 1;

    SEXP Rlist = PROTECT(allocVector(VECSXP, Rlist_len));

    for (const auto& [key, set] : cpp_map_int_set_int) {
        int m = set.size();
        SEXP Rvec = PROTECT(allocVector(INTSXP, m));
        int* ptr = INTEGER(Rvec);
        int j = 0;
        for (const auto& e : set)
            ptr[j++] = e;

        if (key < 0 || key > Rlist_len - 1) {
            error("In Cpp_map_int_set_int_to_Rlist(): key: %d is out of range[0,Rlist_len - 1]\n", key);
        }

        SET_VECTOR_ELT(Rlist, key, Rvec);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return Rlist;
}


/**
 * @brief Converts a C++ vector of integer pairs to an R matrix (unique pointer version)
 *
 * This function takes a unique pointer to a vector of integer pairs and
 * converts it to an R matrix with two columns. Each pair in the input vector
 * becomes a row in the output matrix.
 *
 * @param cpp_vector A unique pointer to a vector of integer pairs to be converted
 * @return SEXP An R matrix (INTSXP) with two columns, containing the data from the input vector
 *
 * @note The function uses PROTECT/UNPROTECT to manage R's garbage collection.
 *       Caller does not need to UNPROTECT the returned SEXP.
 *
 * @warning The function assumes that the input vector is not empty. Behavior is
 *          undefined for empty vectors.
 */
SEXP uptr_vector_of_pairs_to_R_matrix(const std::unique_ptr<std::vector<std::pair<int, int>>>& cpp_vector) {
    int n_rows = cpp_vector->size();
    SEXP R_matrix = PROTECT(allocMatrix(INTSXP, n_rows, 2));
    int* matrix_data = INTEGER(R_matrix);
    for (int i = 0; i < n_rows; ++i) {
        matrix_data[i] = (*cpp_vector)[i].first + 1;           // + 1 to change it to 1-based integers
        matrix_data[i + n_rows] = (*cpp_vector)[i].second + 1;
    }
    UNPROTECT(1);
    return R_matrix;
}

/**
 * @brief Converts a C++ vector of integer pairs to an R matrix
 *
 * This function takes a vector of integer pairs and converts it to an R matrix
 * with two columns. Each pair in the input vector becomes a row in the output matrix.
 *
 * @param cpp_vector A const reference to a vector of integer pairs to be converted
 * @return SEXP An R matrix (INTSXP) with two columns, containing the data from the input vector
 *
 * @note The function uses PROTECT/UNPROTECT to manage R's garbage collection.
 *       Caller does not need to UNPROTECT the returned SEXP.
 *
 * @warning The function assumes that the input vector is not empty. Behavior is
 *          undefined for empty vectors.
 */
SEXP cpp_vector_of_pairs_to_R_matrix(const std::vector<std::pair<int, int>>& cpp_vector) {
    int n_rows = cpp_vector.size();
    SEXP R_matrix = PROTECT(allocMatrix(INTSXP, n_rows, 2));
    int* matrix_data = INTEGER(R_matrix);
    for (int i = 0; i < n_rows; ++i) {
        matrix_data[i] = cpp_vector[i].first + 1;           // + 1 to change it to 1-based integers
        matrix_data[i + n_rows] = cpp_vector[i].second + 1;
    }
    UNPROTECT(1);
    return R_matrix;
}


/**
 * @brief Converts a flat C++ vector to an R matrix.
 *
 * This function takes a flat C++ vector representing a matrix in row-major order
 * and converts it into an R matrix object (SEXP). The function is particularly
 * useful for interfacing C++ functions that return flat matrices with R, where
 * matrix objects are more commonly used.
 *
 * @param flat_matrix A const reference to a std::vector<double> containing the
 *                    matrix elements in row-major order.
 * @param nrow The number of rows in the matrix.
 * @param ncol The number of columns in the matrix.
 *
 * @return SEXP A PROTECT'd R matrix object (REALSXP) containing the data from
 *              the input flat_matrix. The caller is responsible for calling
 *              UNPROTECT(1) after using the returned object.
 *
 * @note This function allocates memory for a new R matrix and copies the data
 *       from the input vector. The input vector is not modified.
 *
 * @warning The function assumes that the size of flat_matrix is equal to
 *          nrow * ncol. No bounds checking is performed, so ensure the input
 *          parameters are consistent to avoid undefined behavior.
 *
 * @see S_shortest_path for an example of how this function is used in
 *      conjunction with R interface functions.
 *
 * Example usage:
 * @code
 * std::vector<double> flat_mat = {1.0, 2.0, 3.0, 4.0};
 * SEXP r_mat = PROTECT(flat_vector_to_R_matrix(flat_mat, 2, 2));
 * // Use r_mat...
 * UNPROTECT(1);
 * @endcode
 */
SEXP flat_vector_to_R_matrix(const std::vector<double>& flat_matrix, int nrow, int ncol) {
    SEXP r_matrix = PROTECT(allocMatrix(REALSXP, nrow, ncol));
    double* matrix_ptr = REAL(r_matrix);

    for (int i = 0; i < nrow * ncol; ++i) {
        matrix_ptr[i] = flat_matrix[i];
    }

    UNPROTECT(1);
    return r_matrix;
}


/**
 * @brief Converts a C++ set of integers to an R integer vector
 *
 * @param set Input set of integers to be converted
 * @return SEXP R integer vector containing all elements from the input set
 *
 * @note The returned SEXP is protected and unprotected within the function
 * @note The order of elements in the resulting R vector follows the natural ordering of the set
 */
SEXP convert_set_to_R(const std::set<int>& set) {
    SEXP Rvec = PROTECT(allocVector(INTSXP, set.size()));
    int* ptr = INTEGER(Rvec);
    int i = 0;
    for (const auto& val : set) {
        ptr[i++] = val;
    }
    UNPROTECT(1);
    return Rvec;
}


/**
 * @brief Converts a C++ map from integers to sets of integers to a named R list
 *
 * @param map_set Input map where keys are integers and values are sets of integers
 * @return SEXP Named R list where names are the keys from the map (as strings)
 *              and values are integer vectors (converted from sets)
 *
 * @note The returned SEXP is protected and unprotected within the function
 * @note List names are created by converting integer keys to strings
 */
SEXP convert_map_set_to_R(const std::map<int, std::set<int>>& map_set) {
    int n = map_set.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    SEXP names = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    for (const auto& pair : map_set) {
        // Convert the set to R vector
        SET_VECTOR_ELT(Rlist, i, convert_set_to_R(pair.second));
        // Set the name as the key
        SET_STRING_ELT(names, i, mkChar(std::to_string(pair.first).c_str()));
        i++;
    }

    setAttrib(Rlist, R_NamesSymbol, names);
    UNPROTECT(2);
    return Rlist;
}

/**
 * @brief Converts a C++ map from integer pairs to vectors of integer sets to a named R list
 *
 * @param map_vec_set Input map where keys are pairs of integers and values are vectors of integer sets
 * @return SEXP Named R list where:
 *              - Names are created by concatenating the pair of keys with underscore ("key1_key2")
 *              - Each value is a list of integer vectors (converted from vector of sets)
 *
 * @note The returned SEXP is protected and unprotected within the function
 * @note Creates a nested list structure to represent the vector of sets
 */
SEXP convert_map_vector_set_to_R(const std::map<std::pair<int,int>, std::vector<std::set<int>>>& map_vec_set) {
    int n = map_vec_set.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    SEXP names = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    for (const auto& pair : map_vec_set) {
        // Create nested list for vector of sets
        SEXP inner_list = PROTECT(allocVector(VECSXP, pair.second.size()));
        for (size_t j = 0; j < pair.second.size(); ++j) {
            SET_VECTOR_ELT(inner_list, j, convert_set_to_R(pair.second[j]));
        }

        // Set the name as "key1_key2"
        std::string name = std::to_string(pair.first.first) + "_" + std::to_string(pair.first.second);
        SET_STRING_ELT(names, i, mkChar(name.c_str()));

        SET_VECTOR_ELT(Rlist, i, inner_list);
        UNPROTECT(1); // inner_list
        i++;
    }

    setAttrib(Rlist, R_NamesSymbol, names);
    UNPROTECT(2);
    return Rlist;
}

/**
 * @brief Converts a C++ map of cell trajectories to a named R list
 *
 * @param cell_traj Input map where keys are cell_t structures and values are sets of size_t indices
 * @return SEXP Named R list where:
 *              - Names are created by concatenating cell components ("lmax_lmin_cell_index")
 *              - Values are integer vectors containing trajectory indices
 *
 * @note The returned SEXP is protected and unprotected within the function
 * @note Converts size_t values to integers for R compatibility
 */
SEXP convert_cell_trajectories_to_R(const std::map<cell_t, std::set<size_t>>& cell_traj) {
    int n = cell_traj.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    SEXP names = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    for (const auto& pair : cell_traj) {
        // Convert set of size_t to vector of integers
        SEXP Rvec = PROTECT(allocVector(INTSXP, pair.second.size()));
        int* ptr = INTEGER(Rvec);
        int j = 0;
        for (const auto& val : pair.second) {
            ptr[j++] = static_cast<int>(val);
        }

        // Create name from cell_t components
        std::string name = std::to_string(pair.first.lmax) + "_" +
                          std::to_string(pair.first.lmin) + "_" +
                          std::to_string(pair.first.cell_index);
        SET_STRING_ELT(names, i, mkChar(name.c_str()));

        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1); // Rvec
        i++;
    }

    setAttrib(Rlist, R_NamesSymbol, names);
    UNPROTECT(2);
    return Rlist;
}


/**
 * @brief Converts a C++ map from integer pairs to vectors of integers to a named R list
 *
 * @param map_vec Input map where keys are pairs of integers and values are vectors of integers
 * @return SEXP Named R list where:
 *              - Names are created by concatenating the pair of keys with underscore ("key1_key2")
 *              - Values are integer vectors
 *
 * @note The returned SEXP is protected and unprotected within the function
 */
SEXP convert_map_vector_to_R(const std::map<std::pair<int,int>, std::vector<int>>& map_vec) {
    int n = map_vec.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    SEXP names = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    for (const auto& pair : map_vec) {
        // Convert vector to R vector
        SEXP Rvec = PROTECT(allocVector(INTSXP, pair.second.size()));
        int* ptr = INTEGER(Rvec);
        for (size_t j = 0; j < pair.second.size(); ++j) {
            ptr[j] = pair.second[j];
        }

        // Set the name as "key1_key2"
        std::string name = std::to_string(pair.first.first) + "_" + std::to_string(pair.first.second);
        SET_STRING_ELT(names, i, mkChar(name.c_str()));

        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1); // Rvec
        i++;
    }

    setAttrib(Rlist, R_NamesSymbol, names);
    UNPROTECT(2);
    return Rlist;
}

/**
 * @brief Converts a named R list to a C++ map from integer pairs to vectors of integers
 *
 * This is an inverse of convert_map_vector_to_R
 *
 * @param Rlist Input named R list where:
 *              - Names are strings in the format "key1_key2" where key1 and key2 are integers
 *              - Values are integer vectors
 * @return std::map<std::pair<int,int>, std::vector<int>> Map where:
 *              - Keys are pairs of integers extracted from list names
 *              - Values are vectors of integers from the list elements
 *
 * @throws std::runtime_error If list names are not in the expected "key1_key2" format
 *                           or if list elements are not integer vectors
 * @note The function handles R object protection internally
 */
std::map<std::pair<int,int>, std::vector<int>> convert_R_to_map_vector(SEXP Rlist) {
    if (!isNewList(Rlist)) {
        Rf_error("Input must be an R list");
    }

    std::map<std::pair<int,int>, std::vector<int>> result;
    int n = Rf_length(Rlist);

    // Get names
    SEXP names = PROTECT(getAttrib(Rlist, R_NamesSymbol));
    if (names == R_NilValue) {
        UNPROTECT(1);
        Rf_error("Input list must have names");
    }

    // Process each element
    for (int i = 0; i < n; i++) {
        // Parse name to get key pair
        std::string name = CHAR(STRING_ELT(names, i));
        size_t underscore_pos = name.find('_');
        if (underscore_pos == std::string::npos) {
            UNPROTECT(1);
            Rf_error("List names must be in format 'key1_key2'");
        }

        // Extract and convert keys
        try {
            int key1 = std::stoi(name.substr(0, underscore_pos));
            int key2 = std::stoi(name.substr(underscore_pos + 1));
            std::pair<int,int> key(key1, key2);

            // Get vector element
            SEXP Rvec = VECTOR_ELT(Rlist, i);
            if (!isInteger(Rvec)) {
                UNPROTECT(1);
                Rf_error("List elements must be integer vectors");
            }

            // Convert R vector to std::vector
            int* ptr = INTEGER(Rvec);
            int vec_length = Rf_length(Rvec);
            std::vector<int> value(ptr, ptr + vec_length);

            // Insert into map
            result[key] = value;
        }
        catch (const std::invalid_argument& e) {
            UNPROTECT(1);
            Rf_error("Failed to parse integer keys from list names");
        }
    }

    UNPROTECT(1);
    return result;
}


/**
 * @brief Converts a C++ map from integer pairs to sets of integers (procells) to a named R list
 *
 * @param procells Input map where keys are pairs of integers (typically max-min pairs)
 *                 and values are sets of integers representing proto-cells
 * @return SEXP Named R list where:
 *              - Names are created by concatenating the pair of keys with underscore ("max_min")
 *              - Values are integer vectors (converted from sets)
 *
 * @note The returned SEXP is protected and unprotected within the function
 * @note List names follow the format "max_min" where max and min are the components of the key pair
 */
SEXP convert_procells_to_R(const std::map<std::pair<int,int>, std::set<int>>& procells) {
    int n = procells.size();
    SEXP Rlist = PROTECT(allocVector(VECSXP, n));
    SEXP names = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    for (const auto& pair : procells) {
        // Convert the set to R vector
        SEXP Rvec = PROTECT(allocVector(INTSXP, pair.second.size()));
        int* ptr = INTEGER(Rvec);
        int j = 0;
        for (const auto& val : pair.second) {
            ptr[j++] = val;
        }

        // Set the name as "max_min"
        std::string name = std::to_string(pair.first.first) + "_" + std::to_string(pair.first.second);
        SET_STRING_ELT(names, i, mkChar(name.c_str()));

        SET_VECTOR_ELT(Rlist, i, Rvec);
        UNPROTECT(1); // Rvec
        i++;
    }

    setAttrib(Rlist, R_NamesSymbol, names);
    UNPROTECT(2);
    return Rlist;
}


/**
 * @brief Converts an R list containing shortest paths data into a C++ map structure
 *
 * @details The function takes an R list with three named components:
 *          - 'i': Integer vector of source vertices (1-based indices)
 *          - 'j': Integer vector of target vertices (1-based indices)
 *          - 'paths': List of integer vectors, each containing a path from i to j (1-based indices)
 *          The function converts these R data structures back into a C++ map where keys are
 *          pairs of vertices (i,j) and values are vectors containing the paths between them.
 *          All indices are converted from R's 1-based to C++'s 0-based indexing.
 *
 * @param s_shortest_paths SEXP object containing the R list with shortest paths data
 *
 * @return std::map with pairs of vertices as keys and vectors of path vertices as values
 *
 * @note The input R list should have the exact structure as produced by S_create_path_graph
 *
 * @throw No explicit exceptions, but R error handling may trigger if input structure is invalid
 */
std::map<std::pair<int,int>, std::vector<int>> shortest_paths_Rlist_to_cpp_map(SEXP s_shortest_paths) {
    std::map<std::pair<int,int>, std::vector<int>> result;

    // Extract the three components from the list
    SEXP i_coords = VECTOR_ELT(s_shortest_paths, 0);  // First element (i coordinates)
    SEXP j_coords = VECTOR_ELT(s_shortest_paths, 1);  // Second element (j coordinates)
    SEXP paths = VECTOR_ELT(s_shortest_paths, 2);     // Third element (paths list)

    // Get pointers to the coordinate arrays
    int* i_ptr = INTEGER(i_coords);
    int* j_ptr = INTEGER(j_coords);

    // Get the length of the arrays (they should all be the same length)
    R_xlen_t n_paths = XLENGTH(i_coords);

    // Iterate through all paths
    for (R_xlen_t idx = 0; idx < n_paths; ++idx) {
        // Get current i,j coordinates (subtract 1 to convert from R's 1-based to C++'s 0-based indexing)
        int i = i_ptr[idx] - 1;
        int j = j_ptr[idx] - 1;

        // Get the current path vector
        SEXP current_path = VECTOR_ELT(paths, idx);
        int* path_ptr = INTEGER(current_path);
        R_xlen_t path_length = XLENGTH(current_path);

        // Create vector for this path
        std::vector<int> path_vec;
        path_vec.reserve(path_length);

        // Copy path values (subtract 1 to convert from R's 1-based to C++'s 0-based indexing)
        for (R_xlen_t k = 0; k < path_length; ++k) {
            path_vec.push_back(path_ptr[k] - 1);
        }

        // Insert into map
        result[std::make_pair(i, j)] = std::move(path_vec);
    }

    return result;
}




/**
 * @brief Converts a 3D vector of doubles to an R list of matrices
 *
 * @param data 3D vector organized as [outer_idx][middle_idx][inner_idx]
 *
 * @return SEXP (List) A list of matrices where:
 *         - Each list element corresponds to outer_idx
 *         - Each matrix has rows corresponding to middle_idx
 *         - Each matrix has columns corresponding to inner_idx
 *
 * @note Memory Management:
 *       - Allocates new R objects
 *       - Caller must handle PROTECT/UNPROTECT
 *       - Returns PROTECTED object
 *
 * @warning
 *       - Assumes non-empty input vector
 *       - Assumes consistent dimensions for all inner vectors
 *       - Caller must UNPROTECT returned object
 */
SEXP convert_vector_vector_vector_double_to_R(
    const std::vector<std::vector<std::vector<double>>>& data) {

    if (data.empty()) return R_NilValue;

    int n_outer = data.size();
    SEXP result = PROTECT(allocVector(VECSXP, n_outer));

    for (int i = 0; i < n_outer; i++) {
        if (data[i].empty()) {
            SET_VECTOR_ELT(result, i, R_NilValue);
            continue;
        }

        int n_rows = data[i].size();
        int n_cols = data[i][0].size();

        // Create matrix for this outer index
        SEXP matrix = PROTECT(allocMatrix(REALSXP, n_rows, n_cols));
        double* ptr = REAL(matrix);

        // Fill matrix
        for (int j = 0; j < n_rows; j++) {
            for (int k = 0; k < n_cols; k++) {
                ptr[j + k * n_rows] = data[i][j][k];  // Column-major order for R
            }
        }

        SET_VECTOR_ELT(result, i, matrix);
        UNPROTECT(1);  // Unprotect matrix after setting it in list
    }

    UNPROTECT(1);

    return result;  // Still protected, caller must UNPROTECT
}
