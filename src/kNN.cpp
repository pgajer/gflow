#include "kNN.h"

#include <vector>
#include <cmath>
#include <stdexcept>

#include <ANN/ANN.h>
#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_kNN(SEXP RX, SEXP Rk);
    SEXP S_kNN_v2(SEXP RX, SEXP Rk);
}

// ----------- C++ helper (RAII-only) -----------
knn_result_t kNN(const std::vector<std::vector<double>>& X, int k) {
    const int n = static_cast<int>(X.size());
    if (n <= 0) throw std::invalid_argument("kNN(): X has zero rows");
    const int d = static_cast<int>(X[0].size());
    if (d <= 0) throw std::invalid_argument("kNN(): X has zero columns");
    for (int i = 1; i < n; ++i) {
        if (static_cast<int>(X[i].size()) != d)
            throw std::invalid_argument("kNN(): inconsistent row lengths");
    }
    if (k <= 0)      throw std::invalid_argument("kNN(): k must be positive");
    if (k > n)       throw std::invalid_argument("kNN(): k cannot exceed n");

    // Build ANN dataset
    ANNpointArray dataPts = annAllocPts(n, d);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < d; ++j)
            dataPts[i][j] = X[i][j];

    ANNkd_tree* kdTree = nullptr;
    try {
        kdTree = new ANNkd_tree(dataPts, n, d);
    } catch (...) {
        annDeallocPts(dataPts);
        annClose();
        throw;
    }

    // Scratch buffers
    std::vector<ANNidx>  nnIdx(static_cast<size_t>(k));
    std::vector<ANNdist> nnDist(static_cast<size_t>(k));

    knn_result_t out;
    out.n = n; out.k = k;
    out.indices.resize(static_cast<size_t>(n) * k);
    out.distances.resize(static_cast<size_t>(n) * k);

    // Exact search
    for (int i = 0; i < n; ++i) {
        kdTree->annkSearch(dataPts[i], k, nnIdx.data(), nnDist.data(), 0.0);
        for (int j = 0; j < k; ++j) {
            const size_t off = static_cast<size_t>(i) * k + j;
            out.indices[off]   = nnIdx[j];
            out.distances[off] = std::sqrt(static_cast<double>(nnDist[j]));
        }
    }

    delete kdTree;
    annDeallocPts(dataPts);
    annClose();

    return out;
}

// ----------- R interfaces (no raw ownership crossing R boundary) -----------
SEXP S_kNN_v2(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = Rf_coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = Rf_coerceVector(Rk, INTSXP));  nprot++;

    if (!Rf_isMatrix(RX)) { UNPROTECT(nprot); Rf_error("S_kNN_v2: RX must be a matrix"); }

    const int k = INTEGER(Rk)[0];
    SEXP dim = Rf_getAttrib(RX, R_DimSymbol);
    const int n = INTEGER(dim)[0];
    const int d = INTEGER(dim)[1];
    if (n <= 0 || d <= 0) { UNPROTECT(nprot); Rf_error("S_kNN_v2: RX has non-positive dims"); }

    const double* X = REAL(RX);
    std::vector<std::vector<double>> Xv(n, std::vector<double>(d));
    for (int j = 0; j < d; ++j)
        for (int i = 0; i < n; ++i)
            Xv[i][j] = X[i + n * j];

    knn_result_t r;
    try { r = kNN(Xv, k); }
    catch (const std::exception& e) { UNPROTECT(nprot); Rf_error("S_kNN_v2: %s", e.what()); }

    SEXP nn_i = PROTECT(Rf_allocMatrix(INTSXP,   n, k)); nprot++;
    SEXP nn_d = PROTECT(Rf_allocMatrix(REALSXP,  n, k)); nprot++;
    for (int j = 0; j < k; ++j)
        for (int i = 0; i < n; ++i) {
            const size_t off = static_cast<size_t>(i) * k + j;
            INTEGER(nn_i)[i + n * j] = r.indices[off];
            REAL(nn_d)[i + n * j]    = r.distances[off];
        }

    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, nn_i);
    SET_VECTOR_ELT(res, 1, nn_d);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("indices"));
    SET_STRING_ELT(names, 1, Rf_mkChar("distances"));
    Rf_setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}

// Legacy direct-ANN version, also RAII-only.
SEXP S_kNN(SEXP RX, SEXP Rk) {
    int nprot = 0;
    PROTECT(RX = Rf_coerceVector(RX, REALSXP)); nprot++;
    PROTECT(Rk = Rf_coerceVector(Rk, INTSXP));  nprot++;

    if (!Rf_isMatrix(RX)) { UNPROTECT(nprot); Rf_error("S_kNN: RX must be a matrix"); }
    const int k = INTEGER(Rk)[0];

    SEXP dim = Rf_getAttrib(RX, R_DimSymbol);
    const int n = INTEGER(dim)[0];
    const int d = INTEGER(dim)[1];
    const double* X = REAL(RX);

    ANNpointArray dataPts = annAllocPts(n, d);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < d; ++j)
            dataPts[i][j] = X[i + n * j];

    ANNkd_tree* kdTree = nullptr;
    try { kdTree = new ANNkd_tree(dataPts, n, d); }
    catch (...) { annDeallocPts(dataPts); annClose(); UNPROTECT(nprot); Rf_error("S_kNN: kd-tree build failed"); }

    std::vector<ANNidx>  nnIdx(static_cast<size_t>(k));
    std::vector<ANNdist> nnDist(static_cast<size_t>(k));

    SEXP nn_i = PROTECT(Rf_allocMatrix(INTSXP,  n, k)); nprot++;
    SEXP nn_d = PROTECT(Rf_allocMatrix(REALSXP, n, k)); nprot++;

    for (int i = 0; i < n; ++i) {
        kdTree->annkSearch(dataPts[i], k, nnIdx.data(), nnDist.data(), 0.0);
        for (int j = 0; j < k; ++j) {
            INTEGER(nn_i)[i + n * j] = nnIdx[j];
            REAL(nn_d)[i + n * j]    = std::sqrt(static_cast<double>(nnDist[j]));
        }
    }

    delete kdTree;
    annDeallocPts(dataPts);
    annClose();

    SEXP res = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(res, 0, nn_i);
    SET_VECTOR_ELT(res, 1, nn_d);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(names, 0, Rf_mkChar("indices"));
    SET_STRING_ELT(names, 1, Rf_mkChar("distances"));
    Rf_setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(nprot);
    return res;
}
