#include <R.h>
#include <Rinternals.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

// in R run: .Call("S_gflow_openmp_diag", PACKAGE = "gflow")

extern "C" SEXP S_gflow_openmp_diag() {
  SEXP ans = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, Rf_mkChar("openmp_compiled"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("_OPENMP"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("omp_get_max_threads"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("omp_get_num_procs"));
  Rf_setAttrib(ans, R_NamesSymbol, nms);

#ifdef _OPENMP
  SET_VECTOR_ELT(ans, 0, Rf_ScalarLogical(1));
  SET_VECTOR_ELT(ans, 1, Rf_ScalarInteger((int)_OPENMP));
  SET_VECTOR_ELT(ans, 2, Rf_ScalarInteger((int)omp_get_max_threads()));
  SET_VECTOR_ELT(ans, 3, Rf_ScalarInteger((int)omp_get_num_procs()));
#else
  SET_VECTOR_ELT(ans, 0, Rf_ScalarLogical(0));
  SET_VECTOR_ELT(ans, 1, Rf_ScalarInteger(0));
  SET_VECTOR_ELT(ans, 2, Rf_ScalarInteger(1));
  SET_VECTOR_ELT(ans, 3, Rf_ScalarInteger(1));
#endif

  UNPROTECT(2);
  return ans;
}
