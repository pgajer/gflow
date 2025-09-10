#ifndef MSR2_EIGEN_UTILS_H_
#define MSR2_EIGEN_UTILS_H_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>   // Rprintf, REprintf

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

// #include <omp.h>
#include "omp_compat.h"

// Undefine conflicting macros after including R headers
#undef length
#undef Rf_eval

//#include <iostream>
#include <string>

#include <Eigen/Dense>
#include <Eigen/Sparse>

SEXP EigenVectorXd_to_SEXP(const Eigen::VectorXd& vec);
SEXP EigenMatrixXd_to_SEXP(const Eigen::MatrixXd& mat);
SEXP EigenSparseMatrix_to_SEXP(const Eigen::SparseMatrix<double>& mat);

#endif // MSR2_EIGEN_UTILS_H_
