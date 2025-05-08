#ifndef MSR2_EIGEN_UTILS_H_
#define MSR2_EIGEN_UTILS_H_

#include <R.h>
#include <Rinternals.h>

// Prevent macro collision with OpenMP
#ifdef match
#undef match
#endif

#include <omp.h>

// Undefine conflicting macros after including R headers
#undef length

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>

SEXP EigenVectorXd_to_SEXP(const Eigen::VectorXd& vec);
SEXP EigenMatrixXd_to_SEXP(const Eigen::MatrixXd& mat);
SEXP EigenSparseMatrix_to_SEXP(const Eigen::SparseMatrix<double>& mat);

#endif // MSR2_EIGEN_UTILS_H_
