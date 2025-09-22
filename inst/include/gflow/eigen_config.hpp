#pragma once

// Hard block any accidental TBB usage from Eigen
#ifdef EIGEN_USE_TBB
#  error "EIGEN_USE_TBB must not be defined for CRAN. Remove any TBB linkage."
#endif

// If OpenMP is not active, force Eigen single-thread mode
#ifndef _OPENMP
#  ifndef EIGEN_DONT_PARALLELIZE
#    define EIGEN_DONT_PARALLELIZE
#  endif
#endif
