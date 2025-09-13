#ifndef CV_DEG0_H_
#define CV_DEG0_H_

void C_cv_deg0_binloss(const int    *rnfolds,
                       const int    *rnreps,
                       const int    *rnNN,
                       const int    *rybinary,
                       const int    *nn_i,
                       const double *nn_d,
                       const double *y,
                       const int    *rK,
                       const int    *rng,
                       const int    *rnrX,
                       const double *rbw,
                       const int    *rminK,
                       const int    *rikernel,
                             double *rbinloss);

void C_cv_deg0_mae(const int    *rnfolds,
                   const int    *rnreps,
                   const int    *rnNN,
                   const int    *rybinary,
                   const int    *nn_i,
                   const double *nn_d,
                   const double *y,
                   const int    *rK,
                   const int    *rng,
                   const int    *rnrX,
                   const double *rbw,
                   const int    *rminK,
                   const int    *rikernel,
				         double *rMAE);
#endif // CV_DEG0_H_
