#ifndef LM_H_
#define LM_H_

void C_flmw(const double *x,
            const double *y,
            const double *w,
            const int    *rnr,
            const int    *rnc,
                  double *beta,
                  int    *status);

void C_flmw_1D(const double *x,
               const double *y,
               const double *w,
               const int    *rn,
                     double *beta,
                     int    *status);

#endif // LM_H_
