#ifndef MSR2_STATS_UTILS_H_
#define MSR2_STATS_UTILS_H_

double mean(const double *x, int n );
double wmean(const double *x, const double *w, int n);
double median(double *a, int n );
double* runif_hcube(int n, int dim, double* L, double w);
int dcmp( void const *e1, void const *e2 );

#endif // MSR2_STATS_UTILS_H_
