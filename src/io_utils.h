#ifndef GFLOW_IO_UTILS_H_
#define GFLOW_IO_UTILS_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

FILE* file_open(const char *file, const char *mode);
void file_close(FILE *fp, const char *file);

void print_double_array(const double *x, int n);
void print_double_array_with_precision(const double *x, int n, int precision);
void write_double_array(const double *x, int n, const char *file);

void print_int_array(const int *x, int n);
void write_int_array(const int *x, int n, const char *file);

void print_2d_double_array_first_n(const double *x, int nr, int nc, int N, int precision);
void print_2d_double_array(const double *x, int nr, int nc);
void write_2d_double_array(const double *x, int nr, int nc, const char *file);
void write_2d_int_array(const int *x, int nr, int nc, const char *file);

#ifdef __cplusplus
}
#endif

#endif // GFLOW_IO_UTILS_H_
