/*
  grid routines
*/

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "msr2.h"

/*!
  \brief Creates a equidistant 2D grid with each edge of the grid having the same length

  \param rdx   A reference to length of each edge.
  \param rx1L  A reference to the left end of the x1 range.
  \param rn1   A reference to number of grid elements along the first coordinate.
  \param rx2L  A reference to the left end of the x2 range.
  \param rn2   A reference to number of grid elements along the second coordinate.
  \param grid  The output grid.
*/
void C_create_ED_grid_2D(const double *rdx,
                         const double *rx1L,
                         const int    *rn1,
                         const double *rx2L,
                         const int    *rn2,
                         double *grid)
{
    double dx  = rdx[0];
    double x1L = rx1L[0];
    int n1     = rn1[0];
    double x2L = rx2L[0];
    int n2     = rn2[0];

    double* x1 = calloc(n1, sizeof(double));
    CHECK_PTR(x1);

    double* x2 = calloc(n2, sizeof(double));
    CHECK_PTR(x2);

    x1[0] = x1L;
    x2[0] = x2L;

    for ( int i = 1; i < n1; i++ )
      x1[i] = x1[i-1] + dx;

    for ( int i = 1; i < n2; i++ )
      x2[i] = x2[i-1] + dx;


    for ( int i1 = 0; i1 < n1; i1++ )
    {
      for ( int i2 = 0; i2 < n2; i2++ )
      {
        grid[    (i1 + i2*n1)*2] = x1[i1];
        grid[1 + (i1 + i2*n1)*2] = x2[i2];
      }
    }

    free(x1);
    free(x2);
}

/*!
  \brief Creates a equidistant 3D grid with each edge of the grid having the same length

  \param rdx   A reference to length of each edge.
  \param rx1L  A reference to the left end of the x1 range.
  \param rn1   A reference to number of grid elements along the first coordinate.
  \param rx2L  A reference to the left end of the x2 range.
  \param rn2   A reference to number of grid elements along the second coordinate.
  \param rx3L  A reference to the left end of the x3 range.
  \param rn3   A reference to number of grid elements along the third coordinate.
  \param grid  The output grid.
*/
void C_create_ED_grid_3D(const double *rdx,
                         const double *rx1L,
                         const int    *rn1,
                         const double *rx2L,
                         const int    *rn2,
                         const double *rx3L,
                         const int    *rn3,
                         double *grid)
{
    double dx  = rdx[0];

    double x1L = rx1L[0];
    int n1     = rn1[0];
    double x2L = rx2L[0];
    int n2     = rn2[0];
    double x3L = rx3L[0];
    int n3     = rn3[0];

    double* x1 = calloc(n1, sizeof(double));
    CHECK_PTR(x1);

    double* x2 = calloc(n2, sizeof(double));
    CHECK_PTR(x2);

    double* x3 = calloc(n3, sizeof(double));
    CHECK_PTR(x3);

    x1[0] = x1L;
    x2[0] = x2L;
    x3[0] = x3L;

    for ( int i = 1; i < n1; i++ )
      x1[i] = x1[i-1] + dx;

    for ( int i = 1; i < n2; i++ )
      x2[i] = x2[i-1] + dx;

    for ( int i = 1; i < n3; i++ )
      x3[i] = x3[i-1] + dx;

    int n12 = n1 * n2;
    for ( int i1 = 0; i1 < n1; i1++ )
    {
      for ( int i2 = 0; i2 < n2; i2++ )
      {
        for ( int i3 = 0; i3 < n3; i3++ )
        {
          grid[    (i1 + i2*n1 + i3*n12)*3] = x1[i1];
          grid[1 + (i1 + i2*n1 + i3*n12)*3] = x2[i2];
          grid[2 + (i1 + i2*n1 + i3*n12)*3] = x3[i3];
        }
      }
    }

    free(x1);
    free(x2);
    free(x3);
}


/*!
  \brief Creates a 2D grid with equal number of points (ENPs) on each axis

  \param rn    A reference to the number of elements along each axis.
  \param rx1L  A reference to the left end of the x1 range.
  \param rx1R  A reference to the right end of the x1 range.
  \param rx2L  A reference to the left end of the x2 range.
  \param rx2R  A reference to the right end of the x2 range.
  \param rf    A reference to the range extension factor.
  \param grid  The output grid.
*/
void C_create_ENPs_grid_2D(const int *rn,
                           const double *rx1L,
                           const double *rx1R,
                           const double *rx2L,
                           const double *rx2R,
                           const double *rf,
                           double *grid)
{
    int n = rn[0];
    double f = rf[0]/2;

    double x1L = rx1L[0];
    double x1R = rx1R[0];
    double x2L = rx2L[0];
    double x2R = rx2R[0];

    double Dx1 = x1R - x1L;
    double Dx2 = x2R - x2L;

    if ( f > 0 )
    {
      double dx1 = Dx1 * f;
      double dx2 = Dx2 * f;
      x1L -= dx1;
      x1R += dx1;
      x2L -= dx2;
      x2R += dx2;

      Dx1 = x1R - x1L;
      Dx2 = x2R - x2L;
    }

    double D1 = Dx1 / (n-1);
    double D2 = Dx2 / (n-1);

    double* x1 = calloc(n, sizeof(double));
    CHECK_PTR(x1);

    double* x2 = calloc(n, sizeof(double));
    CHECK_PTR(x2);

    x1[0] = x1L;
    x2[0] = x2L;

    for ( int i = 1; i < n; i++ )
    {
      x1[i] = x1[i-1] + D1;
      x2[i] = x2[i-1] + D2;
    }

    for ( int i1 = 0; i1 < n; i1++ )
    {
      for ( int i2 = 0; i2 < n; i2++ )
      {
        grid[    (i1 + i2*n)*2] = x1[i1];
        grid[1 + (i1 + i2*n)*2] = x2[i2];
      }
    }

    free(x1);
    free(x2);
}



/*!
  \brief Creates a 3D grid with equal number of points (ENPs) on each axis

  \param rn    A reference to the number of elements along each axis.
  \param rx1L  A reference to the left end of the x1 range.
  \param rx1R  A reference to the right end of the x1 range.
  \param rx2L  A reference to the left end of the x2 range.
  \param rx2R  A reference to the right end of the x2 range.
  \param rx3L  A reference to the left end of the x3 range.
  \param rx3R  A reference to the right end of the x3 range.
  \param rf    A reference to the range extension factor.
  \param grid  The output grid.
*/
void C_create_ENPs_grid_3D(const int *rn,
                           const double *rx1L,
                           const double *rx1R,
                           const double *rx2L,
                           const double *rx2R,
                           const double *rx3L,
                           const double *rx3R,
                           const double *rf,
                           double *grid)
{
    int n = rn[0];
    double f = rf[0]/2;

    double x1L = rx1L[0];
    double x1R = rx1R[0];

    double x2L = rx2L[0];
    double x2R = rx2R[0];

    double x3L = rx3L[0];
    double x3R = rx3R[0];

    double Dx1 = x1R - x1L;
    double Dx2 = x2R - x2L;
    double Dx3 = x3R - x3L;

    if ( f > 0 )
    {
      double dx1 = Dx1 * f;
      double dx2 = Dx2 * f;
      double dx3 = Dx3 * f;

      x1L -= dx1;
      x1R += dx1;

      x2L -= dx2;
      x2R += dx2;

      x3L -= dx3;
      x3R += dx3;

      Dx1 = x1R - x1L;
      Dx2 = x2R - x2L;
      Dx3 = x3R - x3L;
    }

    double D1 = Dx1 / (n-1);
    double D2 = Dx2 / (n-1);
    double D3 = Dx3 / (n-1);

    double* x1 = calloc(n, sizeof(double));
    CHECK_PTR(x1);

    double* x2 = calloc(n, sizeof(double));
    CHECK_PTR(x2);

    double* x3 = calloc(n, sizeof(double));
    CHECK_PTR(x3);

    x1[0] = x1L;
    x2[0] = x2L;
    x3[0] = x3L;

    for ( int i = 1; i < n; i++ )
    {
      x1[i] = x1[i-1] + D1;
      x2[i] = x2[i-1] + D2;
      x3[i] = x3[i-1] + D3;
    }

    int n2 = n*n;
    for ( int i1 = 0; i1 < n; i1++ )
    {
      for ( int i2 = 0; i2 < n; i2++ )
      {
        for ( int i3 = 0; i3 < n; i3++ )
        {
          grid[    (i1 + i2*n + i3*n2)*3] = x1[i1];
          grid[1 + (i1 + i2*n + i3*n2)*3] = x2[i2];
          grid[2 + (i1 + i2*n + i3*n2)*3] = x3[i3];
        }
      }
    }

    free(x1);
    free(x2);
    free(x3);
}



/*!
  \brief Creates an equidistant grid in an x-dimensional bounding box with each edge of the grid having the same length.

  \param rw        A reference to the common length of each edge.
  \param L         An array containing the left end of the range for each dimension.
  \param rdim      A reference to the number of dimensions.
  \param size_d    An array of the number of grid points along each axis.
  \param rTotElts  A reference to the number of number of elements in the grid.
  \param grid      The output grid, a pointer to a pre-allocated memory space to store the coordinates.

  The function calculates the number of grid elements along each coordinate based on the provided common edge length (rw), the left end of the range (rL), and the right end of the range (rR) for each dimension. It then creates an equidistant grid, filling the grid parameter with the coordinates.
*/
void C_create_ED_grid_xD(const double *rw,
                         const double *L,
                         const int    *rdim,
                         const int    *size_d,
                         const int    *rTotElts,
                               double *grid)
{
    double w = rw[0];
    int dim  = rdim[0];
    int totalElements = rTotElts[0];

    double** x = calloc(dim, sizeof(double*));
    CHECK_PTR(x);

    for (int i = 0; i < dim; i++) {
      x[i] = calloc(size_d[i], sizeof(double));
      CHECK_PTR(x[i]);
      x[i][0] = L[i];
      for (int j = 1; j < size_d[i]; j++)
        x[i][j] = x[i][j - 1] + w;
    }

    // Parallelize this loop using OpenMP
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < totalElements; i++) {
      int index = i;
      int indexMultiplier = 1;
      for (int d = 0; d < dim; d++) {
        grid[i * dim + d] = x[d][index / indexMultiplier % size_d[d]];
        indexMultiplier *= size_d[d];
      }
    }

    // Free coordinate arrays
    for (int i = 0; i < dim; i++)
      free(x[i]);
    free(x);
}


