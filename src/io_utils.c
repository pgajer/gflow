#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/*!
  Opens a file using the given mode. Exits the program with an appropriate message if the operation fails.
  \param file File name.
  \param mode Mode in which to open the file.
  \return FILE pointer.
*/
FILE* file_open(const char *file, const char *mode)
{
    FILE* retval = fopen(file, mode);
    if (retval == NULL) {
        Rf_error("Error in %s at line %d: Could not open %s\n", __FILE__, __LINE__, file);
    }
    return retval;
}

/*!
  Closes the file pointed to by fp. Exits the program with an appropriate message if the operation fails.
  \param fp FILE pointer.
  \param file File name.
*/
void file_close(FILE *fp, const char *file)
{
    if (fp == NULL) {
        Rf_error("Error in %s at line %d: File pointer is NULL. Cannot close file.\n", __FILE__, __LINE__);
    }

    if (fclose(fp) != 0) {
        Rf_error("Error in %s at line %d: Could not close %s\n", __FILE__, __LINE__, file ? file : "file");
    }
}

/*!
  Prints values of a double array

  \param x    A 1D double array.
  \param n    The number of elements of x.
*/
void print_double_array(const double *x, int n)
{
    for ( int i = 0; i < (n-1); i++ )
      Rprintf("%.5f, ", x[i]);
    Rprintf("%.5f\n", x[n-1]);
}

/*!
  Prints values of a double array

  \param x    A 1D double array.
  \param n    The number of elements of x.
  \param precision The print precision of doubles.
*/
void print_double_array_with_precision(const double *x, int n, int precision)
{
    // Create a format string for printing doubles with the specified precision
    char formatStr[10]; // Ensure this is large enough to hold the format string
    snprintf(formatStr, sizeof(formatStr), "%%.%df ", precision); // Creates a format string like "%.3f "

    for ( int i = 0; i < (n-1); i++ )
      Rprintf(formatStr, x[i]);
    Rprintf(formatStr, x[n-1]);
    Rprintf("\n");
}


/*!
  Writes values of a double array to a file

  \param x    A 1D double array.
  \param n    The number of elements of x.
  \param file  A file name of the file to which the array will be written in the one-column format.

*/
void write_double_array(const double *x, int n, const char *file)
{
    FILE* fh = file_open(file, "w");

    for ( int i = 0; i < n; i++ )
      fprintf(fh, "%.14lf\n",x[i]);

    file_close(fh, file);
}

/*!
  Prints values of an integer array

  \param x    A 1D int array.
  \param n    The number of elements of x.

*/
void print_int_array(const int *x, int n)
{
    for ( int i = 0; i < (n-1); i++ )
      Rprintf("%d, ", x[i]);
    Rprintf("%d\n", x[n-1]);
}

/*!
  Writes values of an integer array to a file

  \param x    A 1D int array.
  \param n    The number of elements of x.
  \param file  A file name of the file to which the array will be written in the one-column format.

*/
void write_int_array(const int *x, int n, const char *file)
{
    FILE* fh = file_open(file, "w");

    for ( int i = 0; i < n; i++ )
      fprintf(fh, "%d\n", x[i]);

    file_close(fh, file);
}

/*!
  Prints the first N rows of a double 2D array

  \param x     A double 2D array written in a column-major format.
  \param nr    The number of rows of x.
  \param nc    The number of columns of x.
  \param N     The number of rows to print.
  \param precision The print precision of doubles.
*/
void print_2d_double_array_first_n(const double *x, int nr, int nc, int N, int precision) {
    // Create a format string for printing doubles with the specified precision
    char formatStr[10]; // Ensure this is large enough to hold the format string
    snprintf(formatStr, sizeof(formatStr), "%%.%df ", precision); // Creates a format string like "%.3f "

    for (int i = 0; i < nr && i < N; i++) {
        for (int j = 0; j < nc; j++) {
            Rprintf(formatStr, *(x + i + j * nr));
        }
        Rprintf("\n");
    }
}

/*!
  Prints values of a double 2D array

  \param x     A double 2D array written in a column-major format.
  \param nr    The number of rows of x.
  \param nc    The number of columns of x.

*/
void print_2d_double_array(const double *x, int nr, int nc)
{
    //int nc1 = nc - 1;
    for ( int i = 0; i < nr; i++ )
    {
      for ( int j = 0; j < nc; j++ )
      {
        Rprintf("%.3f ",*(x + i + j*nr));
      }
      //Rprintf("%.3f\n",*(x + i + nc1*nr));
      Rprintf("\n");
    }
}

/*!
  Writes values of a double 2D array to a CSV file

  \param x     A double 2D array written in a column-major format.
  \param nr    The number of rows of x.
  \param nc    The number of columns of x.
  \param file  A file name.

*/
void write_2d_double_array(const double *x, int nr, int nc, const char *file)
{
    FILE* fh = file_open(file, "w");

    int nc1 = nc - 1;
    for ( int i = 0; i < nr; i++ )
    {
      for ( int j = 0; j < nc1; j++ )
      {
        fprintf(fh, "%.10lf,",*(x + i + j*nr));
      }
      fprintf(fh, "%.10lf\n",*(x + i + nc1*nr));
    }

    file_close(fh, file);
}

/*!
  Writes values of an int 2D array to a CSV file

  \param x     An int 2D array written in a column-major format.
  \param nr    The number of rows of x.
  \param nc    The number of columns of x.
  \param file  A file name.

*/
void write_2d_int_array(const int *x, int nr, int nc, const char *file)
{
    FILE* fh = file_open(file, "w");

    int nc1 = nc - 1;
    for ( int i = 0; i < nr; i++ )
    {
      for ( int j = 0; j < nc1; j++ )
      {
        fprintf(fh, "%d,",*(x + i + j*nr));
      }
      fprintf(fh, "%d",*(x + i + nc1*nr));
      fprintf(fh, "\n");
    }

    file_close(fh, file);
}

/*!
    \brief Writes a 2D array of doubles to a CSV file.

    Writes a 2D array of doubles along with its associated row lengths to a CSV file.
    The first row of the file stores the number of rows (nX).
    The second row stores the lengths of each row (row_len).
    The subsequent rows store the elements of each row of X.

    \param file The name of the file to which the 2D array will be written.
    \param X Pointer to the 2D array of doubles.
    \param row_len Pointer to an array containing the lengths of each row in X.
    \param nX The number of rows in the 2D array X.
*/
void write_double2_csv(const char* file, const double** X, int* row_len, int nX) {
    FILE* fp = fopen(file, "w");
    if (fp == NULL) {
        return;
    }

    // Write the value of nX to the first row
    fprintf(fp, "%d\n", nX);

    // Write the row lengths to the second row
    for (int i = 0; i < nX; i++) {
        if (i != 0) {
            fprintf(fp, ",");
        }
        fprintf(fp, "%d", row_len[i]);
    }
    fprintf(fp, "\n");

    // Write the remaining rows for elements of X
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < row_len[i]; j++) {
            if (j != 0) {
                fprintf(fp, ",");
            }
            fprintf(fp, "%f", X[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/*!
    \brief Reads a 2D array of doubles from a CSV file.

    Reads a 2D array of doubles and its associated row lengths from a CSV file.
    The first row of the file must contain the number of rows (nX).
    The second row must contain the lengths of each row (row_len).
    The subsequent rows contain the elements of each row of X.

    \param file The name of the file from which the 2D array will be read.
    \param X Pointer to a pointer that will hold the 2D array of doubles read from the file.
    \param row_len Pointer to a pointer that will hold the array containing the lengths of each row in X.
    \param nX Pointer to an integer that will hold the number of rows in the 2D array X.
*/
void read_double2_csv(const char* file, double*** X, int** row_len, int* nX) {
    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
        return;
    }

    // Read the value of nX from the first row
    int numScanned = fscanf(fp, "%d\n", nX);
    if (numScanned != 1) {
        Rf_warning("Failed to read the number of rows in the file.\n");
    }

    // Allocate memory for row_len
    *row_len = (int*) malloc((*nX) * sizeof(int));

    // Read the second row to populate row_len
    for (int i = 0; i < *nX; i++) {
        int rsize = fscanf(fp, "%d,", &((*row_len)[i]));
        if (rsize != 1) {
            Rf_warning("Failed to read the row_len.\n");
        }
    }

    // Allocate memory for X based on nX and row_len
    *X = (double**) malloc((*nX) * sizeof(double*));
    for (int i = 0; i < *nX; i++) {
        (*X)[i] = (double*) malloc((*row_len)[i] * sizeof(double));
    }

    // Read the remaining rows to populate X
    int x;
    for (int i = 0; i < *nX; i++) {
        for (int j = 0; j < (*row_len)[i]; j++) {
            x = fscanf(fp, "%lf,", &((*X)[i][j]));
            if (x != 1) {
                Rf_warning("Failed to read the %d-th element of the %d-th row.\n", j, i);
            }
        }
    }

    fclose(fp);
}


/*!
    \brief Writes a 2D array of ints to a CSV file.

    Writes a 2D array of ints along with its associated row lengths to a CSV file.
    The first row of the file stores the number of rows (nX).
    The second row stores the lengths of each row (row_len).
    The subsequent rows store the elements of each row of X.

    \param file The name of the file to which the 2D array will be written.
    \param X Pointer to the 2D array of ints.
    \param row_len Pointer to an array containing the lengths of each row in X.
    \param nX The number of rows in the 2D array X.
*/
void write_int2_csv(const char* file, const int** X, int* row_len, int nX) {
    FILE* fp = fopen(file, "w");
    if (fp == NULL) {
        return;
    }

    // Write the value of nX to the first row
    fprintf(fp, "%d\n", nX);

    // Write the row lengths to the second row
    for (int i = 0; i < nX; i++) {
        if (i != 0) {
            fprintf(fp, ",");
        }
        fprintf(fp, "%d", row_len[i]);
    }
    fprintf(fp, "\n");

    // Write the remaining rows for elements of X
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < row_len[i]; j++) {
            if (j != 0) {
                fprintf(fp, ",");
            }
            fprintf(fp, "%d", X[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/*!
    \brief Reads a 2D array of ints from a CSV file.

    Reads a 2D array of ints and its associated row lengths from a CSV file.
    The first row of the file must contain the number of rows (nX).
    The second row must contain the lengths of each row (row_len).
    The subsequent rows contain the elements of each row of X.

    \param file The name of the file from which the 2D array will be read.
    \param X Pointer to a pointer that will hold the 2D array of ints read from the file.
    \param row_len Pointer to a pointer that will hold the array containing the lengths of each row in X.
    \param nX Pointer to an integer that will hold the number of rows in the 2D array X.
*/
void read_int2_csv(const char* file, int*** X, int** row_len, int* nX) {
    FILE* fp = fopen(file, "r");
    if (fp == NULL) {
        return;
    }

    // Read the value of nX from the first row
    int numScanned = fscanf(fp, "%d\n", nX);
    if (numScanned != 1) {
        Rf_warning("Failed to read the number of rows in the file.\n");
    }


    // Allocate memory for row_len
    *row_len = (int*) malloc((*nX) * sizeof(int));

    // Read the second row to populate row_len
    for (int i = 0; i < *nX; i++) {
        int rsize = fscanf(fp, "%d,", &((*row_len)[i]));
        if (rsize != 1) {
            Rf_warning("Failed to read the row_len.\n");
        }
    }

    // Allocate memory for X based on nX and row_len
    *X = (int**) malloc((*nX) * sizeof(int*));
    for (int i = 0; i < *nX; i++) {
        (*X)[i] = (int*) malloc((*row_len)[i] * sizeof(int));
    }

    // Read the remaining rows to populate X
    int x;
    for (int i = 0; i < *nX; i++) {
        for (int j = 0; j < (*row_len)[i]; j++) {
            x = fscanf(fp, "%d,", &((*X)[i][j]));
            if (x != 1) {
                Rf_warning("Failed to read the %d-th element of the %d-th row.\n", j, i);
            }
        }
    }

    fclose(fp);
}
