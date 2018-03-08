// Author: Salvatore Filippone salvatore.filippone@cranfield.ac.uk
//

// Computes matrix-vector product. Matrix A is in row-major order
// i.e. A[i, j] is stored in i * COLS + j element of the vector.

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "wtime.h"

// Matrix dimensions.
const int ROWS = 4096;
const int COLS = 4096;


// Simple CPU implementation of matrix-vector product
void MatrixVector(int rows, int cols, const double* A, const double* x, double* restrict y) 
{
  int row,col, idx;
  double t;
  for (row = 0; row < rows; ++row) {
    t=0.0;
    for (col = 0; col < cols; ++col) {
      idx = row * cols + col;
      t = t + A[idx]*x[col];
    }
    y[row] = t;
  }
}

int main(int argc, char** argv) 
{
  
  if (argc < 3) {
    fprintf(stderr,"Usage: %s  rows cols\n",argv[0]);
  }
  int nrows=atoi(argv[1]);
  int ncols=atoi(argv[2]);
  


  
  double* A = (double*) malloc(sizeof(double)*nrows * ncols);
  double* x = (double*) malloc(sizeof(double)*nrows );
  double* y = (double*) malloc(sizeof(double)*nrows );
  int row, col, idx;
  
  srand(12345);
  for ( row = 0; row < nrows; ++row) {
    for ( col = 0; col < ncols; ++col) {
      idx = row * ncols + col;
      A[idx] = 100.0f * ((double) rand()) / RAND_MAX;
    }
    x[row] = 100.0f * ((double) rand()) / RAND_MAX;      
  }
  
  double t1 = wtime();
  MatrixVector(nrows, ncols, A, x, y);
  double t2 = wtime();
  double tmlt = (t2-t1);
  double mflops = (2.0e-6)*nrows*ncols/tmlt;
  
  fprintf(stdout,"Matrix-Vector product of size %d x %d with 1 thread: time %lf  MFLOPS %lf \n",
	  nrows,ncols,tmlt,mflops);
  free(A);
  free(x);
  free(y);
  return 0;
}
