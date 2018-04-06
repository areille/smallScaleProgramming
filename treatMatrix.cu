/*
*   This program reads a matrix with Matrix Market format,
*   then it turns it into 2 objects :
*    - one matrix in CSR format
*    - one matrix in ELLPACK format
*   Then it proceeds the mulctiplication of this matrix with a vector using cuda
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <iostream>
#include "wtime.h"
#include "mmio.h"

#include <cuda_runtime.h> // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h> // For CUDA SDK timers

//Simple dimension: define a 1D block structure
#define BD 256

const dim3 BLOCK_DIM(BD);

// GPU implementation of matrix_vector product: see if you can use
// one thread per row. You'll need to get the addressing right!
// each block of rows.
__global__ void gpuMatrixVectorCSR(int rows, int cols, const float *A,
                                const float *x, float *y)
{
    int tr = threadIdx.x;
    int row = blockIdx.x * blockDim.x + tr;
    if (row < rows)
    {
        // Starting address of indexing within matrix A
        int idxm = row * cols;
        float t = 0.0;
        for (int ic = 0; ic < cols; ic++)
        {
            t += A[idxm] * x[ic];
            idxm++;
        }
        y[row] = t;
    }
}

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
    double *val;
    bool isCsrFormat;

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename] [1 for CSR, 2 for Ellpack]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
        if (strcmp(argv[2], "1") == 0)
        {
            isCsrFormat = true;
        }
        else if (strcmp(argv[2], "2") == 0)
        {
            isCsrFormat = false;
        }
        else
        {
            printf("Second argument should be 1 for CSR or 2 for ELLPACK\n");
            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
        exit(1);

    /* reseve memory for matrices */

    I = (int *)malloc(nz * sizeof(int));
    J = (int *)malloc(nz * sizeof(int));
    val = (double *)malloc(nz * sizeof(double));

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--; /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f != stdin)
        fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i = 0; i < nz; i++)
        fprintf(stdout, "%d %d %20.19g\n", I[i], J[i], val[i]);
    printf("Columns : %d, rows : %d, non-zeros : %d\n\n", M, N, nz);

    // CREATING VECTORS
    double *x = (double *)malloc(sizeof(double) * M);
    double *y = (double *)malloc(sizeof(double) * M);

    for (i = 0; i < M; i++)
    {
        x[i] = 100.0f * ((double)rand()) / RAND_MAX;
    }
}