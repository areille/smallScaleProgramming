/*
*   This program reads a matrix with Matrix Market format,
*   then it turns it into 2 objects :
*    - one matrix in CSR format
*    - one matrix in ELLPACK format
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include "wtime.h"
#include "mmio.h"

// Simple CPU implementation of matrix-vector product
void MatrixVectorCSR(int M, const int *IRP, const int *JA, const double *AS, const double *x, double *y)
{
    int i, j;
    double t;
    for (i = 0; i < M; ++i)
    {
        t = 0.0;
        for (j = IRP[i]; j < IRP[i + 1]; ++j)
        {
            t += AS[j] * x[JA[j]];
        }
        y[i] = t;
    }
}

inline double dmin(double a, double b) { return a < b ? a : b; }

inline int max(int a, int b) { return a > b ? a : b; }
inline int min(int a, int b) { return a < b ? a : b; }

// muylti cpu
const int ntimes = 50;
const int BBS = 1000;

void MatrixVectorCSRParallel(int M, const int *IRP, const int *JA, const double *AS, const double *x, double *y)
{
    int i, j;
    double sum;
#pragma omp parallel default(none) shared(x, y, M, IRP, JA, AS) private(i, j, sum)
    for (i = 0; i < M; ++i)
    {
        sum = 0.0;
#pragma omp for
        for (j = IRP[i]; j < IRP[i + 1]; ++j)
        {
            sum += AS[j] * x[JA[j]];
        }
#pragma omp critical
        y[i] = sum;
    }
}

void MatrixVectorELL(int M, const int MAXNZ, int **JA, double **AS, const double *x, double *y)
{
    int i, j;
    double t;
    for (i = 0; i < M; ++i)
    {
        t = 0.0;
        for (j = 0; j < MAXNZ; ++j)
        {
            t += AS[i][j] * x[JA[i][j]];
        }
        y[i] = t;
    }
}

void MatrixVectorELLParallel(int M, const int MAXNZ, int **JA, double **AS, const double *x, double *y)
{
    int i, j;
    double sum;
#pragma omp parallel default(none) shared(x, y, M, MAXNZ, JA, AS) private(i, j, sum)
    for (i = 0; i < M; ++i)
    {
        sum = 0.0;
#pragma omp for
        for (j = IRP[i]; j < IRP[i + 1]; ++j)
        {
            sum += AS[i][j] * x[JA[i][j]];
        }
#pragma omp critical
        y[i] = sum;
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
    if (isCsrFormat == true)
    {
        /*************************/
        /* CSR FORMAT CALCULATION*/
        /*************************/
        int *IRP = (int *)malloc((M + 1) * sizeof(int));
        // ASSUMING MATLAB FIRST COLUMN INDEXING
        IRP[0] = 1;
        int index = 0;
        int local_row_nz = 1;
        for (i = 0; i < nz; i++)
        {
            if (I[i] == I[i + 1])
            {
                local_row_nz++;
            }
            else
            {
                IRP[index + 1] = IRP[index] + local_row_nz;
                local_row_nz = 1;
                index++;
            }
        }
        // printf("---- CSR FORMAT RESULTS ----\n\n");
        // /****************/
        // /* IRP PRINTING */
        // /****************/
        // printf("IRP = [");
        // for (i = 0; i < M + 1; i++)
        // {
        //     if (i < M)
        //     {
        //         printf("%d, ", IRP[i]);
        //     }
        //     else
        //     {
        //         printf("%d", IRP[i]);
        //     }
        // }
        // printf("]\n");
        // /****************/

        // /****************/
        // /* JA PRINTING */
        // /****************/
        // printf("JA = [");
        // for (i = 0; i < nz; i++)
        // {
        //     if (i < nz - 1)
        //     {
        //         printf("%d, ", J[i] + 1);
        //     }
        //     else
        //     {
        //         printf("%d", J[i] + 1);
        //     }
        // }
        // printf("]\n");
        // /****************/

        // /****************/
        // /* AS PRINTING */
        // /****************/
        // printf("AS = [");
        // for (i = 0; i < nz; i++)
        // {
        //     if (i < nz - 1)
        //     {
        //         printf("%3.2f, ", val[i]);
        //     }
        //     else
        //     {
        //         printf("%3.2f", val[i]);
        //     }
        // }
        // printf("]\n");
        // /****************/

        double tmlt = 1e100;
        for (int try = 0; try < ntimes; try ++)
        {
            double t1 = wtime();
            MatrixVectorCSRParallel(M, IRP, J, val, x, y);
            double t2 = wtime();
            tmlt = dmin(tmlt, (t2 - t1));
        }
        double mflops = (2.0e-6) * nrows / tmlt;
#pragma omp parallel
        {
#pragma omp master
            {
                fprintf(stdout, "Matrix-Vector product (block_unroll_2) of size %d with %d threads: time %lf  MFLOPS %lf \n",
                        M, omp_get_num_threads(), tmlt, mflops);
            }
        }
        free(IRP);
    }
    else
    {
        /*************************/
        /* ELLPACK FORMAT CALCULATION*/
        /*************************/
        int MAXNZ;
        /**************************/
        /*    MAXNZ CALULATION    */
        /**************************/
        MAXNZ = 0;
        int local_row_nz = 1;
        for (i = 0; i < nz; i++)
        {
            if (I[i] == I[i + 1])
            {
                local_row_nz++;
            }
            else
            {
                if (local_row_nz > MAXNZ)
                {
                    MAXNZ = local_row_nz;
                }
                local_row_nz = 1;
            }
        }
        printf("MAXNZ : %d\n", MAXNZ);
        /**************************/

        /***********************/
        /* JA & AS CALCULATION */
        /***********************/
        int **JA = (int **)malloc((M + 1) * sizeof(int *));
        double **AS = (double **)malloc(M * sizeof(double *));
        for (i = 0; i < M; i++)
        {
            JA[i] = (int *)malloc(MAXNZ * sizeof(int));
            AS[i] = (double *)malloc(MAXNZ * sizeof(double));
        }
        int *locals = (int *)malloc(M * sizeof(int));

        i = 0;     // from 0 to nz
        int j = 0; // from 0 to M
        int k = 0; // from 0 to MAXNZ

        local_row_nz = 1;
        int index = 0;
        for (i = 0; i < nz; i++)
        {
            if (I[i] == index && I[i + 1] == index)
            {
                local_row_nz++;
            }
            else
            {
                locals[index] = local_row_nz;
                local_row_nz = 1;
                index++;
            }
        }
        index = 0;
        for (i = 0; i < nz; i++)
        {
            if (locals[index] == MAXNZ)
            {
                if (k < MAXNZ)
                {
                    JA[j][k] = J[i] + 1;
                    AS[j][k] = val[i];
                    k++;
                }
                else
                {
                    k = 0;
                    if (j < M)
                    {
                        j++;
                        JA[j][k] = J[i] + 1;
                        AS[j][k] = val[i];
                        k++;
                        index++;
                    }
                }
            }
            else
            {
                if (k < locals[index])
                {
                    JA[j][k] = J[i] + 1;
                    AS[j][k] = val[i];
                    printf("k : %d", k);
                    k++;
                    // printf("j : %d, k : %d, local : %d\n", j, k, locals[index]);
                }
                else
                {
                    k = 0;
                    if (j < M - 1)
                    {
                        j++;
                        JA[j][k] = J[i] + 1;
                        AS[j][k] = val[i];
                        k++;
                        index++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        // /***************/
        // /* JA PRINTING */
        // /***************/
        // printf("JA = \n");
        // for (i = 0; i < M; i++)
        // {
        //     for (j = 0; j < MAXNZ; j++)
        //     {
        //         printf("%d  ", JA[i][j]);
        //     }
        //     printf("\n");
        // }

        // /***************/
        // /* AS PRINTING */
        // /***************/
        // printf("AS = \n");
        // for (i = 0; i < M; i++)
        // {
        //     for (j = 0; j < MAXNZ; j++)
        //     {
        //         printf("%3.2f  ", AS[i][j]);
        //     }
        //     printf("\n");
        // }

        double tmlt = 1e100;
        for (int try = 0; try < ntimes; try ++)
        {
            double t1 = wtime();
            // MatrixVectorELL(M, MAXNZ, JA, AS, x, y);
            MatrixVectorELLParallel(M, MAXNZ, JA, AS, x, y);
            double t2 = wtime();
            tmlt = dmin(tmlt, (t2 - t1));
        }
        double mflops = (2.0e-6) * nrows / tmlt;
#pragma omp parallel
        {
#pragma omp master
            {
                fprintf(stdout, "Matrix-Vector product (block_unroll_2) of size %d with %d threads: time %lf  MFLOPS %lf \n",
                        M, omp_get_num_threads(), tmlt, mflops);
            }
        }
        free(JA);
        free(AS);
        free(locals);
    }
    free(x);
    free(y);
    free(I);
    free(J);
    free(val);
    return 0;
}
