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
#include "mmio.h"

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

    if (isCsrFormat == true)
    {
        /*************************/
        /* CSR FORMAT CALCULATION*/
        /*************************/
        int *IRP = (int *)malloc((M + 1) * sizeof(int *));
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
        printf("---- CSR FORMAT RESULTS ----\n\n");
        /****************/
        /* IRP PRINTING */
        /****************/
        printf("IRP = [");
        for (i = 0; i < M + 1; i++)
        {
            if (i < M)
            {
                printf("%d, ", IRP[i]);
            }
            else
            {
                printf("%d", IRP[i]);
            }
        }
        printf("]\n");
        /****************/

        /****************/
        /* JA PRINTING */
        /****************/
        printf("JA = [");
        for (i = 0; i < nz; i++)
        {
            if (i < nz - 1)
            {
                printf("%d, ", J[i]+1);
            }
            else
            {
                printf("%d", J[i]+1);
            }
        }
        printf("]\n");
        /****************/

        /****************/
        /* AS PRINTING */
        /****************/
        printf("AS = [");
        for (i = 0; i < nz; i++)
        {
            if (i < nz - 1)
            {
                printf("%3.2f, ", val[i]);
            }
            else
            {
                printf("%3.2f", val[i]);
            }
        }
        printf("]\n");
        /****************/
        free(IRP);
    }
    else
    {
        /*************************/
        /* ELLPACK FORMAT CALCULATION*/
        /*************************/
        printf("hello");
    }
    return 0;
}
