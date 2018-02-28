#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{

    int M, N, nz;
    int i, j;
    nz = 0;

    printf("Number of rows :\n");
    scanf("%d", &M);
    printf("Number of cols :\n");
    scanf("%d", &N);

    int **mat = (int **)malloc(M * sizeof(int *));
    for (int i = 0; i < M; i++)
        mat[i] = (int *)malloc(N * sizeof(int));

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%dth row, %dth col : \n", i, j);
            scanf("%d", &mat[i][j]);
            if (mat[i][j] != 0)
                nz++;
        }
    }

    printf("Matrix : \n");
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            printf("%d ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /**************************/
    /* IRP, JA & AS CALULATION*/
    /**************************/
    int *IRP = (int *)malloc((M + 1) * sizeof(int *));
    IRP[0] = 0;
    int *JA = (int *)malloc(nz * sizeof(int *));
    int *AS = (int *)malloc(nz * sizeof(int *));
    int index = 0;

    for (i = 0; i < M; i++)
    {
        int local_nz = 0;
        for (j = 0; j < N; j++)
        {
            if (mat[i][j] != 0)
            {
                local_nz++;
                if (index < nz)
                {
                    JA[index] = j;
                    AS[index] = mat[i][j];
                    index++;
                }
            }
        }
        IRP[i + 1] = IRP[i] + local_nz;
    }
    /*****************/

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
            printf("%d, ", JA[i]);
        }
        else
        {
            printf("%d", JA[i]);
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
            printf("%d, ", AS[i]);
        }
        else
        {
            printf("%d", AS[i]);
        }
    }
    printf("]\n");
    /****************/

    free(mat);
    free(IRP);
    free(JA);
    free(AS);
    return 0;
}