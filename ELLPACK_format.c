#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int M, N, nz;
    int MAXNZ;
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
    /*    MAXNZ CALULATION    */
    /**************************/
    MAXNZ = 0;
    for (i = 0; i < M; i++)
    {
        int local_nz = 0;
        for (j = 0; j < N; j++)
        {
            if (mat[i][j] != 0)
            {
                local_nz++;
            }
            if (local_nz > MAXNZ)
                MAXNZ = local_nz;
        }
    }
    /**************************/

    /***********************/
    /* JA & AS CALCULATION */
    /***********************/
    int **JA = (int **)malloc(M * sizeof(int *));
    int **AS = (int **)malloc(M * sizeof(int *));
    for (int i = 0; i < M; i++)
    {
        JA[i] = (int *)malloc(MAXNZ * sizeof(int));
        AS[i] = (int *)malloc(MAXNZ * sizeof(int));
    }

    for (i = 0; i < M; i++)
    {
        int index = 0;
        int local_nz = 0;
        for (j = 0; j < N; j++)
        {
            if (mat[i][j] != 0)
            {
                local_nz++;
                if (index < MAXNZ)
                {
                    JA[i][index] = j;
                    AS[i][index] = mat[i][j];
                    index++;
                    if (local_nz < MAXNZ && JA[i][index] == 0)
                    {
                        JA[i][index] = JA[i][index - 1];
                    }
                }
            }
        }
    }

    printf("MAXNZ = %d\n", MAXNZ);

    /***************/
    /* JA PRINTING */
    /***************/
    printf("JA = \n");
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < MAXNZ; j++)
        {
            printf("%d  ", JA[i][j]);
        }
        printf("\n");
    }

    /***************/
    /* AS PRINTING */
    /***************/
    printf("AS = \n");
    for (i = 0; i < M; i++)
    {
        for (j = 0; j < MAXNZ; j++)
        {
            printf("%d  ", AS[i][j]);
        }
        printf("\n");
    }

    return 0;
}