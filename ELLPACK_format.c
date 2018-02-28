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

    return 0;
}