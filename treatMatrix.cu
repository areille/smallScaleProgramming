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
#include <string.h>
#include <ctype.h>

#include <cuda_runtime.h> // For CUDA runtime API
#include <helper_cuda.h>  // For checkCudaError macro
#include <helper_timer.h> // For CUDA SDK timers

//Simple dimension: define a 1D block structure
#define BD 256

const dim3 BLOCK_DIM(BD);

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];
/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode) ((typecode)[0] == 'M')

#define mm_is_sparse(typecode) ((typecode)[1] == 'C')
#define mm_is_coordinate(typecode) ((typecode)[1] == 'C')
#define mm_is_dense(typecode) ((typecode)[1] == 'A')
#define mm_is_array(typecode) ((typecode)[1] == 'A')

#define mm_is_complex(typecode) ((typecode)[2] == 'C')
#define mm_is_real(typecode) ((typecode)[2] == 'R')
#define mm_is_pattern(typecode) ((typecode)[2] == 'P')
#define mm_is_integer(typecode) ((typecode)[2] == 'I')

#define mm_is_symmetric(typecode) ((typecode)[3] == 'S')
#define mm_is_general(typecode) ((typecode)[3] == 'G')
#define mm_is_skew(typecode) ((typecode)[3] == 'K')
#define mm_is_hermitian(typecode) ((typecode)[3] == 'H')

int mm_is_valid(MM_typecode matcode); /* too complex for a macro */

/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode) ((*typecode)[0] = 'M')
#define mm_set_coordinate(typecode) ((*typecode)[1] = 'C')
#define mm_set_array(typecode) ((*typecode)[1] = 'A')
#define mm_set_dense(typecode) mm_set_array(typecode)
#define mm_set_sparse(typecode) mm_set_coordinate(typecode)

#define mm_set_complex(typecode) ((*typecode)[2] = 'C')
#define mm_set_real(typecode) ((*typecode)[2] = 'R')
#define mm_set_pattern(typecode) ((*typecode)[2] = 'P')
#define mm_set_integer(typecode) ((*typecode)[2] = 'I')

#define mm_set_symmetric(typecode) ((*typecode)[3] = 'S')
#define mm_set_general(typecode) ((*typecode)[3] = 'G')
#define mm_set_skew(typecode) ((*typecode)[3] = 'K')
#define mm_set_hermitian(typecode) ((*typecode)[3] = 'H')

#define mm_clear_typecode(typecode) ((*typecode)[0] = (*typecode)[1] = \
                                         (*typecode)[2] = ' ',         \
                                     (*typecode)[3] = 'G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)

/********************* Matrix Market error codes ***************************/

#define MM_COULD_NOT_READ_FILE 11
#define MM_PREMATURE_EOF 12
#define MM_NOT_MTX 13
#define MM_NO_HEADER 14
#define MM_UNSUPPORTED_TYPE 15
#define MM_LINE_TOO_LONG 16
#define MM_COULD_NOT_WRITE_FILE 17

/******************** Matrix Market internal definitions ********************

   MM_matrix_typecode: 4-character sequence

				    ojbect 		sparse/   	data        storage 
						  		dense     	type        scheme

   string position:	 [0]        [1]			[2]         [3]

   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)

 ***********************************************************************/

#define MM_MTX_STR "matrix"
#define MM_ARRAY_STR "array"
#define MM_DENSE_STR "array"
#define MM_COORDINATE_STR "coordinate"
#define MM_SPARSE_STR "coordinate"
#define MM_COMPLEX_STR "complex"
#define MM_REAL_STR "real"
#define MM_INT_STR "integer"
#define MM_GENERAL_STR "general"
#define MM_SYMM_STR "symmetric"
#define MM_HERM_STR "hermitian"
#define MM_SKEW_STR "skew-symmetric"
#define MM_PATTERN_STR "pattern"

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            return MM_PREMATURE_EOF;
    } while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;

    else
        do
        {
            num_items_read = fscanf(f, "%d %d %d", M, N, nz);
            if (num_items_read == EOF)
                return MM_PREMATURE_EOF;
        } while (num_items_read != 3);

    return 0;
}

char *mm_typecode_to_str(MM_typecode matcode)
{
    char buffer[MM_MAX_LINE_LENGTH];
    char *types[4];
    char *mm_strdup(const char *);
    int error = 0;

    /* check for MTX type */
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;
    else
        error = 1;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;

    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer, "%s %s %s %s", types[0], types[1], types[2], types[3]);
    return mm_strdup(buffer);
}

int mm_read_banner(FILE *f, MM_typecode *matcode)
{
    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    mm_clear_typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
               storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p = mtx; *p != '\0'; *p = tolower(*p), p++)
        ; /* convert to lower case */
    for (p = crd; *p != '\0'; *p = tolower(*p), p++)
        ;
    for (p = data_type; *p != '\0'; *p = tolower(*p), p++)
        ;
    for (p = storage_scheme; *p != '\0'; *p = tolower(*p), p++)
        ;

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);

    /* second field describes whether this is a sparse matrix (in coordinate
storgae) or a dense array */

    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else if (strcmp(crd, MM_DENSE_STR) == 0)
        mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    /* third field */

    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    /* fourth field */

    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}

int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
                               double **val_, int **I_, int **J_)
{
    FILE *f;
    MM_typecode matcode;
    int M, N, nz;
    int i;
    double *val;
    int *I, *J;

    if ((f = fopen(fname, "r")) == NULL)
        return -1;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", fname);
        return -1;
    }

    if (!(mm_is_real(matcode) && mm_is_matrix(matcode) &&
          mm_is_sparse(matcode)))
    {
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n",
                mm_typecode_to_str(matcode));
        return -1;
    }

    /* find out size of sparse matrix: M, N, nz .... */

    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0)
    {
        fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
        return -1;
    }

    *M_ = M;
    *N_ = N;
    *nz_ = nz;

    /* reseve memory for matrices */

    I = (int *)malloc(nz * sizeof(int));
    J = (int *)malloc(nz * sizeof(int));
    val = (double *)malloc(nz * sizeof(double));

    *val_ = val;
    *I_ = I;
    *J_ = J;

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--; /* adjust from 1-based to 0-based */
        J[i]--;
    }
    fclose(f);

    return 0;
}

int mm_is_valid(MM_typecode matcode)
{
    if (!mm_is_matrix(matcode))
        return 0;
    if (mm_is_dense(matcode) && mm_is_pattern(matcode))
        return 0;
    if (mm_is_real(matcode) && mm_is_hermitian(matcode))
        return 0;
    if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
                                   mm_is_skew(matcode)))
        return 0;
    return 1;
}

int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
{
    if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;
    /* set return null parameter values, in case we exit with errors */
    *M = *N = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            return MM_PREMATURE_EOF;
    } while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d", M, N) == 2)
        return 0;

    else /* we have a blank line */
        do
        {
            num_items_read = fscanf(f, "%d %d", M, N);
            if (num_items_read == EOF)
                return MM_PREMATURE_EOF;
        } while (num_items_read != 2);

    return 0;
}

int mm_write_mtx_array_size(FILE *f, int M, int N)
{
    if (fprintf(f, "%d %d\n", M, N) != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

/*-------------------------------------------------------------------------*/

/******************************************************************/
/* use when I[], J[], and val[]J, and val[] are already allocated */
/******************************************************************/

int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int I[], int J[],
                         double val[], MM_typecode matcode)
{
    int i;
    if (mm_is_complex(matcode))
    {
        for (i = 0; i < nz; i++)
            if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2 * i], &val[2 * i + 1]) != 4)
                return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
        for (i = 0; i < nz; i++)
        {
            if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]) != 3)
                return MM_PREMATURE_EOF;
        }
    }

    else if (mm_is_pattern(matcode))
    {
        for (i = 0; i < nz; i++)
            if (fscanf(f, "%d %d", &I[i], &J[i]) != 2)
                return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}

int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
                          double *real, double *imag, MM_typecode matcode)
{
    if (mm_is_complex(matcode))
    {
        if (fscanf(f, "%d %d %lg %lg", I, J, real, imag) != 4)
            return MM_PREMATURE_EOF;
    }
    else if (mm_is_real(matcode))
    {
        if (fscanf(f, "%d %d %lg\n", I, J, real) != 3)
            return MM_PREMATURE_EOF;
    }

    else if (mm_is_pattern(matcode))
    {
        if (fscanf(f, "%d %d", I, J) != 2)
            return MM_PREMATURE_EOF;
    }
    else
        return MM_UNSUPPORTED_TYPE;

    return 0;
}

/************************************************************************
mm_read_mtx_crd()  fills M, N, nz, array of values, and return
type code, e.g. 'MCRS'

if matrix is complex, values[] is of size 2*nz,
 (nz pairs of real/imaginary values)
************************************************************************/

int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **I, int **J,
                    double **val, MM_typecode *matcode)
{
    int ret_code;
    FILE *f;

    if (strcmp(fname, "stdin") == 0)
        f = stdin;
    else if ((f = fopen(fname, "r")) == NULL)
        return MM_COULD_NOT_READ_FILE;

    if ((ret_code = mm_read_banner(f, matcode)) != 0)
        return ret_code;

    if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
          mm_is_matrix(*matcode)))
        return MM_UNSUPPORTED_TYPE;

    if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
        return ret_code;

    *I = (int *)malloc(*nz * sizeof(int));
    *J = (int *)malloc(*nz * sizeof(int));
    *val = NULL;

    if (mm_is_complex(*matcode))
    {
        *val = (double *)malloc(*nz * 2 * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);
        if (ret_code != 0)
            return ret_code;
    }
    else if (mm_is_real(*matcode))
    {
        *val = (double *)malloc(*nz * sizeof(double));
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);
        if (ret_code != 0)
            return ret_code;
    }

    else if (mm_is_pattern(*matcode))
    {
        ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *I, *J, *val,
                                        *matcode);
        if (ret_code != 0)
            return ret_code;
    }

    if (f != stdin)
        fclose(f);
    return 0;
}

int mm_write_banner(FILE *f, MM_typecode matcode)
{
    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code != 2)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int I[], int J[],
                     double val[], MM_typecode matcode)
{
    FILE *f;
    int i;

    if (strcmp(fname, "stdout") == 0)
        f = stdout;
    else if ((f = fopen(fname, "w")) == NULL)
        return MM_COULD_NOT_WRITE_FILE;

    /* print banner followed by typecode */
    fprintf(f, "%s ", MatrixMarketBanner);
    fprintf(f, "%s\n", mm_typecode_to_str(matcode));

    /* print matrix sizes and nonzeros */
    fprintf(f, "%d %d %d\n", M, N, nz);

    /* print values */
    if (mm_is_pattern(matcode))
        for (i = 0; i < nz; i++)
            fprintf(f, "%d %d\n", I[i], J[i]);
    else if (mm_is_real(matcode))
        for (i = 0; i < nz; i++)
            fprintf(f, "%d %d %20.16g\n", I[i], J[i], val[i]);
    else if (mm_is_complex(matcode))
        for (i = 0; i < nz; i++)
            fprintf(f, "%d %d %20.16g %20.16g\n", I[i], J[i], val[2 * i],
                    val[2 * i + 1]);
    else
    {
        if (f != stdout)
            fclose(f);
        return MM_UNSUPPORTED_TYPE;
    }

    if (f != stdout)
        fclose(f);

    return 0;
}

/**
*  Create a new copy of a string s.  mm_strdup() is a common routine, but
*  not part of ANSI C, so it is included here.  Used by mm_typecode_to_str().
*
*/
char *mm_strdup(const char *s)
{
    int len = strlen(s);
    char *s2 = (char *)malloc((len + 1) * sizeof(char));
    return strcpy(s2, s);
}

// GPU implementation of matrix_vector product: see if you can use
// one thread per row. You'll need to get the addressing right!
// each block of rows.
__global__ void gpuMatrixVector(int rows, int cols, const float *A,
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

__global__ void MatrixVectorCSRParallel(int M, const int *IRP, const int *JA, const double *AS, const double *x, double *y)
{
    int tr = threadIdx.x;
    int m = blockIdx.x * blockDim.x + tr;
    if (m < M)
    {
        double sum = 0.0;
        for (j = IRP[m]; j < IRP[m + 1]; ++j)
        {
            sum += AS[j] * x[JA[j]]
        }
        y[m] = sum;
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
                if (index <= M)
                {
                    IRP[index + 1] = IRP[index] + local_row_nz;
                    local_row_nz = 1;
                    index++;
                }
            }
        }

        // ----------------------- Host memory initialisation ----------------------- //

        double *h_y_d = new double[M];
        std::cout << "Matrix-vector product: single thread per row version " << std::endl;
        std::cout << "Test case: " << M << " x " << M << std::endl;

        // ---------------------- Device memory initialisation ---------------------- //
        //  Allocate memory space on the device.

        int *d_IRP, *d_J;
        double *d_val, *d_x, *d_y;

        checkCudaErrors(cudaMalloc((void **)&d_IRP, (M+1) * sizeof(int)));
        checkCudaErrors(cudaMalloc((void **)&d_J, nz * sizeof(int)));
        checkCudaErrors(cudaMalloc((void **)&d_val, nz * sizeof(double)));
        checkCudaErrors(cudaMalloc((void **)&d_x, M * sizeof(double)));
        checkCudaErrors(cudaMalloc((void **)&d_y, M * sizeof(double)));

        // Copy matrices from the host (CPU) to the device (GPU).
        
        checkCudaErrors(cudaMemcpy(d_IRP, IRP, (M+1) * sizeof(int), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(d_J, J, nz * sizeof(int), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(d_val, val, nz * sizeof(double), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(d_x, x, M * sizeof(double), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(d_y, y, M * sizeof(double), cudaMemcpyHostToDevice));


        // ------------------------ Calculations on the CPU ------------------------- //
        float flopcnt = 2.e-6 * nrows * ncols;

        // Create the CUDA SDK timer.
        StopWatchInterface *timer = 0;
        sdkCreateTimer(&timer);

        timer->start();
        // CpuMatrixVector(nrows, ncols, h_A, h_x, h_y);
        MatrixVectorCSR(M, IRP, J, val, x, y);

        timer->stop();
        float cpuflops = flopcnt / timer->getTime();
        std::cout << "  CPU time: " << timer->getTime() << " ms."
                  << " GFLOPS " << cpuflops << std::endl;

        // ------------------------ Calculations on the GPU ------------------------- //

        // Calculate the dimension of the grid of blocks (1D) necessary to cover
        // all rows.
        const dim3 GRID_DIM((nrows - 1 + BLOCK_DIM.x) / BLOCK_DIM.x, 1);

        timer->reset();
        timer->start();
        MatrixVectorCSRParallel<<<GRID_DIM, BLOCK_DIM>>>(M, d_IRP, d_J, d_val, d_x, d_y);
        checkCudaErrors(cudaDeviceSynchronize());

        timer->stop();
        float gpuflops = flopcnt / timer->getTime();
        std::cout << "  GPU time: " << timer->getTime() << " ms."
                  << " GFLOPS " << gpuflops << std::endl;

        // Download the resulting vector d_y from the device and store it in h_y_d.
        checkCudaErrors(cudaMemcpy(h_y_d, d_y, nrows * sizeof(float), cudaMemcpyDeviceToHost));

        // Now let's check if the results are the same.
        float reldiff = 0.0f;
        float diff = 0.0f;

        for (int i = 0; i < M; ++i)
        {
            float maxabs = std::max(std::abs(y[i]), std::abs(h_y_d[i]));
            if (maxabs == 0.0)
                maxabs = 1.0;
            reldiff = std::max(reldiff, std::abs(y[i] - h_y_d[i]) / maxabs);
            diff = std::max(diff, std::abs(y[i] - h_y_d[i]));
        }
        std::cout << "Max diff = " << diff << "  Max rel diff = " << reldiff << std::endl;
        // Rel diff should be as close as possible to unit roundoff; float
        // corresponds to IEEE single precision, so unit roundoff is
        // 1.19e-07
        //

        // ------------------------------- Cleaning up ------------------------------ //

        delete timer;

        checkCudaErrors(cudaFree(d_A));
        checkCudaErrors(cudaFree(d_x));
        checkCudaErrors(cudaFree(d_y));

        delete[] h_A;
        delete[] h_x;
        delete[] h_y;
        delete[] h_y_d;
    }
    return 0;
}