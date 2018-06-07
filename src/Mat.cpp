#include "Mat.h"

void initVector(double *vec, int n, double value)
{
    for (int i = 0; i < n; i++)
        vec[i] = value;
}

void initMatrix(double **mat, int nr, int nc, double value)
{
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            mat[i][j] = value;
}

void initDiagonal(double **mat, int nr, double value)
{
    for (int i = 0; i < nr; i++)
        mat[i][i] = value;
}

void addition(double **A, double **B, double **R, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            R[i][j] = A[i][j] + B[i][j];
        }
    }
}

void additionDiagonal(double a, double **B, double **R, int n)
{
    for (int i = 0; i < n; i++)
    {
        R[i][i] = a + B[i][i];
    }
}

void addition(double *a, double *b, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = a[i] + b[i];
    }
}

void subtraction(double *a, double *b, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = a[i] - b[i];
    }
}

void subtraction(double **A, double **B, double **R, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            R[i][j] = A[i][j] - B[i][j];
        }
    }
}

void multiplication(double **P, double **Q, double **R, int nr1, int nc1, int nc2)
{
    for (int i = 0; i < nr1; i++)
    {
        for (int j = 0; j < nc2; j++)
        {
            R[i][j] = 0;
            for (int k = 0; k < nc1; k++)
            {
                R[i][j] += P[i][k] * Q[k][j];
            }
        }
    }
}

void multiplication(double **P, double **Q, double **R, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            R[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                R[i][j] += P[i][k] * Q[k][j];
            }
        }
    }
}

void multiplication(double **A, double *b, double *r, int n)
{
    multiplication(A, b, r, n, n);
}

void multiplication(double **A, double *b, double *r, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        r[i] = 0;
        for (int j = 0; j < nc; j++)
        {
            r[i] += A[i][j] * b[j];
        }
    }
}

void multiplication(double *b, double **A, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = 0;
        for (int j = 0; j < n; j++)
        {
            r[i] += b[j] * A[j][i];
        }
    }
}

void multiplication(double *a, double *b, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = a[i] * b[i];
    }
}

void multiplication(double *a, double *b, double **R, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            R[i][j] = a[i] * b[j];
        }
    }
}

double multiplication(double *a, double *b, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

void multiplication(double a, double *b, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {
        r[i] = a * b[i];
    }
}

void multiplication(double a, double **B, double **R, int n)
{
    multiplication(a, B, R, n, n);
}

void multiplication(double a, double **B, double **R, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            R[i][j] = a * B[i][j];
        }
    }
}

void divition(double a, double **B, double **R, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            R[i][j] = B[i][j] / a;
        }
    }
}

void divition(double a, double *b, double *r, int n)
{
    for (int i = 0; i < n; i++)
    {

        r[i] = b[i] / a;
    }
}

void printMatrix(double **mat, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            printf("%g ", mat[i][j]);
        }
        printf("\n");
    }
}

void printVector(double *vec, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f ", vec[i]);
    }
    printf("\n");
}

void printDiagonal(double **mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("%f ", mat[i][i]);
    }
    printf("\n");
}

double **createMatrix(int nr, int nc)
{
    double **mat = (double **)malloc((nr) * sizeof(double *));
    if (mat == NULL)
        return (NULL);
    mat[0] = (double *)malloc(nr * nc * sizeof(double));
    if (mat[0] == NULL)
        return (NULL);
    for (int i = 1; i < nr; ++i)
        mat[i] = mat[i - 1] + nc;
    return (mat);
}

double **createMatrix(int n)
{
    return createMatrix(n, n);
}

double **createMatrix(double value, int nr, int nc)
{
    double **mat = (double **)malloc((nr) * sizeof(double *));
    if (mat == NULL)
        return (NULL);
    mat[0] = (double *)malloc(nr * nc * sizeof(double));
    if (mat[0] == NULL)
        return (NULL);
    for (int i = 1; i < nr; ++i)
        mat[i] = mat[i - 1] + nc;
    setValue(mat, value, nr, nc);
    return (mat);
}

double **createMatrix(double *a, double *b, double *c, int nr, int nc)
{
    double **mat = createMatrix(nr, nc);
    for (int i = 0; i < nr; ++i)
    {
        for (int j = 0; j < nc; ++j)
        {
            mat[i][j] = 0;
        }
    }
    for (int i = 0; i < nr; ++i)
    {
        mat[i][i] = b[i];
    }
    for (int i = 0; i < nr - 1; ++i)
    {
        mat[i][i + 1] = c[i];
        mat[i + 1][i] = a[i];
    }
    return (mat);
}

void freeMatrix(double **mat)
{
    free(mat[0]);
    free(mat);
}

double **readMatrix(char *cfile, int *nr, int *nc)
{
    double **mat;
    FILE *f1 = fopen(cfile, "rb");
    if (!f1)
        return (NULL);
    fread(nr, sizeof(int), 1, f1);
    fread(nc, sizeof(int), 1, f1);
    mat = createMatrix(*nr, *nc);
    fread(mat[0], sizeof(double), (*nr) * (*nc), f1);
    fclose(f1);
    return (mat);
}

void matToFile(char *filename, double **mat, int nr, int nc)
{
    FILE *file = fopen(filename, "w");
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            fprintf(file, "%f ", mat[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void vectorToFile(const char *filename, double *vec, int n)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        printf("No se pudo crear el archivo %s\n", filename);
        return;
    }
    for (int i = 0; i < n; i++)
    {
        file << i << " " << vec[i] << std::endl;
    }
    file.close();
}

double *createVector(int n)
{
    return (double *)malloc(n * sizeof(double));
}

double *readVector(char *cfile, int *nr)
{
    double *vec;
    FILE *f1 = fopen(cfile, "rb");
    if (!f1)
        return (NULL);
    fread(nr, sizeof(int), 1, f1);
    vec = (double *)malloc((*nr) * sizeof(double));
    if (vec == NULL)
        return (NULL);
    fread(vec, sizeof(double), *nr, f1);
    fclose(f1);
    return (vec);
}

double *readVector(char *cfile, int n)
{
    double *vec;
    FILE *f1 = fopen(cfile, "rb");
    if (!f1)
        return (NULL);
    vec = (double *)malloc(n * sizeof(double));
    if (vec == NULL)
        return (NULL);
    fread(vec, sizeof(double), n, f1);
    fclose(f1);
    return (vec);
}

double *readTxtVector(const char *filename, int &n)
{
    std::ifstream file(filename);
    std::string line;
    double *vec = nullptr;
    if (file.is_open())
    {
        std::getline(file, line);
        n = std::stoi(line);
        vec = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++)
        {
            std::getline(file, line);
            vec[i] = std::stof(line);
        }
        file.close();
    }
    else
    {
        printf("Error al abrir el archivo %s\n", filename);
    }
    return (vec);
}

void readPoints(char *cfile, double *x, double *y, int *n)
{
    int nn;
    x = y = NULL;
    double *vec;
    FILE *f1 = fopen(cfile, "rb");
    if (!f1)
        return;
    fread(n, sizeof(int), 1, f1);
    fread(n, sizeof(int), 1, f1);
    nn = (*n) * 2;
    x = (double *)malloc((*n) * sizeof(double));
    y = (double *)malloc((*n) * sizeof(double));
    vec = (double *)malloc(nn * sizeof(double));
    if (vec == NULL || y == NULL || x == NULL)
        return;
    fread(vec, sizeof(double), nn, f1);
    for (int i = 0; i < (*n); i++)
    {
        x[i] = vec[i * 2];
        y[i] = vec[i * 2 + 1];
    }
    free(vec);
    fclose(f1);
}

double sumVector(double *vec, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res += vec[i];
    }
    return res;
}

double getNorm2(double *vec, int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += vec[i] * vec[i];
    }
    return sqrt(norm);
}

double getNorm1(double **mat, int nr, int nc)
{
    double norm = -1;
    for (int j = 0; j < nc; j++)
    {
        double aux = 0;
        for (int i = 0; i < nr; i++)
        {
            aux += fabs(mat[i][j]);
        }
        if (aux > norm)
        {
            norm = aux;
        }
    }
    return norm;
}

double getNormInf(double **mat, int nr, int nc)
{
    double norm = -1;
    for (int i = 0; i < nr; i++)
    {
        double aux = 0;
        for (int j = 0; j < nc; j++)
        {
            aux += fabs(mat[i][j]);
        }
        if (aux > norm)
        {
            norm = aux;
        }
    }
    return norm;
}

double getNorm2(double **mat, int nr, int nc)
{
    //No es la norma 2
    double norm = 0;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            norm += mat[i][j] * mat[i][j];
        }
    }
    return sqrt(norm);
}

double getNormF(double **mat, int nr, int nc)
{
    double norm = 0;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            norm += mat[i][j] * mat[i][j];
        }
    }
    return sqrt(norm);
}

double maximumNorm(double *vec, int n)
{
    double norm = -1;
    for (int i = 0; i < n; i++)
    {
        if (norm < fabs(vec[i]))
            norm = fabs(vec[i]);
    }
    return norm;
}

double maximumNorm(double **mat, int nr, int nc)
{
    double norm = -1;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            if (norm < fabs(mat[i][j]))
                norm = fabs(mat[i][j]);
        }
    }
    return norm;
}

bool isSymmetric(double **mat, int n)
{
    double **trans = createMatrix(n);
    transposeMatrix(mat, trans, n);
    bool rest = isSymmetric(mat, trans, n);
    freeMatrix(trans);
    return rest;
}

bool isSymmetric(double **mat, double **transposeMat, int n)
{
    double **R = createMatrix(n);
    subtraction(mat, transposeMat, R, n, n);
    double norm = getNorm1(R, n, n);
    freeMatrix(R);
    return norm < sqrt(DBL_EPSILON);
}

void normalizeVector(double *vec, int n)
{
    normalizeVector(vec, vec, n);
}

void normalizeVector(double *src, double *dst, int n)
{
    normalizeVector(src, dst, getNorm2(src, n), n);
}

void normalizeVector(double *vec, double norm, int n)
{
    normalizeVector(vec, vec, norm, n);
}

void normalizeVector(double *src, double *dst, double norm, int n)
{
    for (int i = 0; i < n; i++)
    {
        dst[i] = src[i] / norm;
    }
}

double **identityMatrix(int n)
{
    double **mat = createMatrix(n, n);
    setIdentityMatrix(mat, n);
    return mat;
}

void setIdentityMatrix(double **I, int n)
{
    setValue(I, 0., n, n);
    setDiagonal(I, 1, n);
}

void setValue(double **A, double value, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            A[i][j] = value;
        }
    }
}

void setValue(double **src, double **dst, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            dst[i][j] = src[i][j];
        }
    }
}

void copyDiagonal(double **src, double **dst, int n)
{
    for (int i = 0; i < n; i++)
    {
        dst[i][i] = src[i][i];
    }
}

void copyVector(double *src, double *dst, int n)
{
    memcpy(dst, src, n * sizeof(double));
    //for (int i = 0; i < n; i++)
    //dst[i] = src[i];
}

void setDiagonal(double **A, double value, int n)
{
    for (int i = 0; i < n; i++)
    {
        A[i][i] = value;
    }
}

void transposeMatrix(double **src, double **dst, int nr, int nc)
{
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            dst[j][i] = src[i][j];
        }
    }
}

void transposeMatrix(double **src, double **dst, int n)
{
    transposeMatrix(src, dst, n, n);
}

double **getVandermonde(int n)
{
    double **A = createMatrix(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = pow((i + 1.) / (double)n, j);
        }
    }
    return A;
}

double **getVandermonde(double *x, int n)
{
    double **A = createMatrix(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = pow(x[i], j);
        }
    }
    return A;
}

void createPositiveDefiniteMatrix(double **mat, int n)
{
    double **aux = createMatrix(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mat[i][j] = ((double)rand() / (RAND_MAX));
        }
    }
    transposeMatrix(mat, aux, n);
    addition(mat, aux, mat, n, n);
    multiplication(0.5, mat, mat, n);
    setIdentityMatrix(aux, n);
    multiplication(n, aux, aux, n);
    addition(mat, aux, mat, n, n);
    freeMatrix(aux);
}
