#include "linearAlgebra.hpp"

bool solveDiagonal(double *d, double *x, double *b, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (fabs(d[i]) < tolerance)
        {
            return false;
        }
        x[i] = b[i] / d[i];
    }
    return true;
}

bool solveTriDiagonal(double *a, double *b, double *c, double *d, double *x, int n)
{
    double *ab, *bb, *cb, *db;
    ab = (double *)malloc(n * sizeof(double));
    bb = (double *)malloc(n * sizeof(double));
    cb = (double *)malloc(n * sizeof(double));
    db = (double *)malloc(n * sizeof(double));
    bb[0] = b[0];
    cb[0] = c[0];
    db[0] = d[0];
    for (int i = 1; i < n; i++)
    {
        bb[i] = bb[i - 1] * b[i] - a[i] * cb[i - 1];
        cb[i] = bb[i - 1] * c[i];
        db[i] = bb[i - 1] * d[i] - a[i] * db[i - 1];
    }
    if (fabs(bb[n - 1]) < tolerance)
    {
        free(ab), free(bb), free(cb), free(db);
        return false;
    }
    x[n - 1] = db[n - 1] / bb[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        if (fabs(bb[i]) < tolerance)
        {
            free(ab), free(bb), free(cb), free(db);
            return false;
        }

        x[i] = (db[i] - cb[i] * x[i + 1]) / bb[i];
    }
    free(ab);
    free(bb);
    free(cb);
    free(db);
    return true;
}

bool solveLowerTriangular(double **L, double *x, double *b, int n)
{
    for (int i = 0; i < n; i++)
    {
        if (fabs(L[i][i]) < tolerance)
        {
            return false;
        }
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += L[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / L[i][i];
    }
    return true;
}

bool solveUpperTriangular(double **U, double *x, double *b, int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        if (fabs(U[i][i]) < tolerance)
        {
            return false;
        }
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += U[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / U[i][i];
    }
    return true;
}

bool gaussianElimination(double **A, double *x, double *b, int n)
{
    for (int k = 0; k < n - 1; k++)
    {
        if (fabs(A[k][k]) < tolerance)
        {
            return false;
        }
        for (int i = k + 1; i < n; i++)
        {
            double z = A[i][k] / A[k][k];
            for (int j = k; k < n; j++)
            {
                A[i][j] = A[i][j] - z * A[k][j];
            }
            b[i] = b[i] - z * b[k];
        }
    }
    return true;
}

bool LUDescomposition(double **A, double **L, double **U, int n)
{
    initMatrix(L, n, n, 0);
    initMatrix(U, n, n, 0);
    initDiagonal(U, n, 1);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L[i][k] * U[k][j];
            }
            L[i][j] = A[i][j] - sum;
        }
        for (int j = i + 1; j < n; j++)
        {
            if (fabs(L[i][i]) < tolerance)
            {
                return false;
            }
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = (A[i][j] - sum) / L[i][i];
        }
    }
    return true;
}

double calcError(double **A, double *x, double *b, int nr, int nc)
{
    double sum = 0;
    for (int i = 0; i < nr; i++)
    {
        double partialSum = 0;
        for (int j = 0; j < nc; j++)
        {
            partialSum += A[i][j] * x[j];
        }
        sum += pow(partialSum - b[i], 2);
    }
    return sqrt(sum);
}

double calcError(double *A, double *x, double *b, int nr)
{
    double sum = 0;
    for (int i = 0; i < nr; i++)
    {

        sum += pow(A[i] * x[i] - b[i], 2);
    }
    return sqrt(sum);
}

double calcError(double **A, double **L, double **U, int nr, int nc)
{
    double **R = createMatrix(nr, nc);
    multiplication(L, U, R, nr);
    double sum = 0;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
        {
            sum += pow(A[i][j] - R[i][j], 2);
        }
    }
    return sqrt(sum);
}

double calcError(double *x, double *s, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += pow(x[i] - s[i], 2);
    }
    return sqrt(sum);
}

double powerIteration(double **A, int n, int maxIterations, double tolerance)
{
    double *v = (double *)malloc(sizeof(double) * n);
    double *y = (double *)malloc(sizeof(double) * n);
    double *aux = (double *)malloc(sizeof(double) * n);
    double lamda;
    initVector(v, n, 1);
    int i;
    double error = 100;
    for (i = 0; i < maxIterations && error > tolerance; i++)
    {
        multiplication(A, v, y, n);
        normalizeVector(y, v, n);
        multiplication(v, A, aux, n);
        lamda = multiplication(aux, v, n);
        multiplication(lamda, v, aux, n);
        subtraction(y, aux, aux, n);
        error = getNorm2(aux, n);
    }
    printf("Iteraciones = %d Error = %f Lambda = %f\n", i, error, lamda);
    free(v);
    free(y);
    free(aux);
    return lamda;
}

void inverseIteration(double **A, double *x, double delta, int n, int maxIterations, double tolerance, double &mu, int &iteration, double &error)
{
    double ro;
    double *y = (double *)malloc(sizeof(double) * n);
    double *w = (double *)malloc(sizeof(double) * n);
    double *yLU = (double *)malloc(sizeof(double) * n);
    double **deltaM = identityMatrix(n);
    double **L = createMatrix(n, n);
    double **U = createMatrix(n, n);

    multiplication(delta, deltaM, deltaM, n);
    subtraction(A, deltaM, deltaM, n, n);
    LUDescomposition(deltaM, L, U, n);

    error = 100;
    for (iteration = 0; iteration < maxIterations && error > tolerance; iteration++)
    {
        solveLowerTriangular(L, yLU, x, n);
        solveUpperTriangular(U, y, yLU, n);
        normalizeVector(x, w, getNorm2(y, n), n);
        normalizeVector(y, x, n);
        ro = multiplication(x, w, n);
        mu = delta + ro;
        multiplication(ro, x, y, n);
        subtraction(w, y, y, n);
        error = getNorm2(y, n);
    }

    free(y);
    free(w);
    free(yLU);
    freeMatrix(deltaM);
    freeMatrix(L);
    freeMatrix(U);
}

void jacobiMethod(double **A, int n, int maxIterations, double tolerance, double **V, int &iteration, double &error)
{
    double t, c, s, delta;
    int iMax, jMax;
    double **G = createMatrix(n, n);
    double **GT = createMatrix(n, n);
    double **auxM = createMatrix(n, n);

    for (iteration = 0; iteration < maxIterations; iteration++)
    {
        error = -1;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (fabs(A[i][j]) > error)
                {
                    error = fabs(A[i][j]);
                    iMax = i;
                    jMax = j;
                }
            }
        }
        if (error < tolerance)
        {
            break;
        }
        delta = (A[jMax][jMax] - A[iMax][iMax]) / (2 * A[iMax][jMax]);
        t = sign(delta) / (fabs(delta) + sqrt(1 + delta * delta));
        c = 1 / sqrt(1 + t * t);
        s = c * t;
        setGivensMatrix(G, n, iMax, jMax, c, s);
        transposeMatrix(G, GT, n);
        multiplication(GT, A, auxM, n);
        multiplication(auxM, G, A, n);
        setValue(V, auxM, n, n);
        multiplication(auxM, G, V, n);
    }
    freeMatrix(G);
    freeMatrix(GT);
    freeMatrix(auxM);
}

int sign(double value)
{
    if (value > 0)
        return 1;
    else if (value < 0)
        return -1;
    else
        return 0;
}

void setGivensValues(double **G, int i, int j, double c, double s)
{
    G[i][i] = G[j][j] = c;
    G[i][j] = s;
    G[j][i] = -s;
}

void setGivensMatrix(double **G, int n, int i, int j, double c, double s)
{
    setIdentityMatrix(G, n);
    setGivensValues(G, i, j, c, s);
}

bool solveLU(double **A, double *x, double *b, int n)
{
    double *y = (double *)malloc(sizeof(double) * n);
    double **L = createMatrix(n, n);
    double **U = createMatrix(n, n);
    bool solved = true;

    initMatrix(L, n, n, 0);
    initMatrix(U, n, n, 0);
    initDiagonal(U, n, 1);

    if (LUDescomposition(A, L, U, n))
    {
        if (solveLowerTriangular(L, y, b, n))
        {
            if (!solveUpperTriangular(U, x, y, n))
            {
                solved = false;
            }
        }
        else
        {
            solved = false;
        }
    }
    else
    {
        solved = false;
    }
    free(y);
    freeMatrix(L);
    freeMatrix(U);
    return solved;
}

bool choleskyDecomposition(double **A, double **L, int n)
{
    setValue(L, 0.0, n, n);
    for (int j = 0; j < n; j++)
    {
        double sum = 0;
        for (int k = 0; k < j; k++)
        {
            sum += L[j][k] * L[j][k];
        }
        if ((L[j][j] = A[j][j] - sum) < 0)
        {
            return false;
        }
        if ((L[j][j] = sqrt(L[j][j])) < tolerance)
        {
            return false;
        }
        for (int i = j + 1; i < n; i++)
        {
            sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
        }
    }
    return true;
}

void conjugateGradient(double **A, double *x, double *b, int n, double tolerance, int &k, double &error)
{
    double *r = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *q = (double *)malloc(n * sizeof(double));
    double *aux = (double *)malloc(n * sizeof(double));
    double d, alpha, beta;
    multiplication(x, A, aux, n);
    subtraction(b, aux, r, n);
    copyVector(r, p, n);
    error = sqrt(multiplication(r, r, n) / n);
    for (k = 0; k < n && error >= tolerance; k++)
    {
        multiplication(p, A, q, n);
        d = multiplication(r, r, n);
        alpha = d / multiplication(p, q, n);
        multiplication(alpha, p, aux, n);
        addition(x, aux, x, n);
        multiplication(alpha, q, aux, n);
        subtraction(r, aux, r, n);
        beta = multiplication(r, r, n) / d;
        multiplication(beta, p, aux, n);
        addition(aux, r, p, n);
        error = sqrt(multiplication(r, r, n) / n);
    }
}

double **inverseMatrix(double **A, int n)
{
    double **aux = createMatrix(n, n);
    double **X = createMatrix(n, n);
    double **E = identityMatrix(n);
    for (int i = 0; i < n; i++)
        solveLU(A, aux[i], E[i], n);
    transposeMatrix(aux, X, n);
    freeMatrix(E);
    freeMatrix(aux);
    return X;
}

double conditionNumber(double **A, int n)
{
    double **AI = inverseMatrix(A, n);
    return getNormInf(A, n, n) * getNormInf(AI, n, n);
    freeMatrix(AI);
}

void aproxGradient(double *x, double *g, double (*f)(double *x, int n), double h, int n)
{
    double **matI = identityMatrix(n);
    double *xp = createVector(n);
    double *he = createVector(n);    
    for (int i = 0; i < n; i++)
    {
        multiplication(h, matI[i], he, n);
        addition(he, x, xp, n);
        g[i] = (f(xp, n) - f(x, n)) / h;
    }    
    freeMatrix(matI);
    free(xp);
    free(he);
}

void aproxHessian(double *x, double **H, double (*f)(double *x, int n), double h, int n)
{
    double **matI = identityMatrix(n);
    double *hei = createVector(n);
    double *hej = createVector(n);
    double *xp = createVector(n);
    for (int i = 0; i < n; i++)
    {
        multiplication(h, matI[i], hei, n);
        for (int j = 0; j < n; j++)
        {
            multiplication(h, matI[j], hej, n);
            addition(hei, x, xp, n);
            H[i][j] = -f(xp, n);
            addition(hej, x, xp, n);
            H[i][j] -= f(xp, n);
            addition(xp, hei, xp, n);
            H[i][j] += f(xp, n) + f(x, n);
            H[i][j] /= h * h;
        }
    }
    freeMatrix(matI);
    free(hei);
    free(hej);
    free(xp);
}

void aproxJacobian(double *x, double **J, void (*FR)(double *x, int n, int m, double *res), double h, int n, int m)
{
    double **matI = identityMatrix(n);
    double *xp = createVector(n);
    double *he = createVector(n);
    double *r1 = createVector(m);
    double *r2 = createVector(m);

    for (int j = 0; j < n; j++)
    {
        multiplication(h, matI[j], he, n);
        addition(he, x, xp, n);
        FR(x, n, m, r1);
        FR(xp, n, m, r2);
        for (int i = 0; i < m; i++)
        {
            J[i][j] = (r2[i] - r1[i]) / h;
        }
    }
    freeMatrix(matI);
    free(he);
    free(xp);
    free(r1);
    free(r2);
}
