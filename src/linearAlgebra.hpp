#ifndef SOLVE_ES_H
#define SOLVE_ES_H

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "Mat.h"

const double tolerance = 0.0000001;

bool solveDiagonal(double* d, double* x, double* b, int n);

bool solveTriDiagonal(double *a, double *b, double *c, double *d, double *x, int n);

bool solveLowerTriangular(double** L, double* x, double* b, int n);

bool solveUpperTriangular(double** U, double* x, double* b, int n);

bool gaussianElimination(double** A, double* x, double* b, int n);

bool LUDescomposition(double** A, double** L, double** U, int n);

double calcError(double **A, double *x, double *b, int nr, int nc);

double calcError(double *A, double *x, double *b, int nr);

double calcError(double **A, double **L, double **U, int nr, int nc);

double calcError(double *x, double *s, int n);

double powerIteration(double **A, int n, int maxIterations, double tolerance);

void inverseIteration(double **A, double *x, double delta, int n, int maxIterations, double tolerance, double &mu, int &iteration, double &error);

void jacobiMethod(double **A, int n, int maxIterations, double tolerance, double **V, int &iteration, double &error);

int sign(double value);

void setGivensValues(double **G, int i, int j, double c, double s);

void setGivensMatrix(double **G, int n, int i, int j, double c, double s);

bool solveLU(double **A, double *x, double *b, int n);

bool choleskyDecomposition(double** A, double** L, int n);

void conjugateGradient(double **A, double *x, double *b, int n, double tolerance, int &k, double &error);

double **inverseMatrix(double **A, int n);

double conditionNumber(double **A, int n);

void aproxGradient(double *x, double *g, double (*f)(double *x, int n), double h, int n);

void aproxHessian(double *x, double **H, double (*f)(double *x, int n), double h, int n);

void aproxJacobian(double *x, double **J, void (*FR)(double *x, int n, int m, double *res), double h, int n, int m);

#endif //SOLVE_ES_H