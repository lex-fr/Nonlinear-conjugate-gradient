#ifndef MAT_H
#define MAT_H

#include <string>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cfloat>

void initVector(double *vec, int n, double value);

void initMatrix(double **mat, int nr, int nc, double value);

double sumVector(double *vec, int n);

double getNorm2(double *vec, int n);

double maximumNorm(double *vec, int n);

double getNorm1(double **mat, int nr, int nc);

double getNorm2(double **mat, int nr, int nc);

double getNormF(double **mat, int nr, int nc);

double getNormInf(double **mat, int nr, int nc);

double maximumNorm(double **mat, int nr, int nc);

bool isSymmetric(double **mat, int n);

bool isSymmetric(double **mat, double **tranposeMat, int n);

void normalizeVector(double *vec, int n);

void normalizeVector(double *src, double *dst, int n);

void normalizeVector(double *vec, double norm, int n);

void normalizeVector(double *src, double *dst, double norm, int n);

void initDiagonal(double **mat, int nr, double value);

void additionDiagonal(double a, double **B, double **R, int n);

void addition(double **A, double **B, double **R, int nr, int nc);

void addition(double *a, double *b, double *r, int n);

void subtraction(double *a, double *b, double *r, int n);

void subtraction(double **A, double **B, double **R, int nr, int nc);

void multiplication(double **P, double **Q, double **R, int nr1, int nc1, int nc2);

void multiplication(double **P, double **Q, double **R, int n);

void multiplication(double **A, double *b, double *r, int n);

void multiplication(double **A, double *b, double *r, int nr, int nc);

void multiplication(double *b, double **A, double *r, int n);

void multiplication(double *a, double *b, double *r, int n);

void multiplication(double *a, double *b, double **R, int n);

double multiplication(double *a, double *b, int n);

void multiplication(double a, double *b, double *r, int n);

void multiplication(double a, double **B, double **R, int n);

void multiplication(double a, double **B, double **R, int nr, int nc);

void divition(double a, double **B, double **R, int nr, int nc);

void divition(double a, double *b, double *r, int n);

void printMatrix(double **mat, int nr, int nc);

void printVector(double *vec, int n);

void printDiagonal(double **mat, int n);

double **createMatrix(int nr, int nc);

double **createMatrix(int n);

double **createMatrix(double value, int nr, int nc);

double **createMatrix(double *a, double *b, double *c, int nr, int nc);

void freeMatrix(double **mat);

double **readMatrix(char *cfile, int *nr, int *nc);

void matToFile(char *filename, double **mat, int nr, int nc);

void vectorToFile(const char *filename, double *vec, int n);

double *createVector(int n);

double *readVector(char *cfile, int *nr);

double *readVector(char *cfile, int n);

double *readTxtVector(const char *cfile, int &n);

void readPoints(char *cfile, double *x, double *y, int *n);

double **identityMatrix(int n);

void setIdentityMatrix(double **I, int n);

void setValue(double **A, double value, int nr, int nc);

void setValue(double **src, double **dst, int nr, int nc);

void copyDiagonal(double **src, double **dst, int n);

void copyVector(double *src, double *dst, int n);

void setDiagonal(double **A, double value, int n);

void transposeMatrix(double **src, double **dst, int nr, int nc);

void transposeMatrix(double **src, double **dst, int n);

double **getVandermonde(int n);

double **getVandermonde(double *x, int n);

double polynomialFunction(double x, double *c, int n);

void createPositiveDefiniteMatrix(double **mat, int n);

#endif // MAT_H