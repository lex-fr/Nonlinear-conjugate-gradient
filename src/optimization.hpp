#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <algorithm>
#include <iostream>

#include "Mat.h"
#include "TestFunctions.hpp"
#include "linearAlgebra.hpp"

// Line Method
double getAlphaBacktrackingWolfeStrong(double *x, double fx, double *g, double *dk, double *aux1, double *aux2, int n, ObjetiveFunction &of, double alpha, double ro, double c1, double c2);

double getAlphaQuadraticWolfeStrong(double *x, double fi0, double *g, double *dk, double *aux1, double *aux2, int n, ObjetiveFunction &of, double alpha, double c1, double c2);

bool armijoCondition(double *x, double fiPrima0, double *g, double *aux, double fx, ObjetiveFunction &of, double alpha, double c1, int n);

bool strongWolfeConditions(double *x, double fiPrima0, double *dk, double *aux1, double *aux2, double fx, ObjetiveFunction &of, double alpha, double c1, double c2, int n);

bool armijoWolfeCondition(double *x, double fiPrima0, double *dk, double *aux, double fx, ObjetiveFunction &of, double alpha, double c1, int n);

bool curvatureStrongCondition(double *x, double fiPrima0, double *dk, double *aux1, double *aux2, double fx, ObjetiveFunction &of, double alpha, double c2, int n);

double getAlphaQuadratic(double *x, double fi0, double *g, double *aux, int n, ObjetiveFunction &of, double alpha, double c1);

double getAlphaCubic(double *x, double fi0, double *g, double *aux, int n, ObjetiveFunction &of, double alpha0, double c1);

double *linearConjugateGradient(double *x0, ObjetiveFunction of, int n, double tolg, double tolx, double tolf, const char *filename);

double *nonlinearConjugateGradient(std::ofstream &file, double *x0, ObjetiveFunction of, int n, double tolg, int maxitr, double c1, double c2, double ro,
                                   double (*calcBeta)(double *gk, double *gkp1, double *dk, double *aux, int n));

// Beta

double betaFletcherReeves(double *gk, double *gkp1, double *dk, double *aux, int n);

double betaPolakRibiere(double *gk, double *gkp1, double *dk, double *aux, int n);

double betaHestenesStiefel(double *gk, double *gkp1, double *dk, double *aux, int n);

double betaRMIL(double *gk, double *gkp1, double *dk, double *aux, int n);

double betaAMRI(double *gk, double *gkp1, double *dk, double *aux, int n);

#endif //OPTIMIZATION_H