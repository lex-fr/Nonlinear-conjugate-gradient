#include "optimization.hpp"

double getAlphaBacktrackingWolfeStrong(double *x, double fx, double *g, double *dk, double *aux1, double *aux2, int n, ObjetiveFunction &of, double alpha, double ro, double c1, double c2)
{
    double fiPrima0 = multiplication(g, dk, n);
    while (!strongWolfeConditions(x, fiPrima0, dk, aux1, aux2, fx, of, alpha, c1, c2, n) && alpha > sqrt(DBL_EPSILON))
    {
        alpha *= ro;
    }
    return alpha;
}

double getAlphaQuadraticWolfeStrong(double *x, double fi0, double *g, double *dk, double *aux1, double *aux2, int n, ObjetiveFunction &of, double alpha, double c1, double c2)
{
    double fiPrima0 = multiplication(g, dk, n);
    while (!strongWolfeConditions(x, fiPrima0, dk, aux1, aux2, fi0, of, alpha, c1, c2, n) && alpha > sqrt(DBL_EPSILON))
    {
        multiplication(alpha, dk, aux1, n);
        addition(x, aux1, aux1, n);
        alpha = (-pow(alpha, 2) * fiPrima0) / fabs(2. * (of.f(aux1, n) - fi0 - fiPrima0 * alpha));
    }
    return alpha;
}

bool strongWolfeConditions(double *x, double fiPrima0, double *dk, double *aux1, double *aux2, double fx, ObjetiveFunction &of, double alpha, double c1, double c2, int n)
{
    //printf("Wolfe - alpha: %e \n", alpha);
    return armijoWolfeCondition(x, fiPrima0, dk, aux1, fx, of, alpha, c1, n) && curvatureStrongCondition(x, fiPrima0, dk, aux1, aux2, fx, of, alpha, c2, n);
}

bool armijoWolfeCondition(double *x, double fiPrima0, double *dk, double *aux, double fx, ObjetiveFunction &of, double alpha, double c1, int n)
{
    multiplication(alpha, dk, aux, n);
    addition(x, aux, aux, n);
    //printf("Amr - 1: %e 2: %e cumple: %d \n", of.f(aux, n), fx + c1 * alpha * fiPrima0, of.f(aux, n) <= fx + c1 * alpha * fiPrima0);
    return of.f(aux, n) <= fx + c1 * alpha * fiPrima0;
}

bool curvatureStrongCondition(double *x, double fiPrima0, double *dk, double *aux1, double *aux2, double fx, ObjetiveFunction &of, double alpha, double c2, int n)
{
    multiplication(alpha, dk, aux1, n);
    addition(x, aux1, aux1, n);
    of.g(aux1, aux2, n);
    //printf("Cur - 1: %e 2: %e cumple: %d \n", fabs(multiplication(aux2, dk, n)), -c2 * fiPrima0, fabs(multiplication(aux2, dk, n)) <= -c2 * fiPrima0);
    return fabs(multiplication(aux2, dk, n)) <= -c2 * fiPrima0;
}

double getAlphaQuadratic(double *x, double fi0, double *g, double *aux, int n, ObjetiveFunction &of, double alpha, double c1)
{
    double fiPrima0 = -multiplication(g, g, n);
    double fiAlpha;
    multiplication(-alpha, g, aux, n);
    addition(x, aux, aux, n);
    fiAlpha = of.f(aux, n);
    while (fiAlpha > fi0 + c1 * alpha * fiPrima0)
    {
        alpha = (-pow(alpha, 2) * fiPrima0) / (2 * (fiAlpha - fi0 - fiPrima0 * alpha));
        multiplication(-alpha, g, aux, n);
        addition(x, aux, aux, n);
        fiAlpha = of.f(aux, n);
    }
    return alpha;
}

double getAlphaCubic(double *x, double fi0, double *g, double *aux, int n, ObjetiveFunction &of, double alpha0, double c1)
{
    double tol = sqrt(DBL_EPSILON);
    double alpha1, fiAlpha0, fiAlpha1, fiPrima0 = -multiplication(g, g, n);
    multiplication(-alpha0, g, aux, n);
    addition(x, aux, aux, n);
    fiAlpha0 = of.f(aux, n);
    if (fiAlpha0 <= fi0 + alpha0 * c1 * fiPrima0)
    {
        return alpha0;
    }
    alpha1 = (-pow(alpha0, 2) * fiPrima0) / (2 * (fiAlpha0 - fi0 - fiPrima0 * alpha0));
    multiplication(-alpha1, g, aux, n);
    addition(x, aux, aux, n);
    fiAlpha1 = of.f(aux, n);
    double a, b, c = fiPrima0;
    while (fabs(alpha0 - alpha1) < tol && fiAlpha1 > fi0 + alpha1 * c1 * fiPrima0)
    {
        multiplication(-alpha0, g, aux, n);
        addition(x, aux, aux, n);
        fiAlpha0 = of.f(aux, n);
        a = pow(alpha1, 2) * pow(alpha0, 2) * alpha1 - alpha0 * pow(alpha1, 2) * pow(alpha0, 2);
        b = (pow(alpha0, 3) * (fi0 - fiAlpha1 + fiPrima0 * alpha1) + pow(alpha1, 3) * (fiAlpha0 - fi0 - fiPrima0 * alpha0)) / a;
        a = (pow(alpha0, 2) * (fiAlpha1 - fi0 - fiPrima0 * alpha1) + pow(alpha1, 2) * (fi0 - fiAlpha0 + fiPrima0 * alpha0)) / a;
        alpha0 = alpha1;
        alpha1 = (sqrt(pow(b, 2) - 3 * a * c) - b) / (3 * a);
        multiplication(-alpha1, g, aux, n);
        addition(x, aux, aux, n);
        fiAlpha1 = of.f(aux, n);
    }
    return alpha1;
}

double *nonlinearConjugateGradient(std::ofstream &file, double *x0, ObjetiveFunction of, int n, double tolg, int maxitr, double c1, double c2, double ro,
                                   double (*calcBeta)(double *gk, double *gkp1, double *dk, double *aux, int n))
{
    int i = 0;
    double *xk = createVector(n);
    double *gk = createVector(n);
    double *gkp1 = createVector(n);
    double *dk = createVector(n);
    double *aux1 = createVector(n);
    double *aux2 = createVector(n);
    double alpha, beta, fk;
    copyVector(x0, xk, n);
    of.g(x0, gk, n);
    multiplication(-1, gk, dk, n);
    do
    {
        fk = of.f(xk, n);
        alpha = getAlphaBacktrackingWolfeStrong(xk, fk, gk, dk, aux1, aux2, n, of, 5, ro, c1, c2);        
        multiplication(alpha, dk, aux1, n);
        addition(xk, aux1, xk, n);
        of.g(xk, gkp1, n);
        beta = calcBeta(gk, gkp1, dk, aux1, n);
        multiplication(beta, dk, aux1, n);
        subtraction(aux1, gkp1, dk, n);
        copyVector(gkp1, gk, n);
        i++;
    } while (getNorm2(gkp1, n) >= tolg && i < maxitr);
    file << i << " & " << of.f(xk, n) << " & " << getNorm2(gkp1, n) << " & ";
    free(gk);
    free(gkp1);
    free(dk);
    free(aux1);
    free(aux2);
    return xk;
}

double betaFletcherReeves(double *gk, double *gkp1, double *dk, double *aux, int n)
{
    return multiplication(gkp1, gkp1, n) / multiplication(gk, gk, n);
}

double betaPolakRibiere(double *gk, double *gkp1, double *dk, double *aux, int n)
{
    subtraction(gkp1, gk, aux, n);
    return std::max(multiplication(gkp1, aux, n) / multiplication(gk, gk, n), 0.0);
}

double betaHestenesStiefel(double *gk, double *gkp1, double *dk, double *aux, int n)
{
    subtraction(gkp1, gk, aux, n);
    return multiplication(gkp1, aux, n) / multiplication(aux, dk, n);
}

double betaRMIL(double *gk, double *gkp1, double *dk, double *aux, int n)
{
    subtraction(gkp1, gk, aux, n);
    return std::max(multiplication(gkp1, aux, n) / multiplication(dk, dk, n), 0.0);
}

double betaAMRI(double *gk, double *gkp1, double *dk, double *aux, int n)
{
    multiplication(getNorm2(gkp1, n) / getNorm2(gk, n), gk, aux, n);
    subtraction(gkp1, aux, aux, n);
    return std::max(multiplication(gkp1, aux, n) / multiplication(dk, dk, n), 0.0);
}
