#ifndef TEST_FUNCTIONS_HPP
#define TEST_FUNCTIONS_HPP

#include <cmath>
#include "Mat.h"

struct ObjetiveFunction
{
    double (*f)(double *x, int n);
    void (*g)(double *x, double *g, int n);    
};

namespace tf
{

ObjetiveFunction getSixHumpCamel();
double sixHumpCamel(double *x, int n);
void gSixHumpCamel(double *x, double *g, int n);

ObjetiveFunction getBooth();
double booth(double *x, int n);
void gBooth(double *x, double *g, int n);

ObjetiveFunction getTreccani();
double treccani(double *x, int n);
void gTreccani(double *x, double *g, int n);

ObjetiveFunction getZettl();
double zettl(double *x, int n);
void gZettl(double *x, double *g, int n);

ObjetiveFunction getFreudensteinRoth();
double freudensteinRoth(double *x, int n);
void gFreudensteinRoth(double *x, double *g, int n);

ObjetiveFunction getHimmelblau();
double himmelblau(double *x, int n);
void gHimmelblau(double *x, double *g, int n);

ObjetiveFunction getRosembrock();
double rosembrock(double *x, int n);
void gRosembrock(double *x, double *g, int n);

ObjetiveFunction getWood();
double wood4D(double *x, int n);
void gWood4D(double *x, double *g, int n);

ObjetiveFunction getConvex1();
double convex1(double *x, int n);
void gConvex1(double *x, double *g, int n);

ObjetiveFunction getConvex2();
double convex2(double *x, int n);
void gConvex2(double *x, double *g, int n);


} // namespace tf

#endif // TEST_FUNCTIONS_HPP