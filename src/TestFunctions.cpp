#include "TestFunctions.hpp"

namespace tf
{

ObjetiveFunction getSixHumpCamel()
{
    ObjetiveFunction of;
    of.f = sixHumpCamel;
    of.g = gSixHumpCamel;

    return of;
}

double sixHumpCamel(double *x, int n)
{
    return (4 - 2.1 * pow(x[0], 2) + pow(x[0], 4) / 3.0) * pow(x[0], 2) + x[0] * x[1] + (-4 + 4 * pow(x[1], 2)) * pow(x[1], 2);
}

void gSixHumpCamel(double *x, double *g, int n)
{
    g[0] = 2 * (pow(x[0], 5) - 4.2 * pow(x[0], 3) + 4 * x[0] + 0.5 * x[1]);
    g[1] = x[0] + 16 * pow(x[1], 3) - 8 * x[1];
}

ObjetiveFunction getBooth()
{
    ObjetiveFunction of;
    of.f = booth;
    of.g = gBooth;

    return of;
}

double booth(double *x, int n)
{
    // x* = (1,3) f(x*) = 0
    return pow(x[0] + 2 * x[1] - 7, 2) + pow(2 * x[0] + x[1] - 5, 2);
}

void gBooth(double *x, double *g, int n)
{
    g[0] = 2 * (x[0] + 2 * x[1] - 7) + 2 * 2 * (2 * x[0] + x[1] - 5);
    g[1] = 2 * 2 * (x[0] + 2 * x[1] - 7) + 2 * (2 * x[0] + x[1] - 5);
}

ObjetiveFunction getTreccani()
{
    ObjetiveFunction of;
    of.f = treccani;
    of.g = gTreccani;

    return of;
}

double treccani(double *x, int n)
{
    // x* = { (-2,0), (0,0) }   f(x*) = 0
    return pow(x[0], 4) + 4 * pow(x[1], 3) + 4 * pow(x[0], 2) + pow(x[1], 2);
}
void gTreccani(double *x, double *g, int n)
{
    g[0] = 4 * pow(x[0], 3) + 3 * 4 * pow(x[1], 2) + 2 * 4 * x[0];
    g[1] = 2 * x[1];
}

ObjetiveFunction getZettl()
{
    ObjetiveFunction of;
    of.f = zettl;
    of.g = gZettl;
    return of;
}

double zettl(double *x, int n)
{
    // x* = (-0.029896, 0)  f(x*) = -0.0037912
    return (1.0 / 4.0) * x[0] + pow(pow(x[0], 2) - 2 * x[0] + pow(x[1], 2), 2);
}
void gZettl(double *x, double *g, int n)
{
    g[0] = (1.0 / 4.0) + 2 * (pow(x[0], 2) - 2 * x[0] + pow(x[1], 2)) * (2 * x[0] - 2);
    g[1] = 2 * (pow(x[0], 2) - 2 * x[0] + pow(x[1], 2)) * (2 * x[1]);
}

ObjetiveFunction getFreudensteinRoth()
{
    ObjetiveFunction of;
    of.f = freudensteinRoth;
    of.g = gFreudensteinRoth;
    return of;
}
double freudensteinRoth(double *x, int n)
{
    return pow(x[0] - 13 + ((5 - x[1]) * x[1] - 2) * x[1], 2) + pow(x[0] - 29 + ((x[1] + 1) * x[1] - 14) * x[1], 2);
}
void gFreudensteinRoth(double *x, double *g, int n)
{
    g[0] = 4 * (x[0] + 3 * pow(x[1], 2) - 8 * x[1] - 21);
    g[1] = 4 * (x[0] * (6 * x[1] - 8) + 3 * pow(x[1], 5) - 10 * pow(x[1], 4) + 2 * pow(x[1], 3) - 60 * pow(x[1], 2) + 6 * x[1] + 216);
}

ObjetiveFunction getHimmelblau()
{
    ObjetiveFunction of;
    of.f = himmelblau;
    of.g = gHimmelblau;
    return of;
}

double himmelblau(double *x, int n)
{
    return pow(pow(x[0], 2) + x[1] - 11, 2) + pow(x[0] + pow(x[1], 2) - 7, 2);
}

void gHimmelblau(double *x, double *g, int n)
{
    g[0] = 2 * (2 * x[0] * (pow(x[0], 2) + x[1] - 11) + x[0] + pow(x[1], 2) - 7);
    g[1] = 2 * (pow(x[0], 2) + 2 * x[1] * (x[0] + pow(x[1], 2) - 7) + x[1] - 11);
}

ObjetiveFunction getWood()
{
    ObjetiveFunction of;
    of.f = wood4D;
    of.g = gWood4D;
    return of;
}

double wood4D(double *x, int n)
{
    return 100 * pow(pow(x[0], 2) - x[1], 2) + pow(x[0] - 1, 2) + pow(x[2] - 1, 2) + 90 * pow(pow(x[2], 2) - x[3], 2) + 10.1 * (pow(x[1] - 1, 2) + pow(x[3] - 1, 2)) + 19.8 * (x[1] - 1) * (x[3] - 1);
}

void gWood4D(double *x, double *g, int n)
{
    g[0] = 400 * x[0] * (pow(x[0], 2) - x[1]) + 2 * (x[0] - 1);
    g[1] = -200 * (pow(x[0], 2) - x[1]) + 20.2 * (x[1] - 1) + 19.8 * (x[3] - 1);
    g[2] = 2 * (x[2] - 1) + 360 * x[2] * (pow(x[2], 2) - x[3]);
    g[3] = -180 * (pow(x[2], 2) - x[3]) + 20.2 * (x[3] - 1) + 19.8 * (x[1] - 1);
}

ObjetiveFunction getRosembrock()
{
    ObjetiveFunction of;
    of.f = rosembrock;
    of.g = gRosembrock;
    return of;
}

double rosembrock(double *x, int n)
{
    double result = 0;
    for (int i = 0; i < n - 1; i++)
    {
        result += 100 * pow(x[i + 1] - pow(x[i], 2), 2) + pow(1 - x[i], 2);
    }
    return result;
}

void gRosembrock(double *x, double *g, int n)
{
    g[0] = -400 * x[0] * (x[1] - pow(x[0], 2)) - 2 + 2 * x[0];
    g[n - 1] = 200 * (x[n - 1] - pow(x[n - 2], 2));
    for (int i = 1; i < n - 1; i++)
    {
        g[i] = 200 * (x[i] - pow(x[i - 1], 2)) - 400 * x[i] * (x[i + 1] - pow(x[i], 2)) - 2 + 2 * x[i];
    }
}

ObjetiveFunction getConvex1()
{
    ObjetiveFunction of;
    of.f = convex1;
    of.g = gConvex1;
    return of;
}

double convex1(double *x, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res += exp(x[i]) - x[i];
    }
    return res;
}

void gConvex1(double *x, double *g, int n)
{
    for (int i = 0; i < n; i++)
    {
        g[i] = exp(x[i]) - 1;
    }
}

ObjetiveFunction getConvex2()
{
    ObjetiveFunction of;
    of.f = convex2;
    of.g = gConvex2;
    return of;
}

double convex2(double *x, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res += ((i + 1.) / n) * (exp(x[i]) - x[i]);
    }
    return res;
}

void gConvex2(double *x, double *g, int n)
{
    for (int i = 0; i < n; i++)
    {
        g[i] = ((i + 1.) / n) * (exp(x[i]) - 1);
    }
}

} // namespace tf
