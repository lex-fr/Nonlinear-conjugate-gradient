#include <iostream>
#include <vector>
#include <chrono>

#include "optimization.hpp"
#include "TestFunctions.hpp"

//Config
double tolg = 1e-6;
int maxitr = 30000;
double c1 = 1e-4;
double c2 = 0.1;
double ro = 0.99;

void evalFunction(std::string filename, ObjetiveFunction of, std::vector<std::string> dim, std::vector<std::string> points, std::vector<std::pair<std::string, double (*)(double *, double *, double *, double *, int)>> betas)
{

    std::string dir = "results/";
    for (auto d : dim)
    {
        for (auto p : points)
        {
            std::ofstream file(dir + filename + p + "_" + d + "d.txt");
            file << "$\\beta$ & k & $f(x_k)$ & $||g(x_k)||$ & ms \\\\ \\hline" << std::endl;
            for (auto beta : betas)
            {
                int n;
                double *x0 = readTxtVector(("data/" + p + "_" + d + "d.txt").c_str(), n);
                file << beta.first << " & ";
                auto start = std::chrono::steady_clock::now();
                double *x = nonlinearConjugateGradient(file, x0, of, n, tolg, maxitr, c1, c2, ro, beta.second);
                file << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count() << " \\\\ \\hline" << std::endl;
                free(x0);
                free(x);
            }
            file.close();
        }
    }
}

int main(int argc, char **argv)
{

    std::vector<std::pair<std::string, double (*)(double *, double *, double *, double *, int)>> betas = {{"FR", betaFletcherReeves}, {"PRP", betaPolakRibiere}, {"RMIL", betaRMIL}, {"AMRI", betaAMRI}};

    std::vector<std::string> dim = {"2"};
    std::vector<std::string> points = {"8", "m8", "10", "m10"};
    evalFunction("SixHumpCamel", tf::getSixHumpCamel(), dim, points, betas);
    points = {"10", "25", "50", "100"};
    evalFunction("Booth", tf::getBooth(), dim, points, betas);
    points = {"5", "10", "50", "100"};
    evalFunction("Treccani", tf::getTreccani(), dim, points, betas);
    points = {"5", "10", "20", "50"};
    evalFunction("Zettl", tf::getZettl(), dim, points, betas);
    points = {"2", "10", "20", "50"};
    evalFunction("FreudensteinRoth", tf::getFreudensteinRoth(), dim, points, betas);
    points = {"10", "50", "100", "200"};
    evalFunction("Himmelblau", tf::getHimmelblau(), dim, points, betas);

    dim = {"4"};
    points = {"13", "16", "20", "30"};
    evalFunction("Wood", tf::getWood(), dim, points, betas);
    dim = {"2", "4", "10", "100", "500", "1000", "10000"};
    evalFunction("Rosembrock", tf::getRosembrock(), dim, points, betas);
    evalFunction("Convex1", tf::getConvex1(), dim, points, betas);
    evalFunction("Convex2", tf::getConvex2(), dim, points, betas);

    return 0;
}