#ifndef GAUSS_H
#define GAUSS_H

#include <vector>

// Gauss Bell Curve
double gaussFunction(double x, double sigma = 1., double mean = 0.);

void defineSmoothingWeights(
    double sigma,
    int& numWeights, std::vector<double>& weights
);

#endif
