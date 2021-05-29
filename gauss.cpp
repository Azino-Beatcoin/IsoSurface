#include "gauss.h"
#include "pi.h"

double gaussFunction(
    double x, double sigma /* = 1. */, double mean /* = 0. */
) {
    double sigma2 = sigma*sigma;
    double xm2 = (x - mean)*(x - mean);
    return (
        (1./sqrt(2.*PI*sigma2))*
        exp(-xm2/(2*sigma2))
    );
}

void defineSmoothingWeights(
    double sigma,
    int& numWeights, std::vector<double>& weights
) {
    int sigma3 = int(ceil(3.*sigma));
    numWeights = 2*sigma3 + 1;
    weights.resize(numWeights);
    double s = 0.;
    for (int x = (-sigma3); x <= sigma3; ++x) {
        double g = gaussFunction(
            double(x), sigma
        );
        weights[x + sigma3] = g;
        s += g;
    }
    for (int i = 0; i < numWeights; ++i) {
        weights.at(i) /= s;
    }
}
