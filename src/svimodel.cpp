#include "svimodel.h"
#include <cmath>

SVIModel::SVIModel(double maturity, double a, double b, double ro, double m, double sig) {
    this->a = a;
    this->b = b;
    this->m = m;
    this->ro = ro;
    this->sig = sig;
    t = maturity;
}

double SVIModel::getValue(double x) const {
    return a + b*(ro*(x - m) + sqrt(std::pow(x - m, 2) + std::pow(sig, 2)));
}

double SVIModel::getVol(double x) const {
    double w = this->getValue(x);
    double t = this->getT();
    return sqrt(w / t);
}
