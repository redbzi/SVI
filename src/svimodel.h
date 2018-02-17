#ifndef SVIMODEL_H
#define SVIMODEL_H

#include <vector>
#include <string>

class SVIModel
{
public:
    SVIModel(double maturity, double a, double b, double ro, double m, double sig);
    ~SVIModel(){}

    double              getT()             const {return t;}
    double              getValue(double x) const;
    double              getVol(double x)   const;

private:
    double a;
    double b;
    double ro;
    double m;
    double sig;
    double t;
};

#endif // SVIMODEL_H
