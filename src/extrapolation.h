#ifndef AMDP_RED_EXTRAPOLATION_H
#define AMDP_RED_EXTRAPOLATION_H

#include <vector>
#include "interpolation.h"
#include "reverse-bs.h"
#include <string>
#include <iostream>

class Extrapolation {
public:
    Extrapolation(double &t, double (&range_k)[2], unsigned long Nk);
    Extrapolation(double &t, std::vector<double> &v_k);
    ~Extrapolation(){}

    void eval(std::string s, SVIModel SVI);

    std::vector<double> getK() const { return m_k;}
    std::vector<double> getVol() const { return m_vol;}

private:
    double m_t;
    double m_theta_t;

    std::vector<double> m_k;
    std::vector<double> m_vol;
};


#endif //AMDP_RED_EXTRAPOLATION_H
