#ifndef AMDP_RED_INTERPOLATION_H
#define AMDP_RED_INTERPOLATION_H

#include <vector>
#include "svimodel.h"
#include "black-scholes.h"
#include "reverse-bs.h"

class Interpolation {
public:
    Interpolation(double &t, double (&range_k)[2], unsigned long Nk);
    Interpolation(double &t, std::vector<double> &v_k);
    ~Interpolation(){}

    void eval(SVIModel &SVI_1, SVIModel &SVI_2);
    void eval(SVIModel &SVI_2); //lower extrapolation (uses the maturity T = 0)

    std::vector<double> getK() const { return m_k;}
    std::vector<double> getVol() const { return m_vol;}

private:
    double m_t;
    double m_theta_t;
    double m_alpha_t;

    std::vector<double> m_k;
    std::vector<double> m_vol;
};


#endif //AMDP_RED_INTERPOLATION_H
