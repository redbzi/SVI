#include "interpolation.h"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace Red {
    double S = 100.;
    double r = 0.025;
}

Interpolation::Interpolation(double &t, double (&range_k)[2], unsigned long Nk){
    m_t = t;
    m_k.push_back(0.);
    double step_k = (range_k[1] - range_k[0]) / Nk;
    for(int i(0); i< Nk + 1; ++i){
        if (range_k[0] + step_k*i != 0.)   m_k.push_back(range_k[0] + step_k*i);
    }
    std::sort(m_k.begin(), m_k.end());
}

Interpolation::Interpolation(double &t, std::vector<double> &v_k){
    m_t = t;
    m_k = v_k;
}

void Interpolation::eval(SVIModel &SVI_1, SVIModel &SVI_2){
    //forwards
    double t_1 = SVI_1.getT();
    double t_2 = SVI_2.getT();
    double F_1 = Red::S * exp(Red::r * t_1);
    double F_t = Red::S * exp(Red::r * m_t);
    double F_2 = Red::S * exp(Red::r * t_2);

    //evaluate vol given SVI_1 and SVI_2 over m_k
    std::vector<double> vY_1;
    std::vector<double> vY_2;
    for (int i(0); i < m_k.size(); ++i){
        vY_1.push_back(SVI_1.getVol(m_k[i]));
        vY_2.push_back(SVI_2.getVol(m_k[i]));
    }

    //Black-Scholes to compute call prices for every x in the common x-vector m_k
    std::vector<double> vC_1;
    std::vector<double> vK_1;
    std::vector<double> vC_2;
    std::vector<double> vK_2;

    for(int i(0); i < m_k.size(); ++i){
        //BS for the first slice t_1
        vK_1.push_back(F_1 * exp(m_k[i]));
        vC_1.push_back(call_price(Red::S, vK_1[i], Red::r, vY_1[i], t_1));

        //BS for the second slice t_2
        vK_2.push_back(F_2 * exp(m_k[i]));
        vC_2.push_back(call_price(Red::S, vK_2[i], Red::r, vY_2[i], t_2));
    }

    int index = -1;
    for(int i(0); i < m_k.size(); ++i){
        if (m_k[i] == 0.){
            index = i;
            break;
        }
    }

    if (index == -1) std::cout << "0. not found in m_k" << std::endl;

    double theta_1 = std::pow(vY_1[index], 2) * t_1;
    double theta_2 = std::pow(vY_2[index], 2) * t_2;
    m_theta_t = ((t_2 - m_t) * theta_1 + (m_t - t_1) * theta_2) / (t_2 - t_1);
    m_alpha_t = (sqrt(theta_2) - sqrt(m_theta_t)) / (sqrt(theta_2) - sqrt(theta_1));

    //interpolate call prices
    std::vector<double> vC;
    std::vector<double> vK;

    for(int i(0); i < m_k.size(); ++i){
        vK.push_back(F_t * exp(m_k[i]));
        vC.push_back(vK[i] * (m_alpha_t * vC_1[i] / vK_1[i] + (1 - m_alpha_t) * vC_2[i] / vK_2[i]));
    }

    for (int i(0); i < m_k.size(); ++i){
        m_vol.push_back(solve_vol(Red::S, vK[i], Red::r, m_t, vC[i]));
    }
}


void Interpolation::eval(SVIModel &SVI_2){
    //overloading of eval for the extrapolation down (we use t=0 and t=SVI_2.getT())
    //forwards
    double t_1 = 0.;
    double t_2 = SVI_2.getT();
    double F_1 = Red::S * exp(Red::r * t_1);
    double F_t = Red::S * exp(Red::r * m_t);
    double F_2 = Red::S * exp(Red::r * t_2);

    //evaluate vol given by SVI_2 over m_k
    std::vector<double> vY_2;
    for(int i(0); i < m_k.size(); ++i){
        vY_2.push_back(SVI_2.getVol(m_k[i]));
    }

    //Black-Scholes to compute call prices for every x in the common x-vector m_k
    std::vector<double> vC_1;
    std::vector<double> vK_1;
    std::vector<double> vC_2;
    std::vector<double> vK_2;

    for(int i(0); i < m_k.size(); ++i){
        //BS for the first slice t_1
        vK_1.push_back(F_1 * exp(m_k[i]));
        vC_1.push_back(std::max(Red::S - vK_1[i], 0.));

        //BS for the second slice t_2
        vK_2.push_back(F_2 * exp(m_k[i]));
        vC_2.push_back(call_price(Red::S, vK_2[i], Red::r, vY_2[i], t_2));
    }

    int index = -1;
    for(int i(0); i < m_k.size(); ++i){
        if (m_k[i] == 0.){
            index = i;
            break;
        }
    }

    if (index == -1) std::cout << "0. not found in m_k" << std::endl;

    double theta_1 = 0.;
    double theta_2 = std::pow(vY_2[index], 2) * t_2;
    m_theta_t = ((t_2 - m_t) * theta_1 + (m_t - t_1) * theta_2) / (t_2 - t_1);
    m_alpha_t = (sqrt(theta_2) - sqrt(m_theta_t)) / (sqrt(theta_2) - sqrt(theta_1));

    //interpolate call prices
    std::vector<double> vC;
    std::vector<double> vK;


    for(int i(0); i < m_k.size(); ++i){
        vK.push_back(F_t * exp(m_k[i]));
        vC.push_back(vK[i] * (m_alpha_t * vC_1[i] / vK_1[i] + (1 - m_alpha_t) * vC_2[i] / vK_2[i]));
    }

    for (int i(0); i < m_k.size(); ++i){
        m_vol.push_back(solve_vol(Red::S, vK[i], Red::r, m_t, vC[i]));
    }
}