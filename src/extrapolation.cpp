#include "extrapolation.h"
#include <cstdlib>
#include <cmath>
#include <algorithm>

Extrapolation::Extrapolation(double &t, double (&range_k)[2], unsigned long Nk){
    m_t = t;
    m_k.push_back(0.);
    double step_k = (range_k[1] - range_k[0]) / Nk;
    for(int i(0); i< Nk + 1; ++i){
        if (range_k[0] + step_k*i != 0.)   m_k.push_back(range_k[0] + step_k*i);
    }
    std::sort(m_k.begin(), m_k.end());
}

Extrapolation::Extrapolation(double &t, std::vector<double> &v_k){
    m_t = t;
    m_k = v_k;
}

void Extrapolation::eval(std::string s, SVIModel SVI){
    double t = SVI.getT();
    if (s=="DOWN"){
        if (m_t > t){
            std::cout << "bad input data in Extrapolation::eval DOWN" << std:: endl;
            exit(0);
        }
        Interpolation inter(m_t, m_k);
        inter.eval(SVI);
        m_k = inter.getK();
        m_vol = inter.getVol();
    }
    if (s=="UP"){
        if (m_t < t){
            std::cout << "bad input data in Extrapolation::eval UP" << std:: endl;
            exit(0);
        }
        std::vector<double> v_vol;
        for(int i(0); i < m_k.size(); ++i){
            v_vol.push_back(SVI.getVol(m_k[i]));
        }

        int index = -1;
        for(int i(0); i < m_k.size(); ++i){
            if (m_k[i] == 0.){
                index = i;
                break;
            }
        }
        if (index == -1) std::cout << "0. not found in vX" << std::endl;
        double theta_n = std::pow(v_vol[index], 2) * t;
        m_theta_t = (1 + m_t - t) * theta_n;

        for (int i(0); i < m_k.size(); ++i){
            m_vol.push_back(sqrt((std::pow(v_vol[i], 2) * t + m_theta_t - theta_n)/m_t));
        }

    }
}