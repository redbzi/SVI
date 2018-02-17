#include "reverse-bs.h"
#include "black-scholes.h"
#include <cmath>
#include <iostream>


double solve2(double &S, double &K, double &r, double &T, double &C_target){
    //in case the Newton-Raphson method fails, we use dichotomy
    double n = 0.;
    double m = 1.;
    double epsilon = 0.01;
    double x = 0.5 * (m + n);
    double y = call_price(S, K, r, x, T);

    do {
        if (y < C_target) {
            m = x;
        }

        if (y > C_target) {
            n = x;
        }

        x = 0.5 * (m + n);
        y = call_price(S, K, r, x, T);
    } while (fabs(y-C_target) > epsilon);

    return x;
}

double solve_vol(double &S, double &K, double &r, double &T, double &C_target){
    //main algorithm to reverse BS
    double x = 1.;
    double epsilon = 0.001;
    double C = call_price(S, K, r, x, T);

    while (std::abs(C - C_target) > epsilon) {
        double V = call_vega(S, K, r, x, T);
        C = call_price(S, K, r, x, T);
        x = x - (C - C_target)/V;
        if (x < 0.) break;
    }
    if (x < 0.){
        x = solve2(S, K, r, T, C_target);
    }
    return x;
}

