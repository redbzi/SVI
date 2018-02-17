#ifndef AMDP_RED_BLACK_SCHOLES_H
#define AMDP_RED_BLACK_SCHOLES_H

double norm_pdf(const double &x);
double norm_cdf(const double &x);
double d_j(const int &j, const double &S, const double &K, const double &r, const double &v, const double &T);
double call_price(const double &S, const double &K, const double &r, const double &v, const double &T);
double call_vega(const double &S, const double &K, const double &r, const double &v, const double &T);

#endif //AMDP_RED_BLACK_SCHOLES_H
