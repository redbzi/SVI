#ifndef AMDP_RED_SURFACEVOL_H
#define AMDP_RED_SURFACEVOL_H

#include <vector>
#include "svimodel.h"
#include "interpolation.h"
#include "extrapolation.h"

class SurfaceVol {
public:
    SurfaceVol(double (&range_k)[2], unsigned long Nk, double (&range_t)[2], unsigned long Nt);
    ~SurfaceVol(){}

    void build(std::vector<SVIModel> &v_SVI);

    double getVol(unsigned long i, unsigned long j) const {return m_vol[i][j];}
    double getK(unsigned long j) const {return m_k[j];}
    double getT(unsigned long i) const {return m_t[i];}

    unsigned long getNK() const {return m_Nk;}
    unsigned long getNT() const {return m_Nt;}

private:
    unsigned long                     m_Nk;
    unsigned long                     m_Nt;
    double                            m_range_k[2];
    std::vector<double>               m_k;
    std::vector<double>               m_t;
    std::vector<std::vector<double> > m_vol; //m_vol[i][j] for the value of slice at maturity m_t[i]// , strike m_k[j]

};


#endif //AMDP_RED_SURFACEVOL_H
