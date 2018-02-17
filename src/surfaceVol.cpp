#include "surfaceVol.h"
#include <algorithm>


SurfaceVol::SurfaceVol(double (&range_k)[2], unsigned long Nk, double (&range_t)[2], unsigned long Nt){
    //discretization in k
    m_k.push_back(0.);
    double step_k = (range_k[1] - range_k[0]) / Nk;
    for(int i(0); i< Nk + 1; ++i){
        if (range_k[0] + step_k*i != 0.)   m_k.push_back(range_k[0] + step_k*i);
    }
    std::sort(m_k.begin(), m_k.end());
    m_Nk = m_k.size();

    m_range_k[0] = m_k[0];
    m_range_k[1] = m_k[m_Nk - 1];

    //discretization in T
    double step_t = (range_t[1] - range_t[0]) / Nt;
    for(int i(0); i< Nt + 1; ++i){
        m_t.push_back(range_t[0] + step_t*i);
    }
    m_Nt = m_t.size();

    //allocate memory space
    m_vol.resize(m_Nt);
    for(int i(0); i < m_vol.size(); ++i){
        m_vol[i].resize(m_Nk);
    }
}

void SurfaceVol::build(std::vector<SVIModel> &v_SVI){ //v_SVI must be sorted from the shortest maturity to the largest one
    //[0, inf - 1] define the interval of indeces where "extrapolation down" needs to be done
    //[inf, sup] "interpolation"
    //[sup + 1, m_Nt - 1] "extrapolation up"

    //////////// Extrapolation Down ////////////
    long inf = 0;
    while(m_t[inf] <= v_SVI[0].getT()){
        inf++;
    }
    if(inf > 0)
    {
        for(int j(0); j < inf; ++j)
        {
            Extrapolation ext(m_t[j], m_k);
            ext.eval("DOWN", v_SVI[0]);
            m_vol[j] = ext.getVol();
        }
    }

    //////////// Extrapolation Up ////////////
    long sup = m_Nt - 1;
    while(v_SVI[v_SVI.size() - 1].getT() <= m_t[sup]){
        sup--;
    }
    if(sup < (long) m_Nt - 1)
    {
        for(long j(sup + 1); j < m_Nt; ++j)
        {
            Extrapolation ext(m_t[j], m_k);
            ext.eval("UP", v_SVI[v_SVI.size() - 1]);
            m_vol[j] = ext.getVol();
        }
    }

    //////// Multiple Interpolation /////////
    if(v_SVI.size() >= 2) //if there are some maturities inside the SVIs range
    {
        long mid = inf;
        for(int i(0); i < v_SVI.size() - 1; ++i){ //along SVI slice intervals A_i = [ SVI[i], SVI[i+1] ]
            long index = sup;
            while(v_SVI[i+1].getT() < m_t[index])
            {
                index--;
                if (index == mid - 1) break; //if true, this means we are in an area where slices have already been computed
            }
            for(long j(mid); j < index + 1; ++j){ //along the remaining slices to find
                Interpolation inter(m_t[j], m_k);
                inter.eval(v_SVI[i], v_SVI[i+1]);
                m_vol[j] = inter.getVol();
                mid++;
            }
        }
    }
}