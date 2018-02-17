#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include "gnuplot-iostream.h"
#include "plot.h"

void plotVol2d(SurfaceVol &SV, std::vector<SVIModel> &v_SVI) {
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:0.5]\n") % SV.getK(0) % SV.getK(SV.getNK() - 1);
    gp << "plot";
    for (int i(0); i < v_SVI.size(); ++i) {
        std::vector<std::pair<double, double> > pts_SVI;
        for (unsigned long j(0); j < SV.getNK(); ++j) {
            double k = SV.getK(j);
            double vol = v_SVI[i].getVol(k);
            pts_SVI.push_back(std::make_pair(k, vol));
        }
        gp << gp.file1d(pts_SVI) << boost::format("with lines lw 3 lc rgb \"red\" title 'SVI t%1%',")  %(i + 1);
    }

    for (unsigned long i(0); i < SV.getNT(); ++i) {
        std::vector<std::pair<double, double> > pts_SV;
        for (unsigned long j(0); j < SV.getNK(); ++j) {
            double k = SV.getK(j);
            pts_SV.push_back(std::make_pair(k, SV.getVol(i, j)));
        }
        if (i == SV.getNT() - 1) gp << gp.file1d(pts_SV) << "with lines lc rgb \"gray\" title ''" << std::endl;
        else gp << gp.file1d(pts_SV) << "with lines lc rgb \"gray\" title '',";
    }
}


void plotVar2d(SurfaceVol &SV, std::vector<SVIModel> &v_SVI){
    /////////////////////////////////////////////////////////////////////////////
    /////                   Check no calendar arbitrage                    //////
    /////////////////////////////////////////////////////////////////////////////
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:0.001]\n")  % SV.getK(0) % SV.getK(SV.getNK() - 1);
    gp << "plot";
    for(int i(0); i < v_SVI.size(); ++i){
        std::vector<std::pair<double, double> > pts_SVI;
        for(unsigned long j(0); j < SV.getNK(); ++j){
            double k = SV.getK(j);
            double w = v_SVI[i].getValue(k);
            pts_SVI.push_back(std::make_pair(k, w));
        }
        gp << gp.file1d(pts_SVI) << boost::format("with lines lw 3 lc rgb \"red\" title 'SVI t%1%',")  %(i + 1);
    }

    for(unsigned long i(0); i < SV.getNT(); ++i){
        std::vector<std::pair<double, double> > pts_SV;
        for(unsigned long j(0); j < SV.getNK(); ++j){
            pts_SV.push_back(std::make_pair(SV.getK(j), pow(SV.getVol(i,j), 2) * SV.getT(i)));
        }
        if (i == SV.getNT() - 1)    gp << gp.file1d(pts_SV) << "with lines lc rgb \"gray\" title ''" << std::endl;
        else                gp << gp.file1d(pts_SV) << "with lines lc rgb \"gray\" title '',";
    }
}


void plotVol3d(SurfaceVol &SV, std::vector<SVIModel> &v_SVI){
    unsigned long SV_Nk = SV.getNK();
    unsigned long SV_Nt = SV.getNT();

    Gnuplot gp;
    gp << "set ticslevel 0\n";
    gp << "set pm3d\n";
    gp << "set palette model RGB\n";
    gp << "set palette functions gray, gray, gray\n";
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [%3%:%4%]\n")  % SV.getK(0) % SV.getK(SV.getNK() - 1) % SV.getT(0) % SV.getT(SV.getNT() - 1);
    gp << "set zrange [0:0.6]\n";
    gp << "set hidden3d nooffset\n";
    gp << "set view 50,40\n";
    gp << "splot ";

    //SV
    std::vector<std::vector<double> > x_pts;
    std::vector<std::vector<double> > y_pts;
    std::vector<std::vector<double> > z_pts;
    x_pts.resize(SV_Nt);
    y_pts.resize(SV_Nt);
    z_pts.resize(SV_Nt);
    for(unsigned long i(0); i < SV_Nt; ++i){
        x_pts[i].resize(SV_Nk);
        y_pts[i].resize(SV_Nk);
        z_pts[i].resize(SV_Nk);
        for(unsigned long j(0); j < SV_Nk; ++j){
            x_pts[i][j] = SV.getK(j);
            y_pts[i][j] = SV.getT(i);
            z_pts[i][j] = SV.getVol(i,j);
        }
    }
    gp << gp.file2d(boost::make_tuple(x_pts, y_pts, z_pts)) << "with lines title 'Volatility Surface'" << std::endl;

}