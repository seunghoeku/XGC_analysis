#include "output.hpp"
#include "adios2.h"
#include "sml.hpp"

// Example: 
// 2D array of shape (d0, d1): index at (i,j) is (d1*i + j)
// 3D array of shape (d0, d1, d2): index at (i,j,k) is ((d1+d2)*i + d2*j + k)
#define GET2D(X, d1, i, j) X[d1*i + j]
#define GET3D(X, d1, d2, i, j, k) X[(d1+d2)*i + d2*j + k]

adios2::IO output_io;
adios2::Engine writer;

extern Simulation sml;
extern adios2::ADIOS ad;

void output(HeatLoad &ion, HeatLoad &elec) {
    static bool first = true;

    if(first) {
        output_io = ad.DeclareIO("output");
        output_io.DefineVariable<double>("psi", {N_SIDE, N_PSI}, {0, 0}, {N_SIDE, N_PSI});
        output_io.DefineVariable<double>("ienflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("iptlflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("eenflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("eptlflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});

        output_io.DefineVariable<double>("io.side", {N_SIDE+1}, {0}, {N_SIDE+1});

        writer = output_io.Open("xgc.heatload.bp", adios2::Mode::Write);

        first = false;
    }

    // double psi[N_SIDE][N_PSI];
    std::vector<double> psi (N_SIDE*N_PSI);  //why vector? not just double?

    for(int is=0; is<N_SIDE; is++) {
        for(int i=0; i<N_PSI; i++){
           // psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
           GET2D(psi, N_PSI, is, i) = sml.pmin[is] + sml.dpsi[is] * (double)(i);
        }
    }

    /*double  ienflux[N_SIDE*N_COND*N_PSI];
    double iptlflux[N_SIDE*N_COND*N_PSI];
    double  eenflux[N_SIDE*N_COND*N_PSI];
    double eptlflux[N_SIDE*N_COND*N_PSI];*/

    std::vector<double>  ienflux(N_SIDE*N_COND*N_PSI);
    std::vector<double> iptlflux(N_SIDE*N_COND*N_PSI);
    std::vector<double>  eenflux(N_SIDE*N_COND*N_PSI);
    std::vector<double> eptlflux(N_SIDE*N_COND*N_PSI);

    for(int is=0; is<N_SIDE; is++) {
        for(int ic=0; ic<N_COND; ic++){
            for(int i=0; i<N_PSI; i++){
                // psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
                GET3D(ienflux, N_COND, N_PSI, is, ic, i) = ion.side[is].en[ic][i];
                GET3D(iptlflux, N_COND, N_PSI, is, ic, i) = ion.side[is].ptl[ic][i];
                GET3D(eenflux, N_COND, N_PSI, is, ic, i) = elec.side[is].en[ic][i];
                GET3D(eptlflux, N_COND, N_PSI, is, ic, i) = elec.side[is].ptl[ic][i];
            }
        }
    }


    // save psi, ion.side[0:N_SIDE].en[0:N_COND][0:N_PSI] and ptl[0:NCOND][0:N_PSI] and electron.
    auto var_psi = output_io.InquireVariable<double>("psi");
    auto var_ienflux =  output_io.InquireVariable<double>("ienflux");
    auto var_iptlflux = output_io.InquireVariable<double>("iptlflux");
    auto var_eenflux =  output_io.InquireVariable<double>("eenflux");
    auto var_eptlflux = output_io.InquireVariable<double>("eptlflux");


    writer.BeginStep();
    writer.Put<double>(var_psi, psi.data());
    writer.Put<double>(var_ienflux,  ienflux.data());
    writer.Put<double>(var_iptlflux, iptlflux.data());
    writer.Put<double>(var_eenflux,  eenflux.data());
    writer.Put<double>(var_eptlflux, eptlflux.data());

    writer.EndStep();
}

void output_finalize()
{
    writer.Close();
}