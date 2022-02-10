#include "output.hpp"
#include "adios2.h"
#include "sml.hpp"

// Example: 
// 2D array of shape (d0, d1): index at (i,j) is (d1*i + j)
// 3D array of shape (d0, d1, d2): index at (i,j,k) is ((d1+d2)*i + d2*j + k)
#define GET2D(X, d1, i, j) X[d1*i + j]
#define GET3D(X, d1, d2, i, j, k) X[(d1+d2)*i + d2*j + k]

adios2::ADIOS ad2;
adios2::IO output_io;
adios2::Engine writer;

extern Simulation sml;

void output(HeatLoad &ion, HeatLoad &elec) {
    static bool first = true;

    if(first) {
        output_io = ad2.DeclareIO("output");
        output_io.DefineVariable<double>("psi", {N_SIDE, N_PSI}, {0, 0}, {N_SIDE, N_PSI});
        output_io.DefineVariable<double>("io.side", {N_SIDE+1}, {0}, {N_SIDE+1});

        writer = output_io.Open("xgc.heatload.bp", adios2::Mode::Write);
        writer.BeginStep();

        first = false;
    }

    // double psi[N_SIDE][N_PSI];
    std::vector<double> psi (N_SIDE*N_PSI);

    for(int is=0; is<N_SIDE; is++) {
        for(int i=0; i<N_PSI; i++){
           // psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
           GET2D(psi, N_PSI, is, i) = sml.pmin[is] + sml.dpsi[is] * (double)(i);
        }
    }

    // save psi, ion.side[0:N_SIDE].en[0:N_COND][0:N_PSI] and ptl[0:NCOND][0:N_PSI] and electron.
    auto var_psi = output_io.InquireVariable<double>("psi");

    writer.BeginStep();
    writer.Put<double>(var_psi, psi.data());
    writer.EndStep();
}

void output_finalize()
{
    writer.Close();
}
