#include "adios2.h"
#include "heatload.hpp"

adios2::ADIOS ad2;
adios2::Engine output;
adios2::IO output_io;

void output(HeatLoad &ion, HeatLoad &elec) {
    static bool first = true;

    output_io = ad2.DeclareIO("output");
    
    if(first) {
        output = output_io.Open(filename, adios2::Mode::Write);
        first = false;
    }
    else {
        output = output_io.Open(filename, adios2::Mode::Append);
    }


    double psi[N_SIDE][N_PSI];
    for(int is=0; is<N_SIDE; is++) {
        for(int i=0; i<N_PSI; i++){
           psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
        }
    }

    // save psi, ion.side[0:N_SIDE].en[0:N_COND][0:N_PSI] and ptl[0:NCOND][0:N_PSI] and electron.

}

void output_finalize()
{
    output.Close();
}
