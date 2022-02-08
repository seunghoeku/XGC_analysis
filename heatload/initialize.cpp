#include "sml.hpp"
#include "adios2.h"

extern Simulation sml;

void init() {
    adios2::ADIOS ad;
    adios2::Engine reader;
    adios2::IO reader_io;

    reader_io = ad.DeclareIO("init");
    reader = reader_io.Open("xgc.units.bp", adios2::Mode::Read);

    double eq_x_r = 0.0;
    int ptl_num = 0;

    reader.Get<double>("eq_x_r", &eq_x_r);
    reader.Get<int>("ptl_num", &ptl_num);
    reader.Close();

    printf ("eq_x_r= %g\n", eq_x_r);
    printf ("ptl_num= %d\n", ptl_num);
}