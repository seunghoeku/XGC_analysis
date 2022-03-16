#include "adios2.h"
#include "heatload_calc.hpp"
#include "sml.hpp"

void output(adios2::ADIOS *ad, HeatLoad &ion, HeatLoad &elec, MPI_Comm comm); // output graphs or data for graphs
void output_finalize();
