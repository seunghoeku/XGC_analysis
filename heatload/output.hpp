#include "adios2.h"
#include "sml.hpp"
#include "heatload_calc.hpp"

void output(adios2::ADIOS *ad, HeatLoad &ion, HeatLoad &elec); // output graphs or data for graphs
void output_finalize();
