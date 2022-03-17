#include <assert.h>
#include <chrono>
#include <string>
#include <thread>

#include "adios2.h"
#include "flags.hpp"
#include "heatload_calc.hpp"
#include "load.hpp"
#include "output.hpp"
#include "particles.hpp"
#include "sml.hpp"

void heatload(adios2::ADIOS *ad);
void heatload_init(adios2::ADIOS *ad, MPI_Comm comm, std::string xgcdir, bool read_restart);
void heatload_init2(adios2::ADIOS *ad, std::string xgcdir);
int heatload_step(adios2::ADIOS *ad, int istep, bool ion_only);
void heatload_finalize();
