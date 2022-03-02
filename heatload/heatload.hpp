#include <assert.h>
#include <string>
#include <chrono>
#include <thread>

#include "adios2.h"
#include "flags.hpp"
#include "load.hpp"
#include "particles.hpp"
#include "sml.hpp"
#include "particles.hpp"
#include "heatload_calc.hpp"
#include "output.hpp"

void heatload(adios2::ADIOS * ad);
void heatload_init(adios2::ADIOS *ad);
void heatload_init2(adios2::ADIOS *ad);
int heatload_step(adios2::ADIOS *ad, int istep);
void heatload_finalize(adios2::ADIOS *ad);

