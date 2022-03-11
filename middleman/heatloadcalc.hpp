#include <iostream>
#include <vector>
// #include <chrono>
// #include <thread>

// #include "adios2.h"
// #include "flags.hpp"
// #include "load.hpp"
// #include "particles.hpp"
// #include "sml.hpp"
// #include "particles.hpp"
// #include "heatload_calc.hpp"
// #include "output.hpp"

#include "heatload.hpp"

class Heatload
{
  public:
    Heatload(adios2::ADIOS *ad, std::string xgcdir, MPI_Comm comm);

    void finalize();

    adios2::StepStatus step();

  public:
    adios2::ADIOS *ad;
    adios2::IO io;
    adios2::Engine reader;

    MPI_Comm comm;
    int comm_size;
    int rank;

    std::string xgcdir;
    int istep;
    Simulation sml; // input parameters that controls simulation.
    t_ParticleDB iesc_db;
    t_ParticleDB eesc_db;
};
