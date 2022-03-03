#include "adios2.h"
#include "particles.hpp"

void load_init(adios2::ADIOS *ad, const std::string &filename, MPI_Comm comm);
void load_finalize();

adios2::StepStatus load_data(Particles &idiv, Particles &ediv, t_ParticlesList &iesc, t_ParticlesList &eesc);