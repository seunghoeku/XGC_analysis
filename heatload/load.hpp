#include "adios2.h"
#include "particles.hpp"

void load_init(const std::string &filename);
void load_finalize();

adios2::StepStatus load_data(Particles &idiv, Particles &ediv, t_ParticlesList &iesc,
                             t_ParticlesList &eesc);