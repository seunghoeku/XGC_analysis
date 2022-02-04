#include <vector>

#include "particles.hpp"
#include "flags.hpp"
#include "sml.hpp"
#include "heatload.hpp"

extern Simulation sml;

// get heatload of single species
void heatload_calc(const t_ParticlesList &div, const HeatLoad &flux) {
    

    //inner or outer
    for(int side=0; side <2; side++) {
        //check bounding box


    }


}