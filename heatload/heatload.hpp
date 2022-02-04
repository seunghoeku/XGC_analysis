#ifndef HEATLOAD_HPP
#define HEATLOAD_HPP

#include <vector>
#include "particles.hpp"
#include "flags.hpp"

#define N_COND 2
#define N_PSI 1000

//class for conditions
class Conditions{
public:
    bool b[N_COND];
    Conditions();
    Conditions(struct Particles ptl);
};

inline Conditions::Conditions(){
    // empty constructor
}

inline Conditions::Conditions(struct Particles ptl){

    // get conditions
    Flags fl(ptl.flag);
    // get poloidal angle

    // get energy (normalized with T0)


    // check every ifs
    b[0]=true;  // always satisfied.
    b[1]=fl.outboard; // escaped from outboard

}

//class for heatload
class HeatLoad1{
public:
    double en[N_PSI];
    double ptl[N_PSI];    
};

class HeatLoad{
public:
    HeatLoad1 flux[N_COND];
};

#endif
