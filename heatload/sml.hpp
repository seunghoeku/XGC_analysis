#ifndef SML_HPP
#define SML_HPP

// Simulation parameters
class Simulation {
    public:

    // "Null constructor"

    // Input parameters
    double psix;
    int npsi, ncond;
    
    double rmin[2], rmax[2], zmin[2], zmax[2];
    double pmin[2];
    double dpsi[2];

    double c2_2m;
    
    Simulation();
};

inline Simulation::Simulation() {

};

#endif