#ifndef SML_HPP
#define SML_HPP

// Simulation parameters
class Simulation {
    public:

    // "Null constructor"

    // Input parameters
    double psix, x_r, x_z;
    double axis_r, axis_z;
    double x_theta, dtheta; // X-point angle, delta-angle
    double dt; // dt to get heatload unit
    int npsi, ncond, ntheta;
    
    double rmin[2], rmax[2], zmin[2], zmax[2];
    double pmin[2];
    double pmax[2];
    double dpsi[2];

    double c2_2m[2];

    Simulation();
};

inline Simulation::Simulation() {

};

#endif