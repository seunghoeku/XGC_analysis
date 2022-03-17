#ifndef SML_HPP
#define SML_HPP

// Simulation parameters
class Simulation
{
  public:
    // "Null constructor"

    // Input parameters
    double psix, x_r, x_z;
    double axis_r, axis_z;
    double x_theta, dtheta; // X-point angle, delta-angle
    double dt;              // dt to get heatload unit
    int npsi, ncond, ntheta;

    double rmin[2], rmax[2], zmin[2], zmax[2];
    double pmin[2];
    double pmax[2];
    double dpsi[2];

    double c2_2m[2];

    Simulation();
};

inline Simulation::Simulation()
{
    rmin[0] = 0;
    rmin[1] = 0;

    rmax[0] = 0;
    rmax[1] = 0;

    zmin[0] = 0;
    zmin[1] = 0;

    zmax[0] = 0;
    zmax[1] = 0;

    pmin[0] = 0;
    pmin[1] = 0;

    pmax[0] = 0;
    pmax[1] = 0;

    dpsi[0] = 0;
    dpsi[1] = 0;

    c2_2m[0] = 0;
    c2_2m[1] = 0;
};

#endif