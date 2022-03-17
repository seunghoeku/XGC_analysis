#ifndef HEATLOAD_HPP
#define HEATLOAD_HPP

#include "flags.hpp"
#include "particles.hpp"
#include "sml.hpp"
#include <iostream>
#include <math.h>
#include <vector>

#define N_SIDE 2
#define N_THETA 8
#define N_COND (N_THETA + 1)
#define N_PSI 1000

extern Simulation sml;

// class for conditions
class Conditions
{
  public:
    bool b[N_COND];
    Conditions();
    Conditions(struct Particle ptl);
    inline double get_angle(double x, double y);
};

inline Conditions::Conditions()
{
    // empty constructor
}

inline Conditions::Conditions(struct Particle ptl)
{

    // get conditions
    Flags fl(ptl.flag);

    // get poloidal angle
    double theta;
    theta = get_angle(ptl.r - sml.axis_r, ptl.z - sml.axis_z);
    theta = fmod(theta - sml.x_theta + 2 * M_PI, 2. * M_PI);

    // get energy (normalized with T0)

    // check every ifs
    b[0] = true; // always satisfied.
    for (int i = 1; i <= N_THETA; i++)
        b[i] = false; // default false
    int itheta = theta / sml.dtheta;
    itheta = std::min(itheta, N_THETA); // itheta = 0 - N_THETA-1
    b[itheta + 1] = true;               // set true for the angle segment

    //    b[1]=fl.outboard; // escaped from outboard
}

inline double Conditions::get_angle(double x, double y)
{
    double angle = atan2(y, x);
    if (x < 0.)
    {
        angle += M_PI;
    }
    return (angle);
}

// class for heatload
class HeatLoad1
{
  public:
    double en[N_COND][N_PSI];
    double ptl[N_COND][N_PSI];
};

class HeatLoad
{
  public:
    HeatLoad1 side[N_SIDE];
    int isp; // species index

    HeatLoad();
    HeatLoad(int isp_in);
};

inline HeatLoad::HeatLoad()
{
    for (int is = 0; is < N_SIDE; is++)
    {
        for (int ic = 0; ic < N_COND; ic++)
        {
            for (int i = 0; i < N_PSI; i++)
            {
                side[is].en[ic][i] = 0.0;
                side[is].ptl[ic][i] = 0.0;
            }
        }
    }

} // empty constructor

inline HeatLoad::HeatLoad(int isp_in)
{
    isp = isp_in;

    for (int is = 0; is < N_SIDE; is++)
    {
        for (int ic = 0; ic < N_COND; ic++)
        {
            for (int i = 0; i < N_PSI; i++)
            {
                side[is].en[ic][i] = 0.0;
                side[is].ptl[ic][i] = 0.0;
            }
        }
    }
}
#endif
