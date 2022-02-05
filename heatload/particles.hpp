#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <assert.h>
#include <map>
#include <vector>

// particle data strcuture
struct Particle
{
    double r;
    double z;
    double phi;
    double rho;
    double w1;
    double w2;
    double mu;
    double w0;
    double f0;
    double B;
    double psi;
    float dw;
    long long gid;
    int flag;
    int esc_step;
};

typedef std::vector<struct Particle> Particles;
typedef std::map<long long, Particle> t_ParticlesList;
typedef std::vector<t_ParticlesList> t_ParticleDB;

Particle search(t_ParticleDB &db, int timestep, long long int gid);

//#include "particles.tpp"
#endif