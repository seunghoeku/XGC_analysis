#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <assert.h>
#include <map>
#include <unordered_map>
#include <vector>

// particle data strcuture
// float is used to save memory. Check any underflow happens??
struct Particle
{
    float r;
    float z;
    float phi;
    float rho;
    float w1;
    float w2;
    float mu;
    float w0;
    float f0;
    float B;
    float psi;
    float dw;
    long long gid;
    int flag;
    int esc_step;
};

typedef std::vector<struct Particle> Particles;
// typedef std::map<long long, Particle> t_ParticlesList;
typedef std::unordered_map<long long, Particle> t_ParticlesList;
typedef std::vector<t_ParticlesList> t_ParticleDB;

Particle search(t_ParticleDB &db, int timestep, long long int gid);

//#include "particles.tpp"
#endif