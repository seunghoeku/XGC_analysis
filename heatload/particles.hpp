#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "mpi.h"

#define MMOD 1'000'000
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
// typedef std::unordered_map<long long, Particle> t_ParticlesList;
typedef std::unordered_map<long long, Particle> t_ParticlesListInner;
typedef std::unordered_map<long long, t_ParticlesListInner> t_ParticlesList;
typedef std::vector<t_ParticlesList> t_ParticleDB;

Particle search(t_ParticleDB &db, int timestep, long long int gid);
void add(t_ParticlesList &pmap, Particle ptl);
void ptldb_save(t_ParticleDB &db, std::string filename);
void ptldb_load(t_ParticleDB &db, std::string filename);
void ptldb_print(t_ParticleDB &db, std::string str);
void ptlmap_sync(t_ParticlesList &pmap, MPI_Comm comm);

//#include "particles.tpp"
#endif
