#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "mpi.h"

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
    Particle() : gid(0) {}
};

typedef std::vector<struct Particle> Particles;
typedef std::unordered_map<long long, struct Particle> t_ParticlesList;
// typedef std::vector<t_ParticlesList> t_ParticleDB;
typedef std::map<int, t_ParticlesList> t_ParticleDB;

void set_nbin(uint64_t nbin);
Particle& search(t_ParticleDB &db, int timestep, long long int gid);
void add(t_ParticlesList &pmap, Particle ptl);
void ptldb_save(t_ParticleDB &db, std::string filename, MPI_Comm comm);
void ptldb_load(t_ParticleDB &db, std::string filename, MPI_Comm comm);
void ptldb_print(t_ParticleDB &db, std::string str);
int ptlmap_count(t_ParticlesList &pmap);
void ptlmap_sync(t_ParticlesList &pmap, MPI_Comm comm);
void ptlmap_print(t_ParticlesList &pmap, std::string str);
void ptldb_dump(t_ParticleDB &db, std::string str);
void ptls_shift(Particles &ptls, Particles &ptls_from_right, MPI_Comm comm);
void insert_or_append(t_ParticleDB &db, int timestep, t_ParticlesList &ptls);

//#include "particles.tpp"
#endif
