#include "particles.hpp"

Particle search(t_ParticleDB &db, int timestep, long long gid)
{
    assert(timestep < db.size());

    Particle nullptl;
    nullptl.gid = -1;

    t_ParticlesList ptls = db[timestep];

    if (ptls.find(gid) != ptls.end())
        return ptls[gid];
    else
        return nullptl;
}
