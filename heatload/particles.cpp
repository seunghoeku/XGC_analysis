#include "particles.hpp"

Particle search(t_ParticleDB &db, int timestep, long long gid)
{
    assert(timestep < db.size());

    Particle ptl;
    ptl.gid = -1;

    t_ParticlesList ptls = db[timestep];
    auto it = ptls.find(gid);

    if (it != ptls.end())
        ptl = it->second;

    return ptl;
}
