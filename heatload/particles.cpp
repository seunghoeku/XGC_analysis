#include "particles.hpp"

Particles search(t_ParticleDB &db, int timestep, long long int gid)
{
    assert(timestep < db.size());

    Particles nullptl;
    nullptl.gid = -1;

    t_ParticlesList ptls = db[timestep];

    for (int i = 0; i < ptls.size(); i++)
    {
        if (ptls[i].gid == gid)
            return ptls[i];
    }

    return nullptl;
}