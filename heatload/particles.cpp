#include "particles.hpp"

// Particle search(t_ParticleDB &db, int timestep, long long gid)
// {
//     assert(timestep < db.size());

//     Particle ptl;
//     ptl.gid = -1;

//     t_ParticlesList ptls = db[timestep];
//     auto it = ptls.find(gid);

//     if (it != ptls.end())
//         ptl = it->second;

//     return ptl;
// }

Particle search(t_ParticleDB &db, int timestep, long long gid)
{
    assert(timestep < db.size());

    Particle ptl;
    ptl.gid = -1;

    t_ParticlesList pmap = db[timestep];
    auto it = pmap.find(gid / MMOD);

    if (it != pmap.end())
    {
        auto inner = it->second;
        auto init = inner.find(gid);
        if (init != inner.end())
        {
            ptl = init->second;
        }
    }

    return ptl;
}

void add(t_ParticlesList &pmap, Particle ptl)
{
    // pmap.insert(std::pair<long long, Particle>(ptl.gid, ptl));
    auto it = pmap.find(ptl.gid / MMOD);
    if (it != pmap.end())
    {
        auto inner = it->second;
        inner.insert(std::pair<long long, Particle>(ptl.gid, ptl));
    }
    else
    {
        t_ParticlesListInner inner;
        inner.insert(std::pair<long long, Particle>(ptl.gid, ptl));
        pmap.insert(std::pair<long long, t_ParticlesListInner>(ptl.gid / MMOD, inner));
    }
}