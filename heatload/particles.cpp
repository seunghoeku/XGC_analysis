#include "particles.hpp"
#include "adios2.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

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

void ptldb_save(t_ParticleDB &db, std::string filename)
{
    std::vector<Particles> ParticleList;

    int istep = 0;
    for (auto const &plist : db)
    {
        Particles ptls;
        for (auto const &inner : plist)
        {
            for (auto const &ptl : inner.second)
            {
                ptls.push_back(ptl.second);
            }
        }
        ParticleList.push_back(ptls);
    }
    // LOG << "PTLS: " << ParticleList.size();

    adios2::ADIOS ad;
    adios2::IO io;
    adios2::Engine writer;

    io = ad.DeclareIO("heatload_restart");
    io.DefineVariable<int>("nstep");
    istep = 0;
    for (auto const &ptls : ParticleList)
    {
        // LOG << ptls.size() << " " << sizeof(struct Particle) << " " << ptls.size() * sizeof(struct Particle);
        std::string vname = boost::str(boost::format("ptls%d") % istep);
        int nbytes = ptls.size() * sizeof(struct Particle);
        io.DefineVariable<char>(vname, {(long unsigned int)nbytes}, {0}, {(long unsigned int)nbytes});
        istep++;
    }

    LOG << "Writing: " << filename;
    writer = io.Open(filename, adios2::Mode::Write, MPI_COMM_SELF);
    istep = 0;
    for (auto const &ptls : ParticleList)
    {
        std::string vname = boost::str(boost::format("ptls%d") % istep);
        writer.Put<char>(vname, (char *)ptls.data());
        istep++;
    }
    writer.Put<int>("nstep", ParticleList.size());
    writer.Close();
}

void ptldb_load(t_ParticleDB &db, std::string filename)
{
    adios2::ADIOS ad;
    adios2::IO io;
    adios2::Engine reader;

    int nstep = 0;
    std::vector<Particles> ParticleList;

    if (boost::filesystem::exists(filename))
    {
        io = ad.DeclareIO("heatload_restart");
        LOG << "Loading: " << filename;
        reader = io.Open(filename, adios2::Mode::Read, MPI_COMM_SELF);

        reader.Get<int>("nstep", nstep);
        reader.PerformGets();

        for (int istep = 0; istep < nstep; istep++)
        {
            Particles ptls;
            std::string vname = boost::str(boost::format("ptls%d") % istep);
            auto var = io.InquireVariable<char>(vname);
            ptls.resize(var.Count()[0]);
            reader.Get<char>(vname, (char *)ptls.data());
            reader.PerformGets();
            ParticleList.push_back(ptls);
        }
        reader.Close();
    }
    else
    {
        LOG << "No db file: " << filename;
    }

    for (int istep = 0; istep < nstep; istep++)
    {
        t_ParticlesList pmap;
        Particles ptls = ParticleList[istep];
        for (auto ptl : ptls)
        {
            add(pmap, ptl);
        }
        db.push_back(pmap);
    }
}

void ptldb_print(t_ParticleDB &db, std::string str)
{
    int istep = 0;
    for (auto const &plist : db)
    {
        int nptls = 0;
        for (auto const &inner : plist)
        {
            nptls += inner.second.size();
        }
        LOG << str << " info: step " << istep << " : " << nptls;
        istep++;
    }
}

void ptlmap_sync(t_ParticlesList &pmap, MPI_Comm comm)
{
    int rank;

    MPI_Comm_rank(comm, &rank);

    Particles ptls;

    if (rank == 0)
    {
        for (auto const &inner : pmap)
        {
            for (auto const &ptl : inner.second)
            {
                ptls.push_back(ptl.second);
            }
        }
    }

    int nptls = ptls.size();
    MPI_Bcast(&nptls, 1, MPI_INT, 0, comm);

    ptls.resize(nptls);
    int nbytes = ptls.size() * sizeof(struct Particle);

    MPI_Bcast(ptls.data(), nbytes, MPI_CHAR, 0, comm);

    // Reconstruct pmap by using the data from rank 0
    if (rank > 0)
    {
        for (auto const &ptl : ptls)
        {
            add(pmap, ptl);
        }
    }
}
