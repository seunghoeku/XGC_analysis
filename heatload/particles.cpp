#include "particles.hpp"
#include "adios2.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

// esc_step is 1-based index. zero means no data in the DB
Particle search(t_ParticleDB &db, int timestep, long long gid)
{
    if (timestep > db.size())
    {
        LOG << "[ERROR] timestep is larger than db.size(): " << timestep << " " << db.size();
        struct Particle none;
        none.gid = 0;
        return none;
    }
    else if (timestep - 1 < 0)
    {
        struct Particle none;
        none.gid = 0;
        return none;
    }
    else
    {
        // (2022/03/17) jyc: "find" is slow
        // unorder_map will return gid=0 particle when no match
        return db[timestep - 1][gid];
    }
}

void add(t_ParticlesList &pmap, Particle ptl)
{
    pmap.insert({ptl.gid, ptl});
}

// helper: map-to-vector function
const Particles maptovec(const t_ParticlesList &pmap)
{
    Particles ptls;
    for (auto const &pair : pmap)
    {
        ptls.push_back(pair.second);
    }
    return ptls;
}

void ptldb_save(t_ParticleDB &db, std::string filename)
{
    int istep = 0;

    // vector of particle vector
    std::vector<Particles> ParticleList;

    for (auto const &pmap : db)
    {
        ParticleList.push_back(maptovec(pmap));
    }

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

    LOG << "Writing: " << filename << " nsteps: " << ParticleList.size();
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
    LOG << "Loading done";
}

void ptldb_print(t_ParticleDB &db, std::string str)
{
    int istep = 0;
    for (auto const &pmap : db)
    {
        int nptls = 0;
        nptls = pmap.size();
        LOG << str << " info: step " << istep << " : " << nptls;
        istep++;
    }
}

void ptldb_dump(t_ParticleDB &db, std::string str)
{
    int istep = 0;
    for (auto const &pmap : db)
    {
        for (auto const &pair : pmap)
        {
            LOG << str << " info: step " << istep << " : " << pair.first;
        }
        istep++;
    }
}

int ptlmap_count(t_ParticlesList &pmap)
{
    return pmap.size();
}

void ptlmap_print(t_ParticlesList &pmap, std::string str)
{
    for (auto const &pair : pmap)
    {
        LOG << str << " : " << pair.first << " " << pair.second.gid;
    }
}

void ptlmap_sync(t_ParticlesList &pmap, MPI_Comm comm)
{
    int rank;

    MPI_Comm_rank(comm, &rank);

    Particles ptls;

    if (rank == 0)
    {
        for (auto const &pair : pmap)
        {
            ptls.push_back(pair.second);
        }
    }

    int nptls = ptls.size();
    MPI_Bcast(&nptls, 1, MPI_INT, 0, comm);

    ptls.resize(nptls);
    int nbytes = ptls.size() * sizeof(struct Particle);
    LOG << "ptlmap_sync nbytes: " << nbytes;

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
