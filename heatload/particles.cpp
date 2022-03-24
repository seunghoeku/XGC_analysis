#include "particles.hpp"
#include "adios2.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

#define NPHASE 11
#define GET(X, i, j) X[i * NPHASE + j]

// esc_step is 1-based index. zero means no data in the DB
Particle search(t_ParticleDB &db, int timestep, long long gid)
{
    if (db.count(timestep))
    {
        return db[timestep][gid];
    }
    else
    {
        struct Particle none;
        none.gid = 0;
        return none;
    }
    // return db[timestep][gid];
    // if (timestep > db.size())
    // {
    //     LOG << "[ERROR] timestep is larger than db.size(): " << timestep << " " << db.size();
    //     struct Particle none;
    //     none.gid = 0;
    //     return none;
    // }
    // else if (timestep - 1 < 0)
    // {
    //     struct Particle none;
    //     none.gid = 0;
    //     return none;
    // }
    // else
    // {
    //     // (2022/03/17) jyc: "find" is slow
    //     // unorder_map will return gid=0 particle when no match
    //     return db[timestep - 1][gid];
    // }
}

void insert_or_append(t_ParticleDB &db, int timestep, t_ParticlesList &ptls)
{
    if (db.count(timestep))
    {
        db[timestep].insert(ptls.begin(), ptls.end());
    }
    else
    {
        db.insert({timestep, ptls});
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

template <typename T>
inline void write1d(adios2::IO &io, adios2::Engine &writer, std::string varname, std::vector<T> &val,
                    std::vector<unsigned long int> n)
{
    auto var = io.InquireVariable<T>(varname);
    var.SetShape({n[0]});
    var.SetSelection({{n[1]}, {n[2]}});
    writer.Put<T>(varname, val.data(), adios2::Mode::Sync);
}

template <typename T>
inline void write2d(adios2::IO &io, adios2::Engine &writer, std::string varname, std::vector<T> &val,
                    std::vector<unsigned long int> n)
{
    auto var = io.InquireVariable<T>(varname);
    var.SetShape({n[0], NPHASE});
    var.SetSelection({{n[1], 0}, {n[2], NPHASE}});
    writer.Put<T>(var, val.data(), adios2::Mode::Sync);
}

template <typename T>
inline void read1d(adios2::IO &io, adios2::Engine &reader, std::string varname, std::vector<T> &val,
                   std::vector<unsigned long int> n)
{
    auto var = io.InquireVariable<T>(varname);
    var.SetSelection({{n[0]}, {n[1]}});
    reader.Get<T>(var, val, adios2::Mode::Sync);
}

template <typename T>
inline void read2d(adios2::IO &io, adios2::Engine &reader, std::string varname, std::vector<T> &val,
                   std::vector<unsigned long int> n)
{
    auto var = io.InquireVariable<T>(varname);
    var.SetSelection({{n[0], 0}, {n[1], NPHASE}});
    reader.Get<T>(var, val, adios2::Mode::Sync);
}

void ptldb_save(t_ParticleDB &db, std::string filename, MPI_Comm comm)
{
    int rank, comm_size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    adios2::ADIOS ad;
    adios2::IO io;
    adios2::Engine writer;

    io = ad.DeclareIO("heatload_restart");
    io.DefineVariable<long>("gid", {0}, {0}, {0});
    io.DefineVariable<int>("flag", {0}, {0}, {0});
    io.DefineVariable<int>("step", {0}, {0}, {0});
    io.DefineVariable<float>("dw", {0}, {0}, {0});
    io.DefineVariable<float>("phase", {0, NPHASE}, {0, 0}, {0, NPHASE});
    io.DefineVariable<int>("timestep");

    LOG << "Writing: " << filename << " nsteps: " << db.size();
    writer = io.Open(filename, adios2::Mode::Write, comm);

    for (auto const &pair : db)
    {
        int timestep = pair.first;
        std::vector<long> gid;
        std::vector<int> flag;
        std::vector<int> step;
        std::vector<float> dw;
        std::vector<float> phase;

        LOG << filename << " step,nptl: " << pair.first << " " << pair.second.size();

        for (auto const &pair2 : pair.second)
        {

            const Particle &p = pair2.second;
            gid.push_back(p.gid);
            flag.push_back(p.flag);
            step.push_back(p.esc_step);
            dw.push_back(p.dw);
            phase.push_back(p.r);
            phase.push_back(p.z);
            phase.push_back(p.phi);
            phase.push_back(p.rho);
            phase.push_back(p.w1);
            phase.push_back(p.w2);
            phase.push_back(p.mu);
            phase.push_back(p.w0);
            phase.push_back(p.f0);
            phase.push_back(p.psi);
            phase.push_back(p.B);
        }

        long unsigned int nptl = pair.second.size();
        std::vector<long unsigned int> nptl_list(comm_size);
        std::vector<long unsigned int> displacement_list(comm_size);
        MPI_Allgather(&nptl, 1, MPI_LONG, nptl_list.data(), 1, MPI_LONG, comm);

        long unsigned int ntotal = 0;
        for (int i = 0; i < nptl_list.size(); i++)
        {
            displacement_list[i] = ntotal;
            ntotal += nptl_list[i];
        }

        writer.BeginStep();
        write1d(io, writer, "gid", gid, {ntotal, displacement_list[rank], nptl_list[rank]});
        write1d(io, writer, "flag", flag, {ntotal, displacement_list[rank], nptl_list[rank]});
        write1d(io, writer, "step", step, {ntotal, displacement_list[rank], nptl_list[rank]});
        write1d(io, writer, "dw", dw, {ntotal, displacement_list[rank], nptl_list[rank]});
        write2d(io, writer, "phase", phase, {ntotal, displacement_list[rank], nptl_list[rank]});
        writer.Put<>("timestep", timestep);
        writer.EndStep();
    }
    writer.Close();
}

void ptldb_load(t_ParticleDB &db, std::string filename, MPI_Comm comm)
{
    assert(db.empty());

    int rank, comm_size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    adios2::ADIOS ad;
    adios2::IO io;
    adios2::Engine reader;

    io = ad.DeclareIO("heatload_restart");
    LOG << "Loading: " << filename;
    reader = io.Open(filename, adios2::Mode::Read, comm);
    while (true)
    {
        adios2::StepStatus status = reader.BeginStep();
        if (status == adios2::StepStatus::OK)
        {
            int timestep;
            std::vector<long> gid;
            std::vector<int> flag;
            std::vector<int> step;
            std::vector<float> dw;
            std::vector<float> phase;

            auto var = io.InquireVariable<long>("gid");
            unsigned long int nptl_total = var.Shape()[0];
            // LOG << "nptl_total: " << nptl_total;

            long unsigned int nptl = nptl_total / comm_size;
            long unsigned int offset = nptl * rank;
            if (rank == comm_size - 1)
                nptl = nptl_total - offset;
            // LOG << "offset,nptl= " << offset << " " << nptl;

            // auto var = io.InquireVariable<long>("gid");
            // var.SetSelection({{offset}, {nptl}});
            // reader.Get<long>(var, gid, adios2::Mode::Sync);
            read1d(io, reader, "gid", gid, {offset, nptl});
            read1d(io, reader, "flag", flag, {offset, nptl});
            read1d(io, reader, "step", step, {offset, nptl});
            read1d(io, reader, "dw", dw, {offset, nptl});
            read2d(io, reader, "phase", phase, {offset, nptl});
            reader.Get<int>("timestep", &timestep, adios2::Mode::Sync);
            reader.EndStep();

            t_ParticlesList ptls;
            for (int i = 0; i < gid.size(); i++)
            {
                struct Particle p;
                p.gid = gid[i];
                p.flag = flag[i];
                p.esc_step = step[i];
                p.r = GET(phase, i, 0);
                p.z = GET(phase, i, 1);
                p.phi = GET(phase, i, 2);
                p.rho = GET(phase, i, 3);
                p.w1 = GET(phase, i, 4);
                p.w2 = GET(phase, i, 5);
                p.mu = GET(phase, i, 6);
                p.w0 = GET(phase, i, 7);
                p.f0 = GET(phase, i, 8);
                p.psi = GET(phase, i, 9);
                p.B = GET(phase, i, 10);
                p.dw = dw[i];

                ptls.insert({p.gid, p});
            }
            // LOG << "timestep,ptls.size= " << timestep << " " << ptls.size();

            db.insert({timestep, ptls});
        }
        else
        {
            break;
        }
    }

    reader.Close();
    LOG << "Loading done";
}

void ptldb_print(t_ParticleDB &db, std::string str)
{
    for (auto const &pair : db)
    {
        LOG << str << " info: step " << pair.first << " : " << pair.second.size();
        // for (auto const &pair2 : pair.second)
        // {
        //     LOG << pair2.first << " " << pair2.second.gid;
        //     break;
        // }
    }
}

void ptldb_dump(t_ParticleDB &db, std::string str)
{
    // int istep = 0;
    // for (auto const &pmap : db)
    // {
    //     for (auto const &pair : pmap)
    //     {
    //         LOG << str << " info: step " << istep << " : " << pair.first;
    //     }
    //     istep++;
    // }
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

void ptls_shift(Particles &ptls, Particles &ptls_from_right, MPI_Comm comm)
{
    assert(ptls_from_right.empty());

    int rank, comm_size;
    MPI_Status status;
    int tag_send = 0, tag_recv = 0;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    int left = rank - 1 < 0 ? comm_size - 1 : rank - 1;
    int right = rank + 1 >= comm_size ? 0 : rank + 1;

    int len = ptls.size();
    LOG << "current len: " << len;

    // send to the left and receive from the right
    MPI_Sendrecv_replace(&len, 1, MPI_INT, left, tag_send, right, tag_recv, comm, &status);
    LOG << "new len: " << len;

    ptls_from_right.resize(len);

    int ptls_nbytes = ptls.size() * sizeof(struct Particle);
    int ptls_from_right_nbytes = ptls_from_right.size() * sizeof(struct Particle);

    MPI_Sendrecv(ptls.data(), ptls_nbytes, MPI_CHAR, left, tag_send, ptls_from_right.data(), ptls_from_right_nbytes,
                 MPI_CHAR, right, tag_recv, comm, &status);
}
