#include <assert.h>

#include <string>

#include "flags.hpp"
#include "load.hpp"
#include "sml.hpp"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

#include "cam_timers.hpp"

#define NPHASE 11
#define GET(X, i, j) X[i * NPHASE + j]

adios2::Engine reader;
adios2::IO reader_io;
MPI_Comm reader_comm;
int reader_comm_size;
int reader_comm_rank;

adios2::Engine dup_writer;
adios2::IO dup_io;

template <typename T> inline std::pair<int, int> split_vector(std::vector<T> &vec, int comm_size, int rank)
{
    int nblock = vec.size() / comm_size;
    int offset = nblock * rank;
    if (rank == comm_size - 1)
        nblock = vec.size() - offset;

    return std::make_pair(offset, nblock);
}

template <typename T>
inline void copy_write(adios2::IO &io, adios2::Engine &writer, adios2::Variable<T> &var, std::vector<T> &val)
{
    auto _var = io.InquireVariable<T>(var.Name());
    _var.SetSelection({var.Start(), var.Count()});
    _var.SetShape(var.Shape());
    writer.Put<T>(var.Name(), val.data(), adios2::Mode::Sync);
}

void load_init(adios2::ADIOS *ad, const std::string &filename, MPI_Comm comm)
{
    reader_io = ad->DeclareIO("escaped_ptls"); // same IO name as in XGC
    LOG << "Loading: " << filename;
    reader = reader_io.Open(filename, adios2::Mode::Read, comm);
    reader_comm = comm;
    MPI_Comm_rank(reader_comm, &reader_comm_rank);
    MPI_Comm_size(reader_comm, &reader_comm_size);

    dup_io = ad->DeclareIO("escaped_ptls_dup");
    dup_io.DefineVariable<long>("igid", {0}, {0}, {0});
    dup_io.DefineVariable<long>("egid", {0}, {0}, {0});
    dup_io.DefineVariable<int>("iflag", {0}, {0}, {0});
    dup_io.DefineVariable<int>("eflag", {0}, {0}, {0});
    dup_io.DefineVariable<int>("istep", {0}, {0}, {0});
    dup_io.DefineVariable<int>("estep", {0}, {0}, {0});
    dup_io.DefineVariable<float>("idw", {0}, {0}, {0});
    dup_io.DefineVariable<float>("edw", {0}, {0}, {0});
    dup_io.DefineVariable<float>("ephase", {0, NPHASE}, {0, 0}, {0, NPHASE});
    dup_io.DefineVariable<float>("iphase", {0, NPHASE}, {0, 0}, {0, NPHASE});
    dup_io.DefineVariable<int>("timestep");

    std::string fname = boost::str(boost::format("%s.copy") % filename);
    boost::filesystem::path p(fname);
    dup_writer = dup_io.Open(p.filename().string(), adios2::Mode::Write, reader_comm);
}

void load_finalize()
{
    reader.Close();
    dup_writer.Close();
}

// idiv, ediv (local), iesc, eesc (local)
adios2::StepStatus load_data(Particles &idiv, Particles &ediv, t_ParticlesList &iesc, t_ParticlesList &eesc,
                             int &timestep)
{
    TIMER_START("LOAD_DATA");
    assert(idiv.empty());
    assert(ediv.empty());
    assert(iesc.empty());
    assert(eesc.empty());

    std::vector<long> igid;
    std::vector<long> egid;
    std::vector<int> iflag;
    std::vector<int> eflag;
    std::vector<int> istep;
    std::vector<int> estep;
    std::vector<float> idw;
    std::vector<float> edw;
    std::vector<float> iphase;
    std::vector<float> ephase;

    TIMER_START("ADIOS_STEP");
    adios2::StepStatus status = reader.BeginStep();
    dup_writer.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
        // Inquire variables
        auto var_igid = reader_io.InquireVariable<long>("igid");
        auto var_egid = reader_io.InquireVariable<long>("egid");
        auto var_iflag = reader_io.InquireVariable<int>("iflag");
        auto var_eflag = reader_io.InquireVariable<int>("eflag");
        auto var_istep = reader_io.InquireVariable<int>("istep");
        auto var_estep = reader_io.InquireVariable<int>("estep");
        auto var_idw = reader_io.InquireVariable<float>("idw");
        auto var_edw = reader_io.InquireVariable<float>("edw");
        auto var_iphase = reader_io.InquireVariable<float>("iphase");
        auto var_ephase = reader_io.InquireVariable<float>("ephase");

        var_igid.SetSelection({{0}, {var_igid.Shape()[0]}});
        var_egid.SetSelection({{0}, {var_egid.Shape()[0]}});
        var_iflag.SetSelection({{0}, {var_iflag.Shape()[0]}});
        var_eflag.SetSelection({{0}, {var_eflag.Shape()[0]}});
        var_istep.SetSelection({{0}, {var_istep.Shape()[0]}});
        var_estep.SetSelection({{0}, {var_estep.Shape()[0]}});
        var_idw.SetSelection({{0}, {var_idw.Shape()[0]}});
        var_edw.SetSelection({{0}, {var_edw.Shape()[0]}});
        var_iphase.SetSelection({{0, 0}, {var_iphase.Shape()[0], var_iphase.Shape()[1]}});
        var_ephase.SetSelection({{0, 0}, {var_ephase.Shape()[0], var_ephase.Shape()[1]}});

        int adios_step = reader.CurrentStep();
        auto block_list_igid = reader.BlocksInfo(var_igid, adios_step);
        auto slice = split_vector(block_list_igid, reader_comm_size, reader_comm_rank);
        int offset = slice.first;
        int nblock = slice.second;
        printf("%d: offset,nblock= %d %d\n", reader_comm_rank, offset, nblock);

        // Read table block by block
        for (int i = offset; i < offset + nblock; i++)
        {
            std::vector<long> _igid;
            std::vector<long> _egid;
            std::vector<int> _iflag;
            std::vector<int> _eflag;
            std::vector<int> _istep;
            std::vector<int> _estep;
            std::vector<float> _idw;
            std::vector<float> _edw;
            std::vector<float> _iphase;
            std::vector<float> _ephase;

            auto block = block_list_igid[i];
            int ncount = 1;
            for (auto &d : block.Count)
            {
                ncount *= d;
            }

            if (ncount > 0)
            {
                var_igid.SetBlockSelection(block.BlockID);
                var_egid.SetBlockSelection(block.BlockID);
                var_iflag.SetBlockSelection(block.BlockID);
                var_eflag.SetBlockSelection(block.BlockID);
                var_istep.SetBlockSelection(block.BlockID);
                var_estep.SetBlockSelection(block.BlockID);
                var_idw.SetBlockSelection(block.BlockID);
                var_edw.SetBlockSelection(block.BlockID);
                var_iphase.SetBlockSelection(block.BlockID);
                var_ephase.SetBlockSelection(block.BlockID);

                reader.Get<long>(var_igid, _igid);
                reader.Get<long>(var_egid, _egid);
                reader.Get<int>(var_iflag, _iflag);
                reader.Get<int>(var_eflag, _eflag);
                reader.Get<int>(var_istep, _istep);
                reader.Get<int>(var_estep, _estep);
                reader.Get<float>(var_idw, _idw);
                reader.Get<float>(var_edw, _edw);
                reader.Get<float>(var_iphase, _iphase);
                reader.Get<float>(var_ephase, _ephase);
                reader.Get<int>("timestep", &timestep);
                TIMER_START("ADIOS_PERFORM_GETS");
                reader.PerformGets();
                TIMER_STOP("ADIOS_PERFORM_GETS");

                TIMER_START("_ADIOS_DUP_WRITE");
                copy_write(dup_io, dup_writer, var_igid, _igid);
                copy_write(dup_io, dup_writer, var_egid, _egid);
                copy_write(dup_io, dup_writer, var_iflag, _iflag);
                copy_write(dup_io, dup_writer, var_eflag, _eflag);
                copy_write(dup_io, dup_writer, var_istep, _istep);
                copy_write(dup_io, dup_writer, var_estep, _estep);
                copy_write(dup_io, dup_writer, var_idw, _idw);
                copy_write(dup_io, dup_writer, var_edw, _edw);
                copy_write(dup_io, dup_writer, var_iphase, _iphase);
                copy_write(dup_io, dup_writer, var_ephase, _ephase);
                dup_writer.Put<int>("timestep", timestep);
                TIMER_STOP("_ADIOS_DUP_WRITE");

                igid.insert(igid.end(), _igid.begin(), _igid.end());
                egid.insert(egid.end(), _egid.begin(), _egid.end());
                iflag.insert(iflag.end(), _iflag.begin(), _iflag.end());
                eflag.insert(eflag.end(), _eflag.begin(), _eflag.end());
                istep.insert(istep.end(), _istep.begin(), _istep.end());
                estep.insert(estep.end(), _estep.begin(), _estep.end());
                idw.insert(idw.end(), _idw.begin(), _idw.end());
                edw.insert(edw.end(), _edw.begin(), _edw.end());
                iphase.insert(iphase.end(), _iphase.begin(), _iphase.end());
                ephase.insert(ephase.end(), _ephase.begin(), _ephase.end());
            }
        }
        reader.EndStep();
        dup_writer.EndStep();
    }
    TIMER_STOP("ADIOS_STEP");

    if (status == adios2::StepStatus::OK)
    {
        assert(iphase.size() / igid.size() == NPHASE);
        assert(ephase.size() / egid.size() == NPHASE);

        // Each will sort out esc or div particles and build own DB
#pragma omp parallel for default(none) shared(igid, iflag, istep, iphase, idw, iesc, idiv, std::cerr)
        for (int i = 0; i < igid.size(); i++)
        {
            struct Particle iptl;
            iptl.gid = igid[i];
            iptl.flag = iflag[i];
            iptl.esc_step = istep[i];
            iptl.r = GET(iphase, i, 0);
            iptl.z = GET(iphase, i, 1);
            iptl.phi = GET(iphase, i, 2);
            iptl.rho = GET(iphase, i, 3);
            iptl.w1 = GET(iphase, i, 4);
            iptl.w2 = GET(iphase, i, 5);
            iptl.mu = GET(iphase, i, 6);
            iptl.w0 = GET(iphase, i, 7);
            iptl.f0 = GET(iphase, i, 8);
            iptl.psi = GET(iphase, i, 9);
            iptl.B = GET(iphase, i, 10);
            iptl.dw = idw[i];

            int flag1; // tmp flag
            flag1 = iflag[i];

            Flags fl(flag1); // decode flags

            if (fl.escaped)
            {
#pragma omp critical
                {
                    // add to esc
                    add(iesc, iptl);
                }
            }
            else
            {
#pragma omp critical
                {
                    // add to div
                    idiv.push_back(iptl);
                }
            }
        }
        LOG << "Done with iesc and idiv";

        // populate ediv with local data
#pragma omp parallel for default(none) shared(egid, eflag, estep, ephase, edw, eesc, ediv, std::cerr)
        for (int i = 0; i < egid.size(); i++)
        {
            struct Particle eptl;
            eptl.gid = egid[i];
            eptl.flag = eflag[i];
            eptl.esc_step = estep[i];
            eptl.r = GET(ephase, i, 0);
            eptl.z = GET(ephase, i, 1);
            eptl.phi = GET(ephase, i, 2);
            eptl.rho = GET(ephase, i, 3);
            eptl.w1 = GET(ephase, i, 4);
            eptl.w2 = GET(ephase, i, 5);
            eptl.mu = GET(ephase, i, 6);
            eptl.w0 = GET(ephase, i, 7);
            eptl.f0 = GET(ephase, i, 8);
            eptl.psi = GET(ephase, i, 9);
            eptl.B = GET(ephase, i, 10);
            eptl.dw = edw[i];

            int flag1; // tmp flag
            flag1 = eflag[i];

            Flags fl(flag1); // decode flags

            // save to div or esc
            if (fl.escaped)
            {
#pragma omp critical
                {
                    // add to esc
                    add(eesc, eptl);
                }
            }
            else
            {
#pragma omp critical
                {
                    // add to div
                    ediv.push_back(eptl);
                }
            }
        }
        LOG << "Done with eesc and ediv";
    }

    TIMER_STOP("LOAD_DATA");
    return status;
}
