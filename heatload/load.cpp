#include <assert.h>

#include <string>

#include "flags.hpp"
#include "load.hpp"
#include "sml.hpp"

#define NPHASE 11
#define GET(X, i, j) X[i * NPHASE + j]

adios2::Engine reader;
adios2::IO reader_io;
MPI_Comm reader_comm;
int reader_comm_size;
int reader_comm_rank;

template <typename T> inline std::pair<int, int> split_vector(std::vector<T> &vec, int comm_size, int rank)
{
    int nblock = vec.size() / comm_size;
    int offset = nblock * rank;
    if (rank == comm_size - 1)
        nblock = vec.size() - offset;

    return std::make_pair(offset, nblock);
}

void load_init(adios2::ADIOS *ad, const std::string &filename, MPI_Comm comm)
{
    reader_io = ad->DeclareIO("escaped_ptls"); // same IO name as in XGC
    reader = reader_io.Open(filename, adios2::Mode::Read, comm);
    reader_comm = comm;
    MPI_Comm_rank(reader_comm, &reader_comm_rank);
    MPI_Comm_size(reader_comm, &reader_comm_size);
}

void load_finalize()
{
    reader.Close();
}

adios2::StepStatus load_data(Particles &idiv, Particles &ediv, t_ParticlesList &iesc, t_ParticlesList &eesc)
{
    // Clear vector
    idiv.clear();
    ediv.clear();
    iesc.clear();
    eesc.clear();

    adios2::StepStatus status = reader.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
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
            auto block = block_list_igid[i];
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
            reader.PerformGets();

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

        assert(iphase.size() / igid.size() == NPHASE);
        assert(ephase.size() / egid.size() == NPHASE);

        // Merge to rank 0
        int len = igid.size();
        std::vector<int> len_list(reader_comm_size);
        std::vector<int> displacement_list(reader_comm_size);

        MPI_Allgather(&len, 1, MPI_INT, len_list.data(), 1, MPI_INT, reader_comm);

        int ntotal = 0;
        for (int i = 0; i < len_list.size(); i++)
        {
            displacement_list[i] = ntotal;
            ntotal += len_list[i];
        }

        std::vector<long> igid_total(ntotal);
        std::vector<int> iflag_total(ntotal);
        std::vector<int> istep_total(ntotal);
        std::vector<float> idw_total(ntotal);
        std::vector<float> iphase_total(ntotal * NPHASE);

        MPI_Gatherv(igid.data(), igid.size(), MPI_LONG, igid_total.data(), len_list.data(), displacement_list.data(),
                    MPI_LONG, 0, reader_comm);
        MPI_Gatherv(iflag.data(), iflag.size(), MPI_INT, iflag_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, reader_comm);
        MPI_Gatherv(istep.data(), istep.size(), MPI_INT, istep_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, reader_comm);
        MPI_Gatherv(idw.data(), idw.size(), MPI_FLOAT, idw_total.data(), len_list.data(), displacement_list.data(),
                    MPI_FLOAT, 0, reader_comm);

        for (int i = 0; i < len_list.size(); i++)
        {
            len_list[i] = len_list[i] * NPHASE;
        }

        ntotal = 0;
        for (int i = 0; i < len_list.size(); i++)
        {
            displacement_list[i] = ntotal;
            ntotal += len_list[i];
        }

        MPI_Gatherv(iphase.data(), iphase.size(), MPI_FLOAT, iphase_total.data(), len_list.data(),
                    displacement_list.data(), MPI_FLOAT, 0, reader_comm);

        // Electron
        len = egid.size();
        MPI_Allgather(&len, 1, MPI_INT, len_list.data(), 1, MPI_INT, reader_comm);

        ntotal = 0;
        for (int i = 0; i < len_list.size(); i++)
        {
            // LOG << boost::format("%d %d") % i % len_list[i];
            displacement_list[i] = ntotal;
            ntotal += len_list[i];
        }

        std::vector<long> egid_total(ntotal);
        std::vector<int> eflag_total(ntotal);
        std::vector<int> estep_total(ntotal);
        std::vector<float> edw_total(ntotal);
        std::vector<float> ephase_total(ntotal * NPHASE);

        MPI_Gatherv(egid.data(), egid.size(), MPI_LONG, egid_total.data(), len_list.data(), displacement_list.data(),
                    MPI_LONG, 0, reader_comm);
        MPI_Gatherv(eflag.data(), eflag.size(), MPI_INT, eflag_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, reader_comm);
        MPI_Gatherv(estep.data(), estep.size(), MPI_INT, estep_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, reader_comm);
        MPI_Gatherv(edw.data(), edw.size(), MPI_FLOAT, edw_total.data(), len_list.data(), displacement_list.data(),
                    MPI_FLOAT, 0, reader_comm);

        for (int i = 0; i < len_list.size(); i++)
        {
            len_list[i] = len_list[i] * NPHASE;
        }

        ntotal = 0;
        for (int i = 0; i < len_list.size(); i++)
        {
            displacement_list[i] = ntotal;
            ntotal += len_list[i];
        }

        MPI_Gatherv(ephase.data(), ephase.size(), MPI_FLOAT, ephase_total.data(), len_list.data(),
                    displacement_list.data(), MPI_FLOAT, 0, reader_comm);

        if (reader_comm_rank == 0)
        {
            // populate particles
            for (int i = 0; i < igid_total.size(); i++)
            {
                struct Particle iptl;
                iptl.gid = igid_total[i];
                iptl.flag = iflag_total[i];
                iptl.esc_step = istep_total[i];
                iptl.r = GET(iphase_total, i, 0);
                iptl.z = GET(iphase_total, i, 1);
                iptl.phi = GET(iphase_total, i, 2);
                iptl.rho = GET(iphase_total, i, 3);
                iptl.w1 = GET(iphase_total, i, 4);
                iptl.w2 = GET(iphase_total, i, 5);
                iptl.mu = GET(iphase_total, i, 6);
                iptl.w0 = GET(iphase_total, i, 7);
                iptl.f0 = GET(iphase_total, i, 8);
                iptl.psi = GET(iphase_total, i, 9);
                iptl.B = GET(iphase_total, i, 10);
                iptl.dw = idw_total[i];

                int flag1; // tmp flag
                flag1 = iflag_total[i];

                Flags fl(flag1); // decode flags

                // save to div or esc
                if (fl.escaped)
                {
                    // add to esc
                    iesc.insert(std::pair<long long, Particle>(iptl.gid, iptl));
                }
                else
                {
                    // add to div
                    idiv.push_back(iptl);
                }
            }

            for (int i = 0; i < egid_total.size(); i++)
            {
                struct Particle eptl;
                eptl.gid = egid_total[i];
                eptl.flag = eflag_total[i];
                eptl.esc_step = estep_total[i];
                eptl.r = GET(ephase_total, i, 0);
                eptl.z = GET(ephase_total, i, 1);
                eptl.phi = GET(ephase_total, i, 2);
                eptl.rho = GET(ephase_total, i, 3);
                eptl.w1 = GET(ephase_total, i, 4);
                eptl.w2 = GET(ephase_total, i, 5);
                eptl.mu = GET(ephase_total, i, 6);
                eptl.w0 = GET(ephase_total, i, 7);
                eptl.f0 = GET(ephase_total, i, 8);
                eptl.psi = GET(ephase_total, i, 9);
                eptl.B = GET(ephase_total, i, 10);
                eptl.dw = edw_total[i];

                int flag1; // tmp flag
                flag1 = eflag_total[i];

                Flags fl(flag1); // decode flags

                // save to div or esc
                if (fl.escaped)
                {
                    // add to esc
                    eesc.insert(std::pair<long long, Particle>(eptl.gid, eptl));
                }
                else
                {
                    // add to div
                    ediv.push_back(eptl);
                }
            }
        }
        reader.EndStep();
    }

    return status;
}
