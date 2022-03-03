#include "heatloadcalc.hpp"
#include "adios2.h"
#include "util.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

#define NPHASE 11
#define GET(X, i, j) X[i * NPHASE + j]

void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db); // calculate heatload

Heatload::Heatload(adios2::ADIOS *ad, MPI_Comm comm)
{
    this->ad = ad;

    this->comm = comm;
    MPI_Comm_rank(comm, &this->rank);
    MPI_Comm_size(comm, &this->comm_size);

    this->istep = 0;

    heatload_init2(ad);

    this->io = ad->DeclareIO("escaped_ptls"); // same IO name as in XGC
    this->reader = this->io.Open("xgc.escaped_ptls.bp", adios2::Mode::Read, comm);
}

void Heatload::finalize()
{
    this->reader.Close();
    heatload_finalize(this->ad);
}

adios2::StepStatus Heatload::step()
{
    Particles idiv;
    Particles ediv;
    t_ParticlesList iesc;
    t_ParticlesList eesc;

    adios2::StepStatus status = this->reader.BeginStep();
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
        auto var_igid = this->io.InquireVariable<long>("igid");
        auto var_egid = this->io.InquireVariable<long>("egid");
        auto var_iflag = this->io.InquireVariable<int>("iflag");
        auto var_eflag = this->io.InquireVariable<int>("eflag");
        auto var_istep = this->io.InquireVariable<int>("istep");
        auto var_estep = this->io.InquireVariable<int>("estep");
        auto var_idw = this->io.InquireVariable<float>("idw");
        auto var_edw = this->io.InquireVariable<float>("edw");
        auto var_iphase = this->io.InquireVariable<float>("iphase");
        auto var_ephase = this->io.InquireVariable<float>("ephase");

        auto block_list_igid = this->reader.BlocksInfo(var_igid, this->istep);

        auto slice = split_vector(block_list_igid, this->comm_size, this->rank);
        LOG << boost::format("offset,nblock= %d %d") % slice.first % slice.second;

        int offset = slice.first;
        int nblock = slice.second;

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

            this->reader.Get<long>(var_igid, _igid);
            this->reader.Get<long>(var_egid, _egid);
            this->reader.Get<int>(var_iflag, _iflag);
            this->reader.Get<int>(var_eflag, _eflag);
            this->reader.Get<int>(var_istep, _istep);
            this->reader.Get<int>(var_estep, _estep);
            this->reader.Get<float>(var_idw, _idw);
            this->reader.Get<float>(var_edw, _edw);
            this->reader.Get<float>(var_iphase, _iphase);
            this->reader.Get<float>(var_ephase, _ephase);
            this->reader.PerformGets();

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

        // Merge to rank 0
        int len = igid.size();
        std::vector<int> len_list(this->comm_size);
        std::vector<int> displacement_list(this->comm_size);

        MPI_Allgather(&len, 1, MPI_INT, len_list.data(), 1, MPI_INT, this->comm);

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
                    MPI_LONG, 0, this->comm);
        MPI_Gatherv(iflag.data(), iflag.size(), MPI_INT, iflag_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, this->comm);
        MPI_Gatherv(istep.data(), istep.size(), MPI_INT, istep_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, this->comm);
        MPI_Gatherv(idw.data(), idw.size(), MPI_FLOAT, idw_total.data(), len_list.data(), displacement_list.data(),
                    MPI_FLOAT, 0, this->comm);

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
                    displacement_list.data(), MPI_FLOAT, 0, this->comm);

        // Electron
        len = egid.size();
        MPI_Allgather(&len, 1, MPI_INT, len_list.data(), 1, MPI_INT, this->comm);

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
                    MPI_LONG, 0, this->comm);
        MPI_Gatherv(eflag.data(), eflag.size(), MPI_INT, eflag_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, this->comm);
        MPI_Gatherv(estep.data(), estep.size(), MPI_INT, estep_total.data(), len_list.data(), displacement_list.data(),
                    MPI_INT, 0, this->comm);
        MPI_Gatherv(edw.data(), edw.size(), MPI_FLOAT, edw_total.data(), len_list.data(), displacement_list.data(),
                    MPI_FLOAT, 0, this->comm);

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
                    displacement_list.data(), MPI_FLOAT, 0, this->comm);

        if (rank == 0)
        {
            for (int k = 0; k < igid_total.size(); k++)
            {
                struct Particle iptl;
                iptl.gid = igid_total[k];
                iptl.flag = iflag_total[k];
                iptl.esc_step = istep_total[k];
                iptl.r = GET(iphase_total, k, 0);
                iptl.z = GET(iphase_total, k, 1);
                iptl.phi = GET(iphase_total, k, 2);
                iptl.rho = GET(iphase_total, k, 3);
                iptl.w1 = GET(iphase_total, k, 4);
                iptl.w2 = GET(iphase_total, k, 5);
                iptl.mu = GET(iphase_total, k, 6);
                iptl.w0 = GET(iphase_total, k, 7);
                iptl.f0 = GET(iphase_total, k, 8);
                iptl.psi = GET(iphase_total, k, 9);
                iptl.B = GET(iphase_total, k, 10);
                iptl.dw = idw_total[k];

                int flag1; // tmp flag
                flag1 = iflag_total[k];

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

            for (int k = 0; k < egid_total.size(); k++)
            {
                struct Particle eptl;
                eptl.gid = egid_total[k];
                eptl.flag = eflag_total[k];
                eptl.esc_step = estep_total[k];
                eptl.r = GET(ephase_total, k, 0);
                eptl.z = GET(ephase_total, k, 1);
                eptl.phi = GET(ephase_total, k, 2);
                eptl.rho = GET(ephase_total, k, 3);
                eptl.w1 = GET(ephase_total, k, 4);
                eptl.w2 = GET(ephase_total, k, 5);
                eptl.mu = GET(ephase_total, k, 6);
                eptl.w0 = GET(ephase_total, k, 7);
                eptl.f0 = GET(ephase_total, k, 8);
                eptl.psi = GET(ephase_total, k, 9);
                eptl.B = GET(ephase_total, k, 10);
                eptl.dw = edw_total[k];

                int flag1; // tmp flag
                flag1 = eflag_total[k];

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

            // debug
            std::cout << std::endl;
            std::cout << ">>> Step: " << this->istep << std::endl;
            std::cout << "Num. of escaped ions: " << iesc.size() << std::endl;
            std::cout << "Num. of escaped elec: " << eesc.size() << std::endl;
            std::cout << "Num. of divertor ions: " << idiv.size() << std::endl;
            std::cout << "Num. of divertor elec: " << ediv.size() << std::endl;

            // print first 10 esc particles
            int count = 0;
            t_ParticlesList::iterator it;
            for (it = iesc.begin(); it != iesc.end(); it++)
            {
                printf("iesc gid, rzphi, flag: %lld %f %f %f %d\n", it->second.gid, it->second.r, it->second.z,
                       it->second.phi, it->second.flag);
                count++;
                if (count > 10)
                    break;
            }

            // separate divertor particles and escaped particles
            this->iesc_db.push_back(iesc);
            this->eesc_db.push_back(eesc);

            // Calculate heatload from divertor particles
            HeatLoad ion(1);
            HeatLoad elec(0);

            heatload_calc(idiv, ion, this->iesc_db); // need to send DB
            heatload_calc(ediv, elec, this->eesc_db);
            output(ad, ion, elec);
        }

        this->reader.EndStep();
        this->istep++;
    }
    return status;
}
