
#include <assert.h>

#include <string>

#include "flags.hpp"
#include "load.hpp"
#include "sml.hpp"

#define GET(X, i, j) X[i * 9 + j]

adios2::ADIOS ad;
adios2::Engine reader;
adios2::IO reader_io;

void load_init(const std::string &filename)
{
    reader_io = ad.DeclareIO("headload");
    reader = reader_io.Open(filename, adios2::Mode::Read);
}

void load_finalize()
{
    reader.Close();
}

adios2::StepStatus load_data(t_ParticlesList &idiv, t_ParticlesList &ediv, t_ParticlesList &iesc, t_ParticlesList &eesc)
{
    // Clear vector
    idiv.clear();
    ediv.clear();
    iesc.clear();
    eesc.clear();

    adios2::StepStatus status = reader.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
        // Inquire variables
        auto var_igid = reader_io.InquireVariable<long>("igid");
        auto var_egid = reader_io.InquireVariable<long>("egid");
        auto var_iflag = reader_io.InquireVariable<int>("iflag");
        auto var_eflag = reader_io.InquireVariable<int>("eflag");
        auto var_idw = reader_io.InquireVariable<float>("idw");
        auto var_edw = reader_io.InquireVariable<float>("edw");
        auto var_iphase = reader_io.InquireVariable<float>("iphase");
        auto var_ephase = reader_io.InquireVariable<float>("ephase");

        var_igid.SetSelection({{0}, {var_igid.Shape()[0]}});
        var_egid.SetSelection({{0}, {var_egid.Shape()[0]}});
        var_iflag.SetSelection({{0}, {var_iflag.Shape()[0]}});
        var_eflag.SetSelection({{0}, {var_eflag.Shape()[0]}});
        var_idw.SetSelection({{0}, {var_idw.Shape()[0]}});
        var_edw.SetSelection({{0}, {var_edw.Shape()[0]}});
        var_iphase.SetSelection({{0, 0}, {var_iphase.Shape()[0], var_iphase.Shape()[1]}});
        var_ephase.SetSelection({{0, 0}, {var_ephase.Shape()[0], var_ephase.Shape()[1]}});

        std::vector<long> igid;
        std::vector<long> egid;
        std::vector<int> iflag;
        std::vector<int> eflag;
        std::vector<float> idw;
        std::vector<float> edw;
        std::vector<float> iphase;
        std::vector<float> ephase;

        reader.Get<long>(var_igid, igid);
        reader.Get<long>(var_egid, egid);
        reader.Get<int>(var_iflag, iflag);
        reader.Get<int>(var_eflag, eflag);
        reader.Get<float>(var_idw, idw);
        reader.Get<float>(var_edw, edw);
        reader.Get<float>(var_iphase, iphase);
        reader.Get<float>(var_ephase, ephase);
        reader.EndStep();

        assert(iphase.size() / igid.size() == 9);
        assert(ephase.size() / egid.size() == 9);

        // populate particles
        for (int i = 0; i < igid.size(); i++)
        {
            Particles iptl;
            iptl.gid = igid[i];
            iptl.flag = iflag[i];
            iptl.ph.r = GET(iphase, i, 0);
            iptl.ph.z = GET(iphase, i, 1);
            iptl.ph.phi = GET(iphase, i, 2);
            iptl.ph.rho = GET(iphase, i, 3);
            iptl.ph.w1 = GET(iphase, i, 4);
            iptl.ph.w2 = GET(iphase, i, 5);
            iptl.ph.mu = GET(iphase, i, 6);
            iptl.ph.w0 = GET(iphase, i, 7);
            iptl.ph.f0 = GET(iphase, i, 8);

            int flag1; // tmp flag
            flag1 = iflag[i];

            Flags fl(flag1); // decode flags

            // save to div or esc
            if (fl.escaped)
            {
                // add to esc
                iesc.insert(std::pair<long long, Particles>(iptl.gid, iptl));
            }
            else
            {
                // add to div
                idiv.insert(std::pair<long long, Particles>(iptl.gid, iptl));
            }
        }

        for (int i = 0; i < egid.size(); i++)
        {
            Particles eptl;
            eptl.gid = egid[i];
            eptl.flag = eflag[i];
            eptl.ph.r = GET(ephase, i, 0);
            eptl.ph.z = GET(ephase, i, 1);
            eptl.ph.phi = GET(ephase, i, 2);
            eptl.ph.rho = GET(ephase, i, 3);
            eptl.ph.w1 = GET(ephase, i, 4);
            eptl.ph.w2 = GET(ephase, i, 5);
            eptl.ph.mu = GET(ephase, i, 6);
            eptl.ph.w0 = GET(ephase, i, 7);
            eptl.ph.f0 = GET(ephase, i, 8);

            int flag1; // tmp flag
            flag1 = eflag[i];

            Flags fl(flag1); // decode flags

            // save to div or esc
            if (fl.escaped)
            {
                // add to esc
                eesc.insert(std::pair<long long, Particles>(eptl.gid, eptl));
            }
            else
            {
                // add to div
                ediv.insert(std::pair<long long, Particles>(eptl.gid, eptl));
            }
        }
    }

    return status;
}
