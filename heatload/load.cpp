
#include <assert.h>

#include <string>

#include "adios2.h"

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

adios2::StepStatus load_data(std::vector<long> &igid, std::vector<int> &iflag, std::vector<float> &idw,
                             std::vector<float> &iphase, std::vector<long> &egid, std::vector<int> &eflag,
                             std::vector<float> &edw, std::vector<float> &ephase)
{
    adios2::StepStatus status = reader.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
        // Clear vector
        igid.clear();
        egid.clear();
        iphase.clear();
        ephase.clear();

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

        reader.Get<long>(var_igid, igid);
        reader.Get<long>(var_egid, egid);
        reader.Get<int>(var_iflag, iflag);
        reader.Get<int>(var_eflag, eflag);
        reader.Get<float>(var_idw, idw);
        reader.Get<float>(var_edw, edw);
        reader.Get<float>(var_iphase, iphase);
        reader.Get<float>(var_ephase, ephase);
        reader.EndStep();
    }

    return status;
}
