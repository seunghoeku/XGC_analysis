#include "adios2.h"
#include "heatload_calc.hpp"
#include "mpi.h"
#include "sml.hpp"

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

extern Simulation sml;

void init(adios2::ADIOS *ad, std::string xgcdir)
{
    adios2::Engine reader;
    adios2::IO reader_io;

    reader_io = ad->DeclareIO("output.units"); // same IO name as in XGC
    boost::filesystem::path fname = boost::filesystem::path(xgcdir) / boost::filesystem::path("xgc.units.bp");
    LOG << "Loading: " << fname;
    reader = reader_io.Open(fname.string(), adios2::Mode::Read, MPI_COMM_SELF);

    double sml_e_charge, sml_prot_mass;
    double ptl_mass_au, ptl_charge_eu, ptl_e_mass_au, ptl_e_charge_eu;
    reader.Get<double>("eq_x_psi", sml.psix);
    reader.Get<double>("eq_x_r", sml.x_r);
    reader.Get<double>("eq_x_z", sml.x_z);
    reader.Get<double>("eq_axis_r", sml.axis_r);
    reader.Get<double>("eq_axis_z", sml.axis_z);
    reader.Get<double>("sml_dt", sml.dt);
    reader.Get<double>("sml_e_charge", sml_e_charge);
    reader.Get<double>("sml_prot_mass", sml_prot_mass);
    reader.Get<double>("ptl_mass_au", ptl_mass_au);
    reader.Get<double>("ptl_charge_eu", ptl_charge_eu);
    reader.Get<double>("ptl_e_mass_au", ptl_e_mass_au);
    reader.Get<double>("ptl_e_charge_eu", ptl_e_charge_eu);
    reader.Get<double>("diag_heat_rmin1", sml.rmin[0]);
    reader.Get<double>("diag_heat_rmax1", sml.rmax[0]);
    reader.Get<double>("diag_heat_zmin1", sml.zmin[0]);
    reader.Get<double>("diag_heat_zmax1", sml.zmax[0]);
    reader.Get<double>("diag_heat_pmin1", sml.pmin[0]);
    reader.Get<double>("diag_heat_pmax1", sml.pmax[0]);
    reader.Get<double>("diag_heat_rmin2", sml.rmin[1]);
    reader.Get<double>("diag_heat_rmax2", sml.rmax[1]);
    reader.Get<double>("diag_heat_zmin2", sml.zmin[1]);
    reader.Get<double>("diag_heat_zmax2", sml.zmax[1]);
    reader.Get<double>("diag_heat_pmin2", sml.pmin[1]);
    reader.Get<double>("diag_heat_pmax2", sml.pmax[1]);
    reader.Close();

    sml.npsi = N_PSI;
    sml.ncond = N_COND;

    sml.c2_2m[0] =
        (sml_e_charge * ptl_e_charge_eu) * (sml_e_charge * ptl_e_charge_eu) * 0.5 / (ptl_e_mass_au * sml_prot_mass);
    sml.c2_2m[1] =
        (sml_e_charge * ptl_charge_eu) * (sml_e_charge * ptl_charge_eu) * 0.5 / (ptl_mass_au * sml_prot_mass);

    sml.dpsi[0] = (sml.pmax[0] - sml.pmin[0]) / static_cast<double>(sml.npsi - 1);
    sml.dpsi[1] = (sml.pmax[1] - sml.pmin[1]) / static_cast<double>(sml.npsi - 1);

    // get poloidal angle of X point
    sml.x_theta = atan2(sml.x_z - sml.axis_z, sml.x_r - sml.axis_r);
    if (sml.x_r - sml.axis_r < 0)
    {
        sml.x_theta = sml.x_theta + M_PI;
    }
    sml.ntheta = N_THETA;
    sml.dtheta = 2. * M_PI / (double)N_THETA;
}