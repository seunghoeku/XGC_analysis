#include "output.hpp"

// Example:
// 2D array of shape (d0, d1): index at (i,j) is (d1*i + j)
// 3D array of shape (d0, d1, d2): index at (i,j,k) is ((d1+d2)*i + d2*j + k)
#define GET2D(X, d1, i, j) X[d1 * i + j]
#define GET3D(X, d1, d2, i, j, k) X[(d1 * d2) * i + d2 * j + k]

adios2::IO output_io;
adios2::Engine writer;

extern Simulation sml;

void output(adios2::ADIOS *ad, HeatLoad &ion, HeatLoad &elec, MPI_Comm comm)
{
    // (2022/03/16) TODO: parallel output
    int comm_size, rank;
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &rank);

    static bool first = true;

    if (first)
    {
        output_io = ad->DeclareIO("heatload");
        output_io.DefineVariable<double>("psi", {N_SIDE, N_PSI}, {0, 0}, {N_SIDE, N_PSI});
        output_io.DefineVariable<double>("ienflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("iptlflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("eenflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});
        output_io.DefineVariable<double>("eptlflux", {N_SIDE, N_COND, N_PSI}, {0, 0, 0}, {N_SIDE, N_COND, N_PSI});

        output_io.DefineVariable<double>("io.side", {N_SIDE + 1}, {0}, {N_SIDE + 1});

        writer = output_io.Open("xgc.heatload.bp", adios2::Mode::Write, comm);

        first = false;
    }

    double psi[N_SIDE * N_PSI];

    for (int is = 0; is < N_SIDE; is++)
    {
        for (int i = 0; i < N_PSI; i++)
        {
            // psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
            GET2D(psi, N_PSI, is, i) = sml.pmin[is] + sml.dpsi[is] * (double)(i);
        }
    }

    double ienflux[N_SIDE * N_COND * N_PSI];
    double iptlflux[N_SIDE * N_COND * N_PSI];
    double eenflux[N_SIDE * N_COND * N_PSI];
    double eptlflux[N_SIDE * N_COND * N_PSI];

    for (int is = 0; is < N_SIDE; is++)
    {
        for (int ic = 0; ic < N_COND; ic++)
        {
            for (int i = 0; i < N_PSI; i++)
            {
                // psi[is][i] = sml.pmin[is] + sml.dpsi[is] * (double)(i);
                GET3D(ienflux, N_COND, N_PSI, is, ic, i) = ion.side[is].en[ic][i];
                GET3D(iptlflux, N_COND, N_PSI, is, ic, i) = ion.side[is].ptl[ic][i];
                GET3D(eenflux, N_COND, N_PSI, is, ic, i) = elec.side[is].en[ic][i];
                GET3D(eptlflux, N_COND, N_PSI, is, ic, i) = elec.side[is].ptl[ic][i];
            }
        }
    }

    // save psi, ion.side[0:N_SIDE].en[0:N_COND][0:N_PSI] and ptl[0:NCOND][0:N_PSI] and electron.
    writer.BeginStep();
    writer.Put<double>("psi", psi);
    writer.Put<double>("ienflux", ienflux);
    writer.Put<double>("iptlflux", iptlflux);
    writer.Put<double>("eenflux", eenflux);
    writer.Put<double>("eptlflux", eptlflux);
    writer.EndStep();
}

void output_finalize()
{
    writer.Close();
}
