#include <iostream>
#include <vector>

#include "adios2.h"
#include "mpi.h"

class Diffusion
{
  public:
    Diffusion(adios2::ADIOS *ad, MPI_Comm comm);

    void finalize();
    void output();

    adios2::StepStatus step();

    // helper
    void vec_reduce(std::vector<double> &vec);

  public:
    adios2::ADIOS *ad;
    adios2::IO io;
    adios2::Engine reader;

    adios2::IO output_io;
    adios2::Engine writer;

    MPI_Comm comm;
    int comm_size;
    int rank;

    int ntriangle;
    int istep;

    std::vector<double> i_dr_avg;
    std::vector<double> i_dr_squared_average;
    std::vector<double> i_dE_avg;
    std::vector<double> i_dE_squared_average;
    std::vector<double> i_marker_den;
    std::vector<double> e_dr_avg;
    std::vector<double> e_dr_squared_average;
    std::vector<double> e_dE_avg;
    std::vector<double> e_dE_squared_average;
    std::vector<double> e_marker_den;
};
