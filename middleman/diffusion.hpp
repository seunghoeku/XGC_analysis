#include <iostream>
#include <vector>

#include "adios2.h"
#include "mpi.h"

class Diffusion
{
  public:
    Diffusion(adios2::ADIOS *ad, MPI_Comm comm);

    void finalize();

    adios2::StepStatus step();

  public:
    adios2::ADIOS *ad;
    adios2::IO io;
    adios2::Engine reader;

    MPI_Comm comm;
    int comm_size;
    int rank;

    int ntriangle;
    int istep;

    std::vector<double> i_marker_den;
    std::vector<double> i_dr_avg;
    std::vector<double> i_En_dr_avg;
    std::vector<double> i_dr_std;
    std::vector<double> i_En_dr_std;
};
