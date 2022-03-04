/*********** XGC *************/

#include <assert.h>
#include <chrono>
#include <string>
#include <thread>

#include "adios2.h"
#include "heatload.hpp"

#define GET(X, i, j) X[i * 9 + j]

// MPI color for MPMD mode: Diffusion(3), Heatload(5), Poincare (7), Middleman(9)
#define MY_COLOR 5

// extern "C" void set_test_type(int test_type);

int main(int argc, char *argv[])
{
    int world_rank, world_size;
    int rank, comm_size;
    MPI_Comm comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // MPI comm split for MPMD mode:
    MPI_Comm_split(MPI_COMM_WORLD, MY_COLOR, world_rank, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    if (rank == 0)
        printf("color,world_rank,world_size,rank,comm_size: %d %d %d %d\n", MY_COLOR, world_rank, world_size, rank,
               comm_size);

    // Parse command line arguments
    // set_test_type(0);
    if (argc > 2)
    {
        printf("ERROR: Too many command line arguments. Available options are "
               "'--test', '--update-test', or neither.\n");
        exit(1);
    }
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--test")
        {
            // set_test_type(1);
        }
        else if (std::string(argv[i]) == "--update-test")
        {
            // set_test_type(2);
        }
        else
        {
            printf("ERROR: Unknown command line argument. Available options are "
                   "'--test', '--update-test', or neither.\n");
            exit(1);
        }
    }

    adios2::ADIOS *ad = new adios2::ADIOS("adios2cfg.xml", comm);

    // run actual routine
    heatload_init(ad, comm);
    int istep = 1;
    while (1)
    {
        int ret = heatload_step(ad, istep);

        if (ret == 0)
            // everything is ok
            istep++;
        else if (ret > 0)
            // wait again
            continue;
        else if (ret == -1)
            // no more data
            break;
    }
    heatload_finalize();

    delete ad;
    MPI_Finalize();
    return 0;
}
