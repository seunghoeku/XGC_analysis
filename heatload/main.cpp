/*********** XGC *************/

#include <assert.h>
#include <chrono>
#include <string>
#include <thread>

#include "adios2.h"
#include "heatload.hpp"

#define GET(X, i, j) X[i * 9 + j]

// extern "C" void set_test_type(int test_type);

int main(int argc, char *argv[])
{
    int rank, world_size;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

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

    adios2::ADIOS ad("adios2cfg.xml");
    // run actual routine
    heatload_init(&ad, comm);
    int istep = 1;
    while (1)
    {
        int ret = heatload_step(&ad, istep);

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
    heatload_finalize(&ad);

    MPI_Finalize();
    return 0;
}
