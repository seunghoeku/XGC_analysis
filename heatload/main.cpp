/*********** XGC *************/

#include <assert.h>
#include <string>
#include <chrono>
#include <thread>

#include "adios2.h"
#include "heatload.hpp"

#define GET(X, i, j) X[i * 9 + j]

// extern "C" void set_test_type(int test_type);

Simulation sml; // input parameters that controls simulation. 
adios2::ADIOS * ad;

int main(int argc, char *argv[]) {
    // Parse command line arguments
    // set_test_type(0);
    if (argc > 2) {
        printf("ERROR: Too many command line arguments. Available options are "
            "'--test', '--update-test', or neither.\n");
        exit(1);
    }
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--test") {
        // set_test_type(1);
        } else if (std::string(argv[i]) == "--update-test") {
        // set_test_type(2);
        } else {
        printf("ERROR: Unknown command line argument. Available options are "
                "'--test', '--update-test', or neither.\n");
        exit(1);
        }
    }

    adios2::ADIOS adios("adios2cfg.xml");
    ad = &adios;
    // run actual routine
    heatload();
}

