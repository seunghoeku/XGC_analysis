/*********** XGC *************/

#include <assert.h>
#include <string>

#include "adios2.h"
#include "flags.hpp"
#include "load.hpp"
#include "particles.hpp"
#include "sml.hpp"
#include "particles.hpp"

#define GET(X, i, j) X[i * 9 + j]

void heatload();
void init(); // initialization
void heatload_calc(std::vector<Particles> ediv ); // calculate heatload
void output(); // output graphs or data for graphs

// extern "C" void set_test_type(int test_type);

Simulation sml; // input parameters that controls simulation. 


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

    // run actual routine
    heatload();
    }

void heatload() {

    // init simulation parameters
    init();

    // init adios
    load_init("xgc.escaped_ptls.su455.bp");

    int i = 0;
    while (1) {
        i++;

        t_ParticlesList iptls;
        t_ParticlesList eptls;

        adios2::StepStatus status = load_data(iptls, eptls);
        if (status != adios2::StepStatus::OK)
            break;

        std::cout << ">>> Step: " << i << std::endl;
        std::cout << "Num. of ions: " << iptls.size() << std::endl;
        std::cout << "Num. of eons: " << eptls.size() << std::endl;

        // print first 10 particle
        for (int i = 0; i < 10; i++) {
            printf("iptl gid, rzphi, flag: %ld %f %f %f %d\n", iptls[i].gid, iptls[i].ph.r, iptls[i].ph.z, iptls[i].ph.phi, iptls[i].flag);
        }

        std::vector<Particles> idiv;
        std::vector<Particles> ediv;
        std::vector<Particles> iesc;
        std::vector<Particles> eesc;

        // separate divertor particles and escaped particles

        // store escaped particles to DB

        // Calculate heatload from divertor particles
        double iheat[sml.ncond][sml.npsi];
        double eheat[sml.ncond][sml.npsi];

        heatload_calc(idiv); // need to send DB
        heatload_calc(ediv);
        output();
    }

    load_finalize();
}

void output() {}