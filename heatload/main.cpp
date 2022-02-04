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
void heatload_calc(t_ParticlesList ediv ); // calculate heatload
void output(); // output graphs or data for graphs
void separate(    std::vector<long> igid, std::vector<int> iflag, std::vector<float> idw, 
    std::vector<float> iphase,  std::vector<Particles> idiv,  std::vector<Particles> esc);

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

    t_ParticleDB iesc_db; 

    int i = 0;
    while (1) {
        i++;

        t_ParticlesList idiv;
        t_ParticlesList ediv;
        t_ParticlesList iesc;
        t_ParticlesList eesc;

        adios2::StepStatus status = load_data(idiv, ediv, iesc, eesc);
        if (status != adios2::StepStatus::OK)
            break;

        std::cout << ">>> Step: " << i << std::endl;
        std::cout << "Num. of escaped ions: " << iesc.size() << std::endl;
        std::cout << "Num. of escaped elec: " << eesc.size() << std::endl;
        std::cout << "Num. of divertor ions: " << idiv.size() << std::endl;
        std::cout << "Num. of divertor elec: " << ediv.size() << std::endl;

        // print first 10 esc particles
        int count = 0;
        std::map<long long, Particles>::iterator it;
        for (it = iesc.begin(); it != iesc.end(); it++) {
            printf("iesc gid, rzphi, flag: %lld %f %f %f %d\n", it->second.gid, it->second.ph.r, it->second.ph.z, it->second.ph.phi, it->second.flag);
            count++;
            if (count>10) break;
        }


        // separate divertor particles and escaped particles
        iesc_db.push_back(iesc);
        Particle ptl = search(iesc_db, i-1, 15824414);
        printf ("Found or not? gid=%lld\n", ptl.gid);

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
