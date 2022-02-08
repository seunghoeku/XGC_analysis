/*********** XGC *************/

#include <assert.h>
#include <string>

#include "adios2.h"
#include "flags.hpp"
#include "load.hpp"
#include "particles.hpp"
#include "sml.hpp"
#include "particles.hpp"
#include "heatload.hpp"

#define GET(X, i, j) X[i * 9 + j]

void heatload();
void init(); // initialization
void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db); // calculate heatload
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
    load_init("xgc.escaped_ptls.bp");

    t_ParticleDB iesc_db; 
    t_ParticleDB eesc_db; 

    int i = 0;
    while (1) {
        i++;

        Particles idiv;
        Particles ediv;
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
        std::map<long long, Particle>::iterator it;
        for (it = iesc.begin(); it != iesc.end(); it++) {
            printf("iesc gid, rzphi, flag: %lld %f %f %f %d\n", it->second.gid, it->second.r, it->second.z, it->second.phi, it->second.flag);
            count++;
            if (count>10) break;
        }


        // separate divertor particles and escaped particles
        iesc_db.push_back(iesc);
        eesc_db.push_back(eesc);
        Particle ptl = search(iesc_db, i-1, 15824414);
        printf ("Found or not? gid=%lld\n", ptl.gid);

        // store escaped particles to DB

        // Calculate heatload from divertor particles
        HeatLoad ion;
        HeatLoad elec;

        heatload_calc(idiv, ion,  iesc_db); // need to send DB
        heatload_calc(ediv, elec, eesc_db);
        output();
    }

    load_finalize();
}

void output() {}
