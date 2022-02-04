/*********** XGC *************/

#include <assert.h>
#include <string>

#include "adios2.h"
#include "flags.hpp"
#include "load.hpp"
#include "particles.hpp"
#include "sml.hpp"

#define GET(X, i, j) X[i * 9 + j]

void heatload();
void init(Simulation &sml); // initialization
void heatload_calc(std::vector<Particles> idiv,
                   std::vector<Particles> ediv); // calculate heatload
void output(); // output graphs or data for graphs

// extern "C" void set_test_type(int test_type);

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

  Simulation sml;

  // init simulation parameters
  init(sml);

  // init adios
  load_init("xgc.escaped_ptls.su455.bp");
  // adios2::ADIOS ad;
  // adios2::IO reader_io = ad.DeclareIO("headload");
  // adios2::Engine reader =
  //     reader_io.Open("xgc.escaped_ptls.su455.bp", adios2::Mode::Read);
  std::vector<long> igid;
  std::vector<long> egid;
  std::vector<int> iflag;
  std::vector<int> eflag;
  std::vector<float> idw;
  std::vector<float> edw;
  std::vector<float> iphase;
  std::vector<float> ephase;

  int i = 0;
  while (1) {
    i++;
    adios2::StepStatus status =
        load_data(igid, iflag, idw, iphase, egid, eflag, edw, ephase);
    if (status != adios2::StepStatus::OK)
      break;
    std::cout << ">>> Step: " << i << std::endl;

    std::cout << "Num. of ions: " << igid.size() << std::endl;
    std::cout << "Num. of eons: " << egid.size() << std::endl;
    assert(iphase.size() / igid.size() == 9);
    assert(ephase.size() / egid.size() == 9);

    // print first 10 particle
    for (int i = 0; i < 10; i++) {
      printf("Particle gid, rzphi: %ld %f %f %f\n", igid[i], GET(iphase, i, 0),
             GET(iphase, i, 1), GET(iphase, i, 2));
    }

    std::vector<Particles> idiv;
    std::vector<Particles> ediv;
    std::vector<Particles> iesc;
    std::vector<Particles> eesc;

    // separate divertor particles and escaped particles

    // store escaped particles to DB

    // Calculate heatload from divertor particles
    heatload_calc(idiv, ediv); // need to send DB

    output();
  }

  load_finalize();
}

void output() {}