/*********** XGC *************/

#include <assert.h>
#include <string>

#include "adios2.h"
#include "particles.hpp"
#include "flags.hpp"

#define GET(X, i, j) X[i * 9 + j]

void heatload();
void init(); // initialization
// receving data from XGC1
void load_data(adios2::IO reader_io, adios2::Engine reader,
               std::vector<long> &igid, std::vector<float> &iphase,
               std::vector<long> &egid, std::vector<float> &ephase);
void heatload_calc(std::vector<Particles> idiv, std::vector<Particles> ediv); // calculate heatload
void output();        // output graphs or data for graphs

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
  init();

  // init adios
  adios2::ADIOS ad;
  adios2::IO reader_io = ad.DeclareIO("headload");
  adios2::Engine reader =
      reader_io.Open("xgc.escaped_ptls.su455.bp", adios2::Mode::Read);
  std::vector<long> igid;
  std::vector<long> egid;
  std::vector<float> iphase;
  std::vector<float> ephase;

  int i = 0;
  while (1) {
    adios2::StepStatus status = reader.BeginStep();
    if (status != adios2::StepStatus::OK)
      break;
    std::cout << ">>> Step:" << i << std::endl;
    load_data(reader_io, reader, igid, iphase, egid, ephase);
    reader.EndStep();
    i++;

    std::cout << "Num. of ions: " << igid.size() << std::endl;
    std::cout << "Num. of eons: " << egid.size() << std::endl;
    assert(iphase.size() / igid.size() == 9);
    assert(ephase.size() / egid.size() == 9);

    // print first 10 particle
    for (int i = 0; i < 10; i++) {
      printf("Particle gid, rzphi: %ld %f %f %f\n", igid[i], GET(iphase, i, 0),
             GET(iphase, i, 1), GET(iphase, i, 2));
    }
  }
  reader.Close();

  std::vector<Particles> idiv;
  std::vector<Particles> ediv;
  std::vector<Particles> iesc;
  std::vector<Particles> eesc;

  //separate divertor particles and escaped particles

  // store escaped particles to DB

  // Calculate heatload from divertor particles 
  heatload_calc(idiv,ediv); // need to send DB

  output();
}

void init() {}

void output() {}