/*********** XGC *************/

#include <assert.h>
#include <chrono>
#include <string>
#include <thread>

#include "adios2.h"
#include "heatload.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/program_options.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

#define GET(X, i, j) X[i * 9 + j]

// MPI color for MPMD mode: Diffusion(3), Heatload(5), Poincare (7), Middleman(9)
#define MY_COLOR 5

namespace po = boost::program_options;

// extern "C" void set_test_type(int test_type);

void init_log(int rank)
{
    std::string fmt = boost::str(boost::format("%d: %%Message%%") % rank);

    // Output message to console
    boost::log::add_console_log(std::cout, boost::log::keywords::format = fmt, boost::log::keywords::auto_flush = true);

    boost::log::add_common_attributes();
}

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

    init_log(rank);

    std::string xgcdir = ".";
    int maxstep = 0;
    bool freshstart = false;
    bool ion_only = false;

    po::options_description desc("Allowed options");
    desc.add_options()
        // clang-format off
        ("help,h", "produce help message")
        ("xgcdir,w", po::value(&xgcdir), "XGC directory")
        ("maxstep,s", po::value<int>(&maxstep)->default_value(0), "max steps")
        ("freshstart,f", po::bool_switch(&freshstart), "fresh start (no restart)")
        ("ion_only,i", po::bool_switch(&ion_only), "Ion only")
        // clang-format on
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (rank == 0)
    {
        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return 1;
        }
        printf("color,world_rank,world_size,rank,comm_size: %d %d %d %d %d\n", MY_COLOR, world_rank, world_size, rank,
               comm_size);
        LOG << "XGC directory: " << xgcdir;
    }

    adios2::ADIOS *ad = new adios2::ADIOS("adios2cfg.xml", comm);

    // run actual routine
    heatload_init(ad, comm, xgcdir, !freshstart);
    int istep = 1;
    while (1)
    {
        int ret = heatload_step(ad, istep, ion_only);

        if (ret == 0)
            // everything is ok
            ;
        else if (ret > 0)
            // wait again
            continue;
        else if (ret == -1)
            // no more data
            break;

        istep++;
        if ((maxstep > 0) && (istep > maxstep))
            break;
    }

    heatload_finalize();

    delete ad;
    MPI_Finalize();
    return 0;
}
