#include <chrono>
#include <iostream>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#include <assert.h>

#include "adios2.h"
#include "mpi.h"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>
#include <boost/program_options.hpp>

#include "diffusion.hpp"

#define LOG BOOST_LOG_TRIVIAL(debug)

// MPI color for MPMD mode: Diffusion(3), Heatload(5), Poincare (7), Middleman(9)
#define MY_COLOR 3

namespace po = boost::program_options;

static void show_usage(std::string name)
{
    std::cerr << "Usage: " << name << std::endl;
}

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

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce help message")("xgcdir,w", po::value(&xgcdir), "XGC directory")(
        "maxstep,s", po::value<int>(&maxstep)->default_value(0), "max steps");

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

    Diffusion diffusion(ad, xgcdir, comm);

    int istep = 1;
    while (true)
    {
        adios2::StepStatus status;

        status = diffusion.step();

        if (status == adios2::StepStatus::EndOfStream)
        {
            std::cout << "Input stream terminated. Exit loop" << std::endl;
            break;
        }
        else if (status == adios2::StepStatus::NotReady)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (status == adios2::StepStatus::OtherError)
        {
            std::cout << "Input stream had errors. Exit loop" << std::endl;
            break;
        }

        istep++;
        if ((maxstep > 0) && (istep > maxstep))
            break;
    }

    diffusion.finalize();

    delete ad;
    MPI_Finalize();
    return 0;
}
