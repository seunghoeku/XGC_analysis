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

#include "diffusion.hpp"
#include "heatloadcalc.hpp"

#define LOG BOOST_LOG_TRIVIAL(debug)

// MPI color for MPMD mode: Diffusion(3), Heatload(5), Poincare (7), Middleman(9)
#define MY_COLOR 9

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

    if (rank == 0)
        printf("world_rank,world_size,rank,comm_size: %d %d %d %d\n", world_rank, world_size, rank, comm_size);

    if (argc < 1)
    {
        if (rank == 0)
            show_usage(argv[0]);
        return 1;
    }

    init_log(rank);

    adios2::ADIOS *ad = new adios2::ADIOS("adios2cfg.xml", comm);

    // Comm for heatload
    MPI_Group comm_group;
    MPI_Comm_group(comm, &comm_group);

    std::vector<int> ranks(comm_size);
    for (int i = 0; i < comm_size; i++)
        ranks[i] = (i + 1) % comm_size;

    MPI_Group heatload_group;
    MPI_Group_incl(comm_group, comm_size, ranks.data(), &heatload_group);

    MPI_Comm heatload_comm;
    MPI_Comm_create(comm, heatload_group, &heatload_comm);

    Diffusion diffusion(ad, comm);
    Heatload heatload(ad, heatload_comm);

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

        status = heatload.step();

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
        if (istep > 2)
            break;
    }

    diffusion.finalize();
    // heatload.finalize();

    delete ad;
    MPI_Finalize();
    return 0;
}
