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

#define LOG BOOST_LOG_TRIVIAL(debug)

#define NCOL 11
#define GET(X, i, j) X[i * NCOL + j]

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

int read_mesh(adios2::ADIOS &ad, adios2::IO &io)
{
    int n_t;
    adios2::Engine reader;
    io = ad.DeclareIO("'diagnosis.mesh");
    reader = io.Open("xgc.mesh.bp", adios2::Mode::Read, MPI_COMM_SELF);
    auto var = io.InquireVariable<int>("n_t");
    reader.Get<int>(var, &n_t);
    reader.Close();

    return n_t;
}

template <typename T>
std::vector<T>& split_vector(std::vector<T> &vec, int world_size, int rank)
{
    int nblock = vec.size() / world_size;
    int offset = nblock * rank;
    if (rank == world_size - 1)
        nblock = vec.size() - offset;
    
    std::vector<T> sub;
    for (int i=offset; i<offset+nblock; i++) 
    {
        sub.push_back(vec[i]);
    }

    return sub;
}

int main(int argc, char *argv[])
{
    int rank, world_size;
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &world_size);

    if (rank == 0)
        printf("rank,world_size: %d %d\n", rank, world_size);

    if (argc < 1)
    {
        if (rank == 0)
            show_usage(argv[0]);
        return 1;
    }

    init_log(rank);

    long unsigned int nphi = 0;    // number of planes (will be set after reading)
    int np_per_plane = world_size; // number of PEs per plane

    adios2::ADIOS ad(comm);
    adios2::IO io;
    adios2::Engine reader;

    int n_t = read_mesh(ad, io);
    if (rank == 0)
        LOG << "num. of triangles: " << n_t;

    std::vector<double> i_marker_den(n_t);
    std::vector<double> i_dr_avg(n_t);
    std::vector<double> i_En_dr_avg(n_t);
    std::vector<double> i_dr_std(n_t);
    std::vector<double> i_En_dr_std(n_t);

    io = ad.DeclareIO("reader");

    int istep = 0;
    reader = io.Open("xgc.tracer_diag.bp", adios2::Mode::Read, comm);
    while (true)
    {
        adios2::StepStatus status = reader.BeginStep();

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

        if (rank == 0)
            printf("Reading step: %d\n", istep);

        auto var_table = io.InquireVariable<double>("table");
        auto block_list = reader.BlocksInfo(var_table, istep);

        int nblock = block_list.size() / world_size;
        int offset = nblock * rank;
        if (rank == world_size - 1)
            nblock = block_list.size() - offset;
        LOG << boost::format("offset,nblock= %d %d") % offset % nblock;

        // Read table block by block
        for (int i = offset; i < offset + nblock; i++)
        {
            auto block = block_list[i];
            var_table.SetBlockSelection(block.BlockID);
            std::vector<double> table;
            reader.Get<double>(var_table, table);
            reader.PerformGets();

            // Process each row
            int nrow = table.size() / NCOL;
            LOG << "table id,nrow: " << block.BlockID << " " << nrow;
            for (int k = 0; k < nrow; k++)
            {
                int itri = int(GET(table, k, 0));
                double _i_dr_average = GET(table, k, 1);
                double _i_dr_squared_average = GET(table, k, 2);
                double _i_dE_average = GET(table, k, 3);
                double _i_dE_squared_average = GET(table, k, 4);
                double _i_marker_den = GET(table, k, 5);

                double _e_dr_average = GET(table, k, 6);
                double _e_dr_squared_average = GET(table, k, 7);
                double _e_dE_average = GET(table, k, 8);
                double _e_dE_squared_average = GET(table, k, 9);
                double _e_marker_den = GET(table, k, 10);

                i_dr_avg[itri] += _i_dr_average;
                i_marker_den[itri] += _i_marker_den;

                // LOG << boost::format("%d: %d %g %g") % block.BlockID % itri % _i_dr_average % _i_marker_den;
            }
        }

        // Merge all to rank 0
        std::vector<double> vec_list[] = {i_marker_den, i_dr_avg, i_En_dr_avg, i_dr_std, i_En_dr_std};

        for (auto &vec : vec_list)
        {
            if (rank == 0)
            {
                MPI_Reduce(MPI_IN_PLACE, vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            else
            {
                MPI_Reduce(vec.data(), vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
        }

        reader.EndStep();
        istep++;
    }

    reader.Close();

    MPI_Barrier(comm);
    MPI_Finalize();
    return 0;
}
