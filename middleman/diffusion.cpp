#include "diffusion.hpp"
#include "util.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

#define NCOL 11
#define GET(X, i, j) X[i * NCOL + j]

inline int read_mesh(adios2::ADIOS *ad, adios2::IO &io)
{
    int n_t;
    adios2::Engine reader;
    io = ad->DeclareIO("'diagnosis.mesh");
    reader = io.Open("xgc.mesh.bp", adios2::Mode::Read, MPI_COMM_SELF);
    auto var = io.InquireVariable<int>("n_t");
    reader.Get<int>(var, &n_t);
    reader.Close();

    return n_t;
}

Diffusion::Diffusion(adios2::ADIOS *ad, MPI_Comm comm)
{
    this->ad = ad;

    this->comm = comm;
    MPI_Comm_rank(comm, &this->rank);
    MPI_Comm_size(comm, &this->comm_size);

    this->ntriangle = read_mesh(ad, this->io);
    this->istep = 0;

    this->i_marker_den.resize(this->ntriangle);
    this->i_dr_avg.resize(this->ntriangle);
    this->i_En_dr_avg.resize(this->ntriangle);
    this->i_dr_std.resize(this->ntriangle);
    this->i_En_dr_std.resize(this->ntriangle);

    this->reader = this->io.Open("xgc.tracer_diag.bp", adios2::Mode::Read, this->comm);
}

void Diffusion::finalize()
{
    this->reader.Close();
}

adios2::StepStatus Diffusion::step()
{
    adios2::StepStatus status = this->reader.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
        auto var_table = this->io.InquireVariable<double>("table");
        auto block_list = reader.BlocksInfo(var_table, this->istep);

        auto slice = split_vector(block_list, this->comm_size, this->rank);
        LOG << boost::format("offset,nblock= %d %d") % slice.first % slice.second;

        int offset = slice.first;
        int nblock = slice.second;

        // Read table block by block
        for (int i = offset; i < offset + nblock; i++)
        {
            auto block = block_list[i];
            var_table.SetBlockSelection(block.BlockID);
            std::vector<double> table;
            this->reader.Get<double>(var_table, table);
            this->reader.PerformGets();

            // Process each row
            int nrow = table.size() / NCOL;
            //LOG << "table id,nrow: " << block.BlockID << " " << nrow;
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

                this->i_dr_avg[itri] += _i_dr_average;
                this->i_marker_den[itri] += _i_marker_den;

                // LOG << boost::format("%d: %d %g %g") % block.BlockID % itri % _i_dr_average % _i_marker_den;
            }

            // Merge all to rank 0
            std::vector<double> vec_list[] = {this->i_marker_den, this->i_dr_avg, this->i_En_dr_avg, this->i_dr_std,
                                              this->i_En_dr_std};

            for (auto &vec : vec_list)
            {
                if (this->rank == 0)
                {
                    MPI_Reduce(MPI_IN_PLACE, vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, 0, this->comm);
                }
                else
                {
                    MPI_Reduce(vec.data(), vec.data(), vec.size(), MPI_DOUBLE, MPI_SUM, 0, this->comm);
                }
            }
        }

        // Save
        if (this->rank == 0)
        {
            // save output
        }

        this->reader.EndStep();
        this->istep++;
    }
    return status;
}
