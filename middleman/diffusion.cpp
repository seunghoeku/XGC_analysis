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
    io = ad->DeclareIO("diagnosis.mesh");
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

    this->i_dr_avg.resize(this->ntriangle);
    this->i_dr_squared_average.resize(this->ntriangle);
    this->i_dE_avg.resize(this->ntriangle);
    this->i_dE_squared_average.resize(this->ntriangle);
    this->i_marker_den.resize(this->ntriangle);
    this->e_dr_avg.resize(this->ntriangle);
    this->e_dr_squared_average.resize(this->ntriangle);
    this->e_dE_avg.resize(this->ntriangle);
    this->e_dE_squared_average.resize(this->ntriangle);
    this->e_marker_den.resize(this->ntriangle);

    this->io = ad->DeclareIO("tracer_diag");
    this->reader = this->io.Open("xgc.tracer_diag.bp", adios2::Mode::Read, this->comm);
}

void Diffusion::finalize()
{
    this->reader.Close();
    if (this->rank == 0)
        this->writer.Close();
}

adios2::StepStatus Diffusion::step()
{
    adios2::StepStatus status = this->reader.BeginStep();
    if (status == adios2::StepStatus::OK)
    {
        auto var_table = this->io.InquireVariable<double>("table");
        auto block_list = reader.BlocksInfo(var_table, this->istep);

        auto slice = split_vector(block_list, this->comm_size, this->rank);
        LOG << boost::format("Diffusion offset,nblock= %d %d") % slice.first % slice.second;

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

            // Process each row:
            // Each row of the "table" contains the following info:
            // triangle, i_dr_average, i_dr_squared_average, i_dE_average, i_dE_squared_average, i_marker_den,
            // e_dr_average, e_dr_squared_average, e_dE_average, e_dE_squared_average, e_marker_den
            int nrow = table.size() / NCOL;
            // LOG << "table id,nrow: " << block.BlockID << " " << nrow;
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

                // LOG << boost::format("%d: %d %g %g") % block.BlockID % itri % _i_dr_average % _i_marker_den;
                this->i_dr_avg[itri] += _i_dr_average;
                this->i_dr_squared_average[itri] += _i_dr_squared_average;
                this->i_dE_avg[itri] += _i_dE_average;
                this->i_dE_squared_average[itri] += _i_dE_squared_average;
                this->e_marker_den[itri] += _e_marker_den;
                this->e_dr_avg[itri] += _e_dr_average;
                this->e_dr_squared_average[itri] += _e_dr_squared_average;
                this->e_dE_avg[itri] += _e_dE_average;
                this->e_dE_squared_average[itri] += _e_dE_squared_average;
                this->e_marker_den[itri] += _e_marker_den;
            }

            // Merge all to rank 0
            std::vector<double> vec_list[] = {
                this->i_dr_avg,
                this->i_dr_squared_average,
                this->i_dE_avg,
                this->i_dE_squared_average,
                this->i_marker_den,
                this->e_dr_avg,
                this->e_dr_squared_average,
                this->e_dE_avg,
                this->e_dE_squared_average,
                this->e_marker_den,
            };

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
            this->output();
        }

        this->reader.EndStep();
        this->istep++;
    }
    return status;
}

void Diffusion::output()
{
    static bool first = true;

    if (first)
    {
        this->output_io = ad->DeclareIO("diffusion");
        long unsigned int ntri = this->ntriangle;

        this->output_io.DefineVariable<double>("i_dr_avg", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("i_dr_squared_average", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("i_dE_avg", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("i_dE_squared_average", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("i_marker_den", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("e_dr_avg", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("e_dr_squared_average", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("e_dE_avg", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("e_dE_squared_average", {ntri}, {0}, {ntri});
        this->output_io.DefineVariable<double>("e_marker_den", {ntri}, {0}, {ntri});

        this->writer = output_io.Open("xgc.diffusion.bp", adios2::Mode::Write, MPI_COMM_SELF);

        first = false;
    }

    this->writer.BeginStep();
    this->writer.Put<double>("i_dr_avg", this->i_dr_avg.data());
    this->writer.Put<double>("i_dr_squared_average", this->i_dr_squared_average.data());
    this->writer.Put<double>("i_dE_avg", this->i_dE_avg.data());
    this->writer.Put<double>("i_dE_squared_average", this->i_dE_squared_average.data());
    this->writer.Put<double>("i_marker_den", this->i_marker_den.data());
    this->writer.Put<double>("e_dr_avg", this->e_dr_avg.data());
    this->writer.Put<double>("e_dr_squared_average", this->e_dr_squared_average.data());
    this->writer.Put<double>("e_dE_avg", this->e_dE_avg.data());
    this->writer.Put<double>("e_dE_squared_average", this->e_dE_squared_average.data());
    this->writer.Put<double>("e_marker_den", this->e_marker_den.data());
    this->writer.EndStep();
}