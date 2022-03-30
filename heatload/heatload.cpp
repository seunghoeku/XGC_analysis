#include "heatload.hpp"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#include "cam_timers.hpp"

void init(adios2::ADIOS *ad, std::string xgcdir);                         // initialization
void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db); // calculate heatload

#define LOG BOOST_LOG_TRIVIAL(debug)

Simulation sml; // input parameters that controls simulation.
t_ParticleDB iesc_db;
t_ParticleDB eesc_db;
MPI_Comm heatload_comm;
int heatload_comm_size;
int heatload_comm_rank;

void heatload_init(adios2::ADIOS *ad, MPI_Comm comm, std::string xgcdir, bool read_restart)
{
    TIMER_START("INIT");
    // init adios
    boost::filesystem::path fname = boost::filesystem::path(xgcdir) / boost::filesystem::path("xgc.escaped_ptls.bp");
    load_init(ad, fname.string(), comm);

    // init simulation parameters
    init(ad, xgcdir);

    heatload_comm = comm;
    MPI_Comm_rank(heatload_comm, &heatload_comm_rank);
    MPI_Comm_size(heatload_comm, &heatload_comm_size);

    if (read_restart)
    {
        ptldb_load(iesc_db, "heatload_iesc_db.bp", heatload_comm);
        ptldb_load(eesc_db, "heatload_eesc_db.bp", heatload_comm);
    }
    TIMER_STOP("INIT");
}

void heatload_init2(adios2::ADIOS *ad, std::string xgcdir)
{
    // init simulation parameters
    init(ad, xgcdir);

    // init adios
    // load_init(ad, "xgc.escaped_ptls.bp");
}

int heatload_step(adios2::ADIOS *ad, int istep, bool ion_only)
{
    TIMER_START("STEP");
    Particles idiv;
    Particles ediv;
    t_ParticlesList iesc;
    t_ParticlesList eesc;

    // idiv, ediv (local), iesc, eesc (local)
    int timestep;
    adios2::StepStatus status = load_data(idiv, ediv, iesc, eesc, timestep);
    LOG << "Done with load_data: xgc timestep = " << timestep;
    if (status == adios2::StepStatus::EndOfStream)
    {
        std::cout << "Input stream terminated. Exit loop" << std::endl;
        return -1;
    }
    else if (status == adios2::StepStatus::NotReady)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        return 1;
    }
    else if (status == adios2::StepStatus::OtherError)
    {
        std::cout << "Input stream had errors. Exit loop" << std::endl;
        return -1;
    }

    LOG << ">>> Step: " << istep;
    LOG << "Num. of escaped ions: " << ptlmap_count(iesc);
    LOG << "Num. of escaped elec: " << ptlmap_count(eesc);

    // store escaped particles to DB
    insert_or_append(iesc_db, timestep, iesc);
    insert_or_append(eesc_db, timestep, eesc);
    ptldb_print(iesc_db, "iesc_db");
    ptldb_print(eesc_db, "eesc_db");

    // Get how many div particles each one has
    int idiv_len = idiv.size();
    int ediv_len = ediv.size();
    std::vector<int> idiv_len_list(heatload_comm_size);
    std::vector<int> ediv_len_list(heatload_comm_size);
    MPI_Allgather(&idiv_len, 1, MPI_INT, idiv_len_list.data(), 1, MPI_INT, heatload_comm);
    MPI_Allgather(&ediv_len, 1, MPI_INT, ediv_len_list.data(), 1, MPI_INT, heatload_comm);

    // Calculate heatload from divertor particles
    HeatLoad ion(1);
    HeatLoad elec(0);

    // div from the next neighbor
    Particles idiv2;
    Particles ediv2;
    Particles &current_idiv = idiv;
    Particles &current_ediv = ediv;
    Particles &next_idiv = idiv2;
    Particles &next_ediv = ediv2;

    int iround = 0;
    do
    {
        LOG << ">>> Step,round: " << istep << " " << iround;
        LOG << "Num. of divertor ions: " << current_idiv.size();
        LOG << "Num. of divertor elec: " << current_ediv.size();

        heatload_calc(current_idiv, ion, iesc_db); // need to send DB
        // Debug
        // ptldb_dump(iesc_db, "iesc_db");
        // auto p = search(iesc_db, 0, 820040877);
        // LOG << "Found? " << p.gid;
        if (!ion_only)
        {
            heatload_calc(current_ediv, elec, eesc_db);
        }

        iround++;

        // Get div from the next neighbor
        // No need to do in the last round
        if (iround < heatload_comm_size)
        {
            next_idiv.clear();
            next_ediv.clear();
            ptls_shift(current_idiv, next_idiv, heatload_comm);
            ptls_shift(current_ediv, next_ediv, heatload_comm);

            // swap reference
            current_idiv = &current_idiv == &idiv ? idiv2 : idiv;
            current_ediv = &current_ediv == &ediv ? ediv2 : ediv;
            next_idiv = &next_idiv == &idiv2 ? idiv : idiv2;
            next_ediv = &next_ediv == &ediv2 ? ediv : ediv2;
        }
    } while (iround < heatload_comm_size);

    output(ad, ion, elec, heatload_comm);

#ifdef CAM_TIMERS
    GPTLprint_memusage("STEP MEMUSAGE");
#endif
    TIMER_STOP("STEP");
    return 0;
}

void heatload_finalize()
{
    TIMER_START("FINALIZE");
    load_finalize();
    output_finalize(heatload_comm);

    ptldb_save(iesc_db, "heatload_iesc_db.bp", heatload_comm);
    ptldb_save(eesc_db, "heatload_eesc_db.bp", heatload_comm);
    TIMER_STOP("FINALIZE");
}

void heatload(adios2::ADIOS *ad)
{
    int i = 0;
    while (true)
    {
        i++;

        Particles idiv;
        Particles ediv;
        t_ParticlesList iesc;
        t_ParticlesList eesc;

        // idiv, ediv (local), iesc, eesc (global)
        int timestep;
        adios2::StepStatus status = load_data(idiv, ediv, iesc, eesc, timestep);
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

        // Sync iesc and iesc with rank 0
        ptlmap_sync(iesc, heatload_comm);
        ptlmap_sync(eesc, heatload_comm);

        LOG << ">>> Step: " << i;
        LOG << "Num. of escaped ions: " << iesc.size();
        LOG << "Num. of escaped elec: " << eesc.size();
        LOG << "Num. of divertor ions: " << idiv.size();
        LOG << "Num. of divertor elec: " << ediv.size();

        // separate divertor particles and escaped particles
        insert_or_append(iesc_db, timestep, iesc);
        insert_or_append(eesc_db, timestep, eesc);

        // store escaped particles to DB

        // Calculate heatload from divertor particles
        HeatLoad ion(1);
        HeatLoad elec(0);

        heatload_calc(idiv, ion, iesc_db); // need to send DB
        heatload_calc(ediv, elec, eesc_db);
        output(ad, ion, elec, heatload_comm);
    }

    load_finalize();
    output_finalize(heatload_comm);
}
