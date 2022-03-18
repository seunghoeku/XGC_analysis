#include "heatload.hpp"
#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

void init(adios2::ADIOS *ad, std::string xgcdir); // initialization
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
    // init simulation parameters
    init(ad, xgcdir);

    // init adios
    boost::filesystem::path fname = boost::filesystem::path(xgcdir) / boost::filesystem::path("xgc.escaped_ptls.bp");
    load_init(ad, fname.string(), comm);

    heatload_comm = comm;
    MPI_Comm_rank(heatload_comm, &heatload_comm_rank);
    MPI_Comm_size(heatload_comm, &heatload_comm_size);

    if (read_restart)
    {
        ptldb_load(iesc_db, "heatload_iesc_db.bp");
        ptldb_load(eesc_db, "heatload_eesc_db.bp");
    }
}

void heatload_init2(adios2::ADIOS *ad, std::string xgcdir)
{
    // init simulation parameters
    init(ad, xgcdir);

    // init adios
    // load_init(ad, "xgc.escaped_ptls.bp");
}

void test()
{
    LOG << "heatload_step";
    t_ParticlesList pmap;
    for (int i = 0; i < 10; i++)
    {
        Particle ptl;
        ptl.gid = i * 10;
        add(pmap, ptl);
    }
    ptlmap_print(pmap, "pmap");

    t_ParticleDB db;
    db.push_back(pmap);

    for (int i = 0; i < 50; i++)
    {
        Particle p = search(db, 0, i);
        LOG << i << " " << p.gid;
    }
}

int heatload_step(adios2::ADIOS *ad, int istep, bool ion_only)
{
    Particles idiv;
    Particles ediv;
    t_ParticlesList iesc;
    t_ParticlesList eesc;

    // idiv, ediv (local), iesc, eesc (global)
    adios2::StepStatus status = load_data(idiv, ediv, iesc, eesc);
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

    // Sync iesc and iesc with rank 0
    ptlmap_sync(iesc, heatload_comm);
    ptlmap_sync(eesc, heatload_comm);

    LOG << ">>> Step: " << istep;
    LOG << "Num. of escaped ions: " << ptlmap_count(iesc);
    LOG << "Num. of escaped elec: " << ptlmap_count(eesc);
    LOG << "Num. of divertor ions: " << idiv.size();
    LOG << "Num. of divertor elec: " << ediv.size();

    // separate divertor particles and escaped particles
    iesc_db.push_back(iesc);
    eesc_db.push_back(eesc);
    ptldb_print(iesc_db, "iesc_db");
    ptldb_print(eesc_db, "eesc_db");

    // store escaped particles to DB

    // Calculate heatload from divertor particles
    HeatLoad ion(1);
    HeatLoad elec(0);

    heatload_calc(idiv, ion, iesc_db); // need to send DB
    if (!ion_only)
    {
        heatload_calc(ediv, elec, eesc_db);
    }
    output(ad, ion, elec, heatload_comm);

    return 0;
}

void heatload_finalize()
{
    load_finalize();
    output_finalize(heatload_comm);

    if (heatload_comm_rank == 0)
    {
        ptldb_save(iesc_db, "heatload_iesc_db.bp");
        ptldb_save(eesc_db, "heatload_eesc_db.bp");
    }
}

void heatload(adios2::ADIOS *ad)
{
    int i = 0;
    while (1)
    {
        i++;

        Particles idiv;
        Particles ediv;
        t_ParticlesList iesc;
        t_ParticlesList eesc;

        // idiv, ediv (local), iesc, eesc (global)
        adios2::StepStatus status = load_data(idiv, ediv, iesc, eesc);
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
        iesc_db.push_back(iesc);
        eesc_db.push_back(eesc);

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
