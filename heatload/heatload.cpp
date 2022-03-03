#include "heatload.hpp"

void init(adios2::ADIOS *ad);                                             // initialization
void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db); // calculate heatload

Simulation sml; // input parameters that controls simulation.
t_ParticleDB iesc_db;
t_ParticleDB eesc_db;
MPI_Comm heatload_comm;
int heatload_comm_size;
int heatload_comm_rank;

void heatload_init(adios2::ADIOS *ad, MPI_Comm comm)
{
    // init simulation parameters
    init(ad);

    // init adios
    load_init(ad, "xgc.escaped_ptls.bp", comm);

    heatload_comm = comm;
    MPI_Comm_rank(heatload_comm, &heatload_comm_rank);
    MPI_Comm_size(heatload_comm, &heatload_comm_size);
}

void heatload_init2(adios2::ADIOS *ad)
{
    // init simulation parameters
    init(ad);

    // init adios
    // load_init(ad, "xgc.escaped_ptls.bp");
}

int heatload_step(adios2::ADIOS *ad, int istep)
{
    Particles idiv;
    Particles ediv;
    t_ParticlesList iesc;
    t_ParticlesList eesc;

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

    if (heatload_comm_rank == 0)
    {
        std::cout << std::endl;
        std::cout << ">>> Step: " << istep << std::endl;
        std::cout << "Num. of escaped ions: " << iesc.size() << std::endl;
        std::cout << "Num. of escaped elec: " << eesc.size() << std::endl;
        std::cout << "Num. of divertor ions: " << idiv.size() << std::endl;
        std::cout << "Num. of divertor elec: " << ediv.size() << std::endl;

        // print first 10 esc particles
        int count = 0;
        t_ParticlesList::iterator it;
        for (it = iesc.begin(); it != iesc.end(); it++)
        {
            printf("iesc gid, rzphi, flag: %lld %f %f %f %d\n", it->second.gid, it->second.r, it->second.z,
                   it->second.phi, it->second.flag);
            count++;
            if (count > 10)
                break;
        }

        // separate divertor particles and escaped particles
        iesc_db.push_back(iesc);
        eesc_db.push_back(eesc);
        Particle ptl = search(iesc_db, istep - 1, 15824414);
        printf("Found or not? gid=%lld\n", ptl.gid);

        // store escaped particles to DB

        // Calculate heatload from divertor particles
        HeatLoad ion(1);
        HeatLoad elec(0);

        heatload_calc(idiv, ion, iesc_db); // need to send DB
        heatload_calc(ediv, elec, eesc_db);
        output(ad, ion, elec);
    }

    return 0;
}

void heatload_finalize(adios2::ADIOS *ad)
{
    load_finalize();
    output_finalize();
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

        std::cout << std::endl;
        std::cout << ">>> Step: " << i << std::endl;
        std::cout << "Num. of escaped ions: " << iesc.size() << std::endl;
        std::cout << "Num. of escaped elec: " << eesc.size() << std::endl;
        std::cout << "Num. of divertor ions: " << idiv.size() << std::endl;
        std::cout << "Num. of divertor elec: " << ediv.size() << std::endl;

        // print first 10 esc particles
        int count = 0;
        t_ParticlesList::iterator it;

        for (it = iesc.begin(); it != iesc.end(); it++)
        {
            printf("iesc gid, rzphi, flag: %lld %f %f %f %d\n", it->second.gid, it->second.r, it->second.z,
                   it->second.phi, it->second.flag);
            count++;
            if (count > 10)
                break;
        }

        // separate divertor particles and escaped particles
        iesc_db.push_back(iesc);
        eesc_db.push_back(eesc);
        Particle ptl = search(iesc_db, i - 1, 15824414);
        printf("Found or not? gid=%lld\n", ptl.gid);

        // store escaped particles to DB

        // Calculate heatload from divertor particles
        HeatLoad ion(1);
        HeatLoad elec(0);

        heatload_calc(idiv, ion, iesc_db); // need to send DB
        heatload_calc(ediv, elec, eesc_db);
        output(ad, ion, elec);
    }

    load_finalize();
    output_finalize();
}
