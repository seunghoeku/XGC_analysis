#include <ctime>
#include <iostream>
#include <stdio.h>
#include <vector>

#include "flags.hpp"
#include "heatload_calc.hpp"
#include "particles.hpp"
#include "sml.hpp"
#ifdef USE_OMP
#include <omp.h>
#endif

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup.hpp>

#define LOG BOOST_LOG_TRIVIAL(debug)

extern Simulation sml;
extern Particle search(t_ParticleDB &db, int timestep, long long gid);

static long int _progress_step_current = 0;
void progress_step(long int total)
{
    ++_progress_step_current;
    // if (_progress_step_current % 1'000 == 0)
    //     std::cerr << ".";
    // if ((_progress_step_current - 1) % 50'000 == 0)
    //     fprintf(stderr, "\n%ld/%ld ", _progress_step_current - 1, total);
    if ((_progress_step_current) % 50'000 == 0)
        LOG << _progress_step_current << "/" << total << " processed";
}

// get heatload of single species
void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db, int show_progress = 0)
{

    printf("\nHeatload calc particle size: %ld\n", div.size());
    // reset progress bar
    _progress_step_current = 0;
    std::time_t start = std::time(nullptr);
#pragma omp parallel for default(none) shared(sml, div, db, sp, show_progress, std::cerr)
    for (int i = 0; i < div.size(); i++)
    {
        if (show_progress)
        {
#pragma omp critical
            {
                progress_step(div.size());
            }
        }

        struct Particle p = div[i]; // particle that hit divertor
        double en = sml.c2_2m[sp.isp] * p.rho * p.rho * p.B * p.B + p.mu * p.B;
        double wp = p.dw * p.w0;

        struct Particle p_esc = search(db, p.esc_step, p.gid); // particle info when it escaped.
        // printf("%lld %d\n", p.gid, p.esc_step);
        if (p_esc.gid > 0)
        {
            Conditions cond(p_esc); // get conditions from particle info when escaped.

            // check inner or outer
            for (int side = 0; side < 2; side++)
            {
                // check bounding box
                if (sml.rmin[side] < p.r && p.r < sml.rmax[side] && sml.zmin[side] < p.z && p.z < sml.zmax[side])
                {

                    double pn = (p.psi - sml.pmin[side]) / sml.dpsi[side];
                    int ip = (int)(pn);
                    if (ip >= 0 && ip < sml.npsi - 1)
                    {
                        double ws = 1.0 - pn + (double)(ip);

                        for (int icond = 0; icond < N_COND; icond++)
                        {
                            if (cond.b[icond])
                            {
#pragma omp critical(spupdate)
                                {
#ifdef USE_OMP
                                    printf("%d: thread rank %d\n", i, omp_get_thread_num());
#endif
                                    sp.side[side].ptl[icond][ip] = wp * ws;
                                    sp.side[side].ptl[icond][ip + 1] = wp * ws;

                                    sp.side[side].en[icond][ip] = en * wp * ws;
                                    sp.side[side].en[icond][ip + 1] = en * wp * ws;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("\nWall time (seconds): %f\n", std::difftime(std::time(nullptr), start));
}
