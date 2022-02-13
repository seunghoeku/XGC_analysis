#include <vector>
#include <iostream>
#include <stdio.h>

#include "particles.hpp"
#include "flags.hpp"
#include "sml.hpp"
#include "heatload.hpp"
#ifdef USE_OMP
#include <omp.h>
#endif

extern Simulation sml;
extern Particle search(t_ParticleDB &db, int timestep, long long gid);

// get heatload of single species
void heatload_calc(const Particles &div, HeatLoad &sp, t_ParticleDB &db) {
    
    printf ("\nHeatload calc particle size: %d\n", div.size());
    #pragma omp parallel for default(none) shared(sml, div, db, sp, std::cerr)
    for(int i=0; i<div.size(); i++) {
        // printf("%d: thread rank %d\n", i, omp_get_thread_num());
        if ((i+1)%100==0) std::cerr << ".";
        if ((i+1)%5000==0) std::cerr << std::endl << i+1;

        struct Particle p = div[i]; // particle that hit divertor
        double en = sml.c2_2m[sp.isp] * p.rho * p.rho * p.B * p.B + p.mu*p.B;
        double wp = p.dw * p.w0;

        struct Particle p_esc = search(db, p.esc_step, p.gid); // particle info when it escaped.
        Conditions cond(p_esc); // get conditions from particle info when escaped.

        //check inner or outer
        for(int side=0; side <2; side++) {
            //check bounding box
            if( sml.rmin[side] < p.r && p.r < sml.rmax[side] &&  sml.zmin[side] < p.z && p.z < sml.zmax[side]) {
                
                double pn=(p.psi-sml.pmin[side])/sml.dpsi[side];
                int ip=(int)(pn);
                if(ip>=0 && ip<sml.npsi-1){
                    double ws=1.0 - pn + (double)(ip);

                    for(int icond=0; icond<N_COND; icond++){
                        if(cond.b[icond]){
#ifdef USE_OMP
                            #pragma omp critical(spupdate)
                            {
                                printf("%d: thread rank %d\n", i, omp_get_thread_num());
#else
                            {
#endif
                                sp.side[side].ptl[icond][ip]   =  wp * ws;
                                sp.side[side].ptl[icond][ip+1] =  wp * ws;

                                sp.side[side].en[icond][ip]   = en * wp * ws;
                                sp.side[side].en[icond][ip+1] = en * wp * ws;
                            }
                        }
                    }
                }

            }
        }

    }
}
