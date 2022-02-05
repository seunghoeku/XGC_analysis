#include <vector>

#include "particles.hpp"
#include "flags.hpp"
#include "sml.hpp"
#include "heatload.hpp"

extern Simulation sml;

// get heatload of single species
void heatload_calc(const Particles &div, HeatLoad &sp) {
    
    for(int i=0; i<div.size(); i++) {

        struct Particle p = div[i];
        double en = sml.c2_2m * p.rho * p.rho * p.B * p.B + p.mu*p.B;
        double wp = p.dw * p.w0;

        Conditions cond(p);

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
