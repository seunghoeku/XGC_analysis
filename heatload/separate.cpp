#include <vector>

#include "particles.hpp"
#include "flags.hpp"
#include "sml.hpp"

using namespace std;
extern Simulation sml;

void separate(    std::vector<long> gid, std::vector<int> flag, std::vector<float> dw, 
    std::vector<float> phase,  std::vector<Particles> div,  std::vector<Particles> esc) 
{

    int n; // size of input vector

    n = gid.size(); // size of gid should be the same with flag, dw, phase

    for(int i=0; i<n; i++){
        // check flag 
        int flag1 ; // tmp flag     
        flag1 = flag[i];

        Flags fl(flag1); // decode flags

        // save to div or esc
        if(fl.escaped) {
            // add to esc

        } else {
            // add to div
        }
    }


}