#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <cmath>
using namespace std;


// ==================================================================================
/* THE BOLTZMANN FACTOR FUNCTION */
// ==================================================================================
double get_boltzmann_factor(System &system, double e_i, double e_f, string movetype) {
    double bf, MAXVALUE=1e10; // we won't care about huge bf's

    if (system.constants.ensemble == "uvt") {
            
    }
    else if (system.constants.ensemble == "nvt") {

    }
    else if (system.constants.ensemble == "npt") {
        if (movetype == "volume") {
            // Frenkel Smit p118
            bf= exp(-( (e_f - e_i)
            + system.constants.pres * system.constants.ATM2REDUCED * 
            (system.pbc.volume - system.pbc.old_volume)
            - (system.stats.count_movables + 1) * system.constants.temp * 
            log(system.pbc.volume/system.pbc.old_volume))/system.constants.temp);
            
            system.stats.volume_change_bf_sum += bf;
            //printf("old_vol: %f; vol: %f\n", system.pbc.old_volume, system.pbc.volume);
        }
    }
    else if (system.constants.ensemble == "nve") {

    }
  
    if (std::isinf(bf)) {
        bf = MAXVALUE; 
        //printf("bf is INFINITE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111\n");
    }
    //printf("bf: %f\n",bf);
    return bf;
}
