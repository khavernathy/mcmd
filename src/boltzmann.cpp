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
    double energy_delta = e_f - e_i;

    if (system.constants.ensemble == "uvt") {
        if (movetype == "add") {
            bf = system.pbc.volume * system.constants.pres * 
            system.constants.ATM2REDUCED/(system.constants.temp * 
            (double)(system.stats.count_movables)) *
                exp(-energy_delta/system.constants.temp);
            system.stats.insert_bf_sum += bf;
        } 
        else if (movetype == "remove") {
            bf = system.constants.temp * 
            ((double)(system.stats.count_movables) + 1.0)/
            (system.pbc.volume*system.constants.pres*system.constants.ATM2REDUCED) *
                exp(-energy_delta/system.constants.temp);
            system.stats.remove_bf_sum += bf;
        }
        else if (movetype == "displace") {
            bf = exp(-energy_delta/(system.constants.temp));
            system.stats.displace_bf_sum += bf;
        }
    }
    else if (system.constants.ensemble == "nvt") {
        if (movetype == "displace") {
            bf = exp(-energy_delta/(system.constants.temp));
            system.stats.displace_bf_sum += bf;
        }
    }
    else if (system.constants.ensemble == "npt") {
        if (movetype == "volume") {
            // Frenkel Smit p118
            bf= exp(-( (energy_delta)
            + system.constants.pres * system.constants.ATM2REDUCED * 
            (system.pbc.volume - system.pbc.old_volume)
            - (system.stats.count_movables + 1) * system.constants.temp * 
                log(system.pbc.volume/system.pbc.old_volume))/system.constants.temp);
            system.stats.volume_change_bf_sum += bf;
        }
        else if (movetype == "displace") {
            bf = exp(-energy_delta/system.constants.temp);
            system.stats.displace_bf_sum += bf;
        }
    }
    else if (system.constants.ensemble == "nve") {
        if (movetype == "displace") {
            bf = pow((system.constants.total_energy - e_f),3.0*system.stats.count_movables/2.0) / pow((system.constants.total_energy - e_i),(3.0*system.stats.count_movables/2.0));
            system.stats.displace_bf_sum += bf;
        }
    }
  
    if (std::isinf(bf)) {
        bf = MAXVALUE; 
        //printf("bf is INFINITE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111\n");
    }
    //printf("bf: %f\n",bf);
    return bf;
}
