#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <cmath>
using namespace std;


// ==================================================================================
/* THE BOLTZMANN FACTOR FUNCTION */
// ==================================================================================
double get_boltzmann_factor(System &system, double e_i, double e_f, int_fast8_t movetype) {
    double bf, MAXVALUE=1e4; // we won't care about huge bf's
    double energy_delta = e_f - e_i;
    double fugacity;

    if (system.proto.size() == 1) fugacity = system.constants.pres;
    else fugacity = system.proto[system.constants.currentprotoid].fugacity;

    if (system.constants.ensemble == ENSEMBLE_UVT) {
        if (movetype == MOVETYPE_INSERT) {
            bf = system.pbc.volume * fugacity * 
            system.constants.ATM2REDUCED/(system.constants.temp * 
            (double)(system.stats.count_movables)) *
                exp(-energy_delta/system.constants.temp);
            if (bf < MAXVALUE) system.stats.insert_bf_sum += bf;
            else system.stats.insert_bf_sum += MAXVALUE;
        } 
        else if (movetype == MOVETYPE_REMOVE) {
            bf = system.constants.temp * 
            ((double)(system.stats.count_movables) + 1.0)/
            (system.pbc.volume* fugacity *system.constants.ATM2REDUCED) *
                exp(-energy_delta/system.constants.temp);
            if (bf < MAXVALUE) system.stats.remove_bf_sum += bf;
            else system.stats.remove_bf_sum += MAXVALUE;
        }
        else if (movetype == MOVETYPE_DISPLACE) {
            bf = exp(-energy_delta/(system.constants.temp));
            if (bf < MAXVALUE) system.stats.displace_bf_sum += bf;
            else system.stats.displace_bf_sum += MAXVALUE;
        }
    }
    else if (system.constants.ensemble == ENSEMBLE_NVT) {
        if (movetype == MOVETYPE_DISPLACE) {
            bf = exp(-energy_delta/(system.constants.temp));
            if (bf < MAXVALUE) system.stats.displace_bf_sum += bf;
            else system.stats.displace_bf_sum += MAXVALUE;
        }
    }
    else if (system.constants.ensemble == ENSEMBLE_NPT) {
        if (movetype == MOVETYPE_VOLUME) {
            // Frenkel Smit p118
            bf= exp(-( (energy_delta)
            + system.constants.pres * system.constants.ATM2REDUCED * 
            (system.pbc.volume - system.pbc.old_volume)
            - (system.stats.count_movables + 1) * system.constants.temp * 
                log(system.pbc.volume/system.pbc.old_volume))/system.constants.temp);
            if (bf < MAXVALUE) system.stats.volume_change_bf_sum += bf;
            else system.stats.volume_change_bf_sum += MAXVALUE;
        }
        else if (movetype == MOVETYPE_DISPLACE) {
            bf = exp(-energy_delta/system.constants.temp);
            if (bf < MAXVALUE) system.stats.displace_bf_sum += bf;
            else system.stats.displace_bf_sum += MAXVALUE;
        }
    }
    else if (system.constants.ensemble == ENSEMBLE_NVE) {
        if (movetype == MOVETYPE_DISPLACE) {
            double exponent = 3.0*system.stats.count_movables/2.0;
            //printf("exponent = %f\n", exponent);

            //printf("num %f\n", pow( (system.constants.total_energy - e_f ), 3.0*system.stats.count_movables/2.0));
            //printf("denom %f\n", pow( (system.constants.total_energy - e_i), 3.0*system.stats.count_movables/2.0));
            bf = pow(
            (system.constants.total_energy - e_f) , exponent) / 
                pow(
                    (system.constants.total_energy - e_i) , exponent);
            if (bf < MAXVALUE) system.stats.displace_bf_sum += bf;
            else system.stats.displace_bf_sum += MAXVALUE;
            //printf("bf %f \n", bf);
        }
    }
  
    if (std::isinf(bf)) {
        bf = MAXVALUE; 
        //printf("bf is INFINITE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111\n");
    }
    //printf("bf: %f\n",bf);
    return bf;
}
