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
    double bf, MAXVALUE=1e4; // we won't care about huge bf's for averaging
    double energy_delta = e_f - e_i;
    double fugacity;

    // determine fugacity
    if (system.proto.size() == 1 && system.constants.sorbate_name.size() == 0 && system.constants.fugacity_single == 0) fugacity = system.constants.pres;
    else if (system.constants.fugacity_single == 1) fugacity = system.proto[0].fugacity;
    else fugacity = system.proto[system.constants.currentprotoid].fugacity;

    // check for N bias (loading bias)
    if (system.constants.bias_uptake_switcher) {
        if ((double)system.stats.count_movables >= system.constants.bias_uptake)
            //fugacity *= 3; // big, but not too big, for biased insertions
            system.constants.bias_uptake_switcher=0; // done biasing if we exceed the bias-N
            printf("Loading bias deactivated! (Desired N reached).\n");
        }       
    }
      //printf("fugac: %f\n", fugacity);

    if (system.constants.ensemble == ENSEMBLE_UVT) {
        if (movetype == MOVETYPE_INSERT) {
            bf = system.pbc.volume * fugacity * 
            system.constants.ATM2REDUCED/(system.constants.temp * 
            (double)(system.stats.count_movables)) *
                exp(-energy_delta/system.constants.temp) *
                (double)system.proto.size(); // bias for multisorbate (thus no change for single)
            if (bf < MAXVALUE) system.stats.insert_bf_sum += bf;
            else system.stats.insert_bf_sum += MAXVALUE;
        } 
        else if (movetype == MOVETYPE_REMOVE) {
            bf = system.constants.temp * 
            ((double)(system.stats.count_movables) + 1.0)/
            (system.pbc.volume* fugacity *system.constants.ATM2REDUCED) *
                exp(-energy_delta/system.constants.temp) /
                (double)system.proto.size(); // bias for multisorbate (thus no change for single)
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
            /*
            system.pbc.old_volume = 200000;
            system.pbc.volume = 210000;
            system.constants.pres = 1.;
            energy_delta = -13.5;
            system.constants.temp = 300;
            system.stats.count_movables = 50;
*/
            // Frenkel Smit p118
            bf= exp(-( (energy_delta)
            + system.constants.pres * system.constants.ATM2REDUCED * 
            (system.pbc.volume - system.pbc.old_volume)
            - (system.stats.count_movables + 1) * system.constants.temp * 
                log(system.pbc.volume/system.pbc.old_volume))/system.constants.temp);

//            printf("\n\nvolume boltz: %f; v_old %f; v_new %f; vdiff %f; Udiff %f\n\n", bf, system.pbc.old_volume, system.pbc.volume, system.pbc.volume - system.pbc.old_volume, energy_delta);
            
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
 
    if (bf > MAXVALUE) {
        bf = MAXVALUE;
    } 
    else if (std::isinf(bf)) {
        printf("GOT INF! bf = %f\n", bf);
        //bf = MAXVALUE; 
    }
    //printf("bf: %f initial %f final %f\n",bf, e_i, e_f);
    return bf;
}
