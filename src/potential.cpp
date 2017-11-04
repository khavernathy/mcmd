#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

//#include "distance.cpp"
#include "mixing.cpp"
#include "lj.cpp"
#include "commy.cpp"
#include "coulombic.cpp"
#include "polar.cpp"
#include "tt.cpp"
#include "pairs.cpp"

// =================== MAIN FUNCTION ======================
// ---------------POTENTIAL OF ENTIRE SYSTEM --------------

double getTotalPotential(System &system) {
    system.checkpoint("started getTotalPotential()");
    int_fast8_t model = system.constants.potential_form;

    // initializers
    double total_potential=0;
    double total_rd=0.0; double total_es = 0.0; double total_polar=0.0;
    system.constants.auto_reject=0;

// =========================================================================
if (system.molecules.size() > 0) { // don't bother with 0 molecules!

    // REPULSION DISPERSION.
    if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJPOLAR || model == POTENTIAL_LJESPOLAR) {
        total_rd = lj(system);
    } else if (model == POTENTIAL_COMMY || model == POTENTIAL_COMMYES || model == POTENTIAL_COMMYESPOLAR) {
        total_rd = commy(system);
    } else if (model == POTENTIAL_TT || model == POTENTIAL_TTES || model == POTENTIAL_TTESPOLAR) {
        total_rd = tt(system);
    }
    
    // ELECTROSTATIC
    if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_COMMYES || model == POTENTIAL_COMMYESPOLAR || model == POTENTIAL_TTES || model == POTENTIAL_TTESPOLAR) {

        if (system.constants.mode=="md" || (!system.constants.auto_reject_option || !system.constants.auto_reject)) { // these only run if MD, or if no bad contact was discovered in MC
            if (system.constants.ewald_es)
                total_es = coulombic_ewald(system); // using ewald method for es
            else
                total_es = coulombic(system); // plain old coloumb
        }
    }
    
    // POLARIZATION
    if (model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR || model == POTENTIAL_COMMYESPOLAR || model == POTENTIAL_TTESPOLAR) {

        if (system.constants.mode=="md" || (!system.constants.auto_reject_option || !system.constants.auto_reject)) { // these only run if MD, or no bad contact was discovered in MC

                total_polar = polarization(system);
        }

    }
}
// ==========================================================================

    total_potential = total_rd + total_es + total_polar;

    // save values to system vars
    system.stats.rd.value = total_rd;
    system.stats.es.value = total_es;
    system.stats.polar.value = total_polar;
    system.stats.potential.value = total_potential;

    system.checkpoint("ended getTotalPotential()");
	return total_potential;
}
