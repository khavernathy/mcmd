#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

//#include "distance.cpp"
#include "lj.cpp"
#include "commy.cpp"
#include "coulombic.cpp"
#include "polar.cpp"
#include "pairs.cpp"

// =================== MAIN FUNCTION ======================
// ---------------POTENTIAL OF ENTIRE SYSTEM --------------
double getTotalPotential(System &system) {
    int_fast8_t model = system.constants.potential_form;
    // compute all interaction distances
    //make_pairs(system);

    // initializers
    double total_potential=0;
    double total_rd=0.0; double total_es = 0.0; double total_polar=0.0;
    system.constants.auto_reject=0;

// =========================================================================
if (system.molecules.size() > 0) { // don't bother with 0 molecules!
    system.checkpoint("starting RD energy calculation.");
    // REPULSION DISPERSION.
    if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJPOLAR || model == POTENTIAL_LJESPOLAR) {
        total_rd = lj(system);
    } else if (model == POTENTIAL_COMMY || model == POTENTIAL_COMMYES || model == POTENTIAL_COMMYESPOLAR) {
        total_rd = commy(system);
    }
    if (system.constants.mode=="md" || (!system.constants.auto_reject_option || !system.constants.auto_reject)) { // these only run if no bad contact was discovered in MC
    system.checkpoint("starting ES energy calculation");
    // ELECTROSTATIC
    if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_COMMYES || model == POTENTIAL_COMMYESPOLAR) {
        if (system.constants.ewald_es)
            total_es = coulombic_ewald(system); // using ewald method for es
        else
            total_es = coulombic(system); // plain old coloumb
    }
    system.checkpoint("starting Polarization energy calculation.");
    // POLARIZATION
    if (model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR || model == POTENTIAL_COMMYESPOLAR) {
        total_polar = polarization(system); // yikes
    }
    }
}
// ==========================================================================

    total_potential = total_rd + total_es + total_polar;

    // save values to vars
    system.stats.rd.value = total_rd;
    system.stats.es.value = total_es;
    system.stats.polar.value = total_polar;
    system.stats.potential.value = total_potential;

//    printf("MC STEP %i ::: rd %f es %f pol %f tot %f\n", system.stats.MCstep, total_rd, total_es, total_polar, total_potential);
	return total_potential;
}
