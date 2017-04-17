#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

#include <distance.cpp>
#include <lj.cpp>
#include <coulombic.cpp>
#include <polar.cpp>
#include <pairs.cpp>

// =================== MAIN FUNCTION ======================
// ---------------POTENTIAL OF ENTIRE SYSTEM --------------
double getTotalPotential(System &system) {
    int_fast8_t model = system.constants.potential_form;
    // compute all interaction distances
    //make_pairs(system);

    // initializers
    double total_potential=0;
    double total_rd=0.0; double total_es = 0.0; double total_polar=0.0;


// =========================================================================
    // REPULSION DISPERSION.
    if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR) {
        total_rd = lj(system);
    }
    // ELECTROSTATIC
    if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR) {
        if (system.constants.ewald_es)
            total_es = coulombic_ewald(system); // using ewald method for es
        else
            total_es = coulombic(system); // plain old coloumb
    }
    // POLARIZATION
    if (model == POTENTIAL_LJESPOLAR) {
        total_polar = polarization(system); // yikes
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
