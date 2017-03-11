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

// =================== MAIN FUNCTION ======================
// ---------------POTENTIAL OF ENTIRE SYSTEM --------------
double * getTotalPotential(System &system, string model) {

    // initializers 
    double total_potential=0;    	
    double total_rd=0.0; double total_es = 0.0; double total_polar=0.0;


// =========================================================================
    // REPULSION DISPERSION.
    if (model == "lj" || model == "ljes" || model == "ljespolar") {
        total_rd = lj(system);
    }
    // ELECTROSTATIC
    if (model == "ljes" || model == "ljespolar") {
        if (system.constants.ewald_es == "on") 
            total_es = coulombic_ewald(system); // using ewald method for es
        else 
            total_es = coulombic(system); // plain old coloumb
    }
    // POLARIZATION 
    if (model == "ljespolar") {
        total_polar = polarization(system); // yikes
    } 
// ==========================================================================


    // save values to stats vars
    total_potential = total_rd + total_es + total_polar;
    
    system.stats.rdU = total_rd; system.stats.current_rd_sum += total_rd;
    system.stats.esU = total_es; system.stats.current_es_sum += total_es;
    system.stats.polarU = total_polar; system.stats.current_polar_sum += total_polar;
    system.stats.totalU += total_potential; system.stats.current_energy_sum += total_potential;


	static double output[4];
        output[0] = total_rd;
        output[1] = total_es;
        output[2] = total_polar;
	return output;
}
