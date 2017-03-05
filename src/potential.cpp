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
	//long unsigned int size = system.atoms.size();
    //printf("STARTING POTENTIAL FUNCTION AGAIN.\n\n");

    // initializers     	
    double total_rd=0.0; double total_es = 0.0; double total_polar=0.0;

    // get all the potentials.
    if (model == "lj" || model == "ljes" || model == "ljespolar") {
        // REPULSION / DISPERSION ENERGY
        total_rd = lj(system);

        if (model == "ljes" || model == "ljespolar") {
            // ELECTROSTATIC ENERGY
            if (system.constants.ewald_es == "on") 
                total_es = coulombic_ewald(system); // using ewald method for es
            else 
                total_es = coulombic(system); // plain old coloumb
            
            if (model == "ljespolar") {
                // POLARIZATION ENERGY
                total_polar = polarization(system); // yikes
            } // end if polar
        } // end if ES
    } // end if LJ
    // POTENTIALS CALCULATED.

 
	static double output[4];
        output[0] = total_rd;
        output[1] = total_es;
        output[2] = total_polar;
	return output;
}
