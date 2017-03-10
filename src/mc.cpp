#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>
#include <potential.cpp>
#include <add_molecule.cpp>
#include <remove_molecule.cpp>
#include <displace_molecule.cpp>
#include <boltzmann.cpp>

// PHAST2 NOT INCLUDED YET

/* (RE)DEFINE THE BOX LENGTHS */
void defineBox(System &system) { // takes input in A
	// easy 90 90 90 systems
	if (system.pbc.alpha == 90 && system.pbc.beta == 90 && system.pbc.gamma == 90) {
		// assumes x_length, y_length, z_length are defined already in system.
		// i.e. the volume-change function does that before calling this function

		system.pbc.basis[0][0] = system.pbc.x_length;
		system.pbc.basis[1][1] = system.pbc.y_length;   
		system.pbc.basis[2][2] = system.pbc.z_length;

		system.pbc.x_max = system.pbc.x_length/2.0;
		system.pbc.x_min = -system.pbc.x_max;
		system.pbc.y_max = system.pbc.y_length/2.0;
		system.pbc.y_min = -system.pbc.y_max;
		system.pbc.z_max = system.pbc.z_length/2.0;
		system.pbc.z_min = -system.pbc.z_max;

		system.pbc.cutoff = system.pbc.x_max; // crude method but good for cubes (and only cubes!!)
		// need to make basis vector system later.
		system.constants.ewald_alpha = 3.5/system.constants.cutoff; // update ewald_alpha if we have a vol change

		system.pbc.volume = system.pbc.x_length * system.pbc.y_length * system.pbc.z_length; // in A^3
	}
	// universal definitions
	else {
		// this could forseeably be a problem if, for some weird reason, someone wants to do NPT with a weird box.
		system.pbc.calcVolume();
		system.pbc.calcRecip();
		system.pbc.calcCutoff();
		system.pbc.calcBoxVertices();
		system.pbc.calcPlanes();
	}
}

/* CHANGE VOLUME OF SYSTEM BY BOLTZMANN PROB -- FOR NPT */
void changeVolumeMove(System &system) {
	// increasing volume will increase energy supplied by external pressure but also allow more room for lower potential among atoms
	system.stats.volume_attempts++;
	// generate small randam distance change for volume adjustment
	double ranf = (double)rand() / (double)RAND_MAX; // for boltz check	
    double ranv = (double)rand() / (double)RAND_MAX; // for volume change
    double old_energy, new_energy;

        double* energies = getTotalPotential(system,system.constants.potential_form);
        old_energy=energies[0]+energies[1]+energies[2]+energies[3];

    double old_side = system.pbc.x_length; // in A; save just in case, to reset if move rejected.
    system.pbc.old_volume = system.pbc.volume;

    // change the volume to test new energy.
    double new_volume = exp(log(system.pbc.volume) + (ranv-0.5)*system.constants.volume_change);// mpmc default = 2.5

    //double new_volume = 1e-30*new_volume_A3;
    double new_side = pow(new_volume, 1.0/3.0);
    system.pbc.x_length = new_side;
    system.pbc.y_length = new_side;
    system.pbc.z_length = new_side;
    defineBox(system);
    //printf("OLDV before energy call: %f; NEWV before energy call: %f\n", system.pbc.old_volume, system.pbc.volume);        
        double* potentials = getTotalPotential(system,system.constants.potential_form);
        new_energy=potentials[0]+potentials[1]+potentials[2]+potentials[3];

    //printf("OLDV before bf call: %f; NEWV before bf call: %f\n", system.pbc.old_volume, system.pbc.volume); 
    double boltzmann_factor = get_boltzmann_factor(system, old_energy, new_energy, "volume");

    //printf("ranf < bf? %f < %f ?\n", ranf, boltzmann_factor);
	if (ranf < boltzmann_factor) {
		// accept move
        //printf("ACCEPTED\n");
		system.stats.volume_change_accepts++;
	} else {
        //printf("REJECTED\n");
		// reject move (move volume back)
        system.pbc.x_length = old_side;
        system.pbc.y_length = old_side;
        system.pbc.z_length = old_side;
        defineBox(system);
	}
}


// ================== MAIN MC FUNCTION. HANDLES MOVE TYPES AND BOLTZMANN ACCEPTANCE =================
// ACCORDING TO DIFFERENT ENSEMBLES
void runMonteCarloStep(System &system, string model) {
    system.checkpoint("Entered runMonteCarloStep().");

	// VOLUME MOVE (only NPT)
	if (system.constants.ensemble == "npt") {
		double VCP = system.constants.vcp_factor/(double)system.stats.count_movables; // Volume Change Probability
		double ranf = (double)rand() / (double)RAND_MAX; // between 0 and 1
        	if (ranf < VCP) {
                    system.checkpoint("doing a volume change move.");
                	changeVolumeMove(system);
                    system.checkpoint("done with volume change move.");
			return; // we tried a volume change, so exit MC step.
		}
	}

	// ADD / REMOVE (only uVT)
	// We'll choose 0-0.5 for add; 0.5-1 for remove (equal prob.)
	if (system.constants.ensemble == "uvt") {
		double IRP = system.constants.insert_factor;  
		double ranf = (double)rand() / (double)RAND_MAX; // between 0 and 1
		if (ranf < IRP) {
            //system.checkpoint("doing an add or remove.");
			// we're going to do an add/remove now.
			double ranf2 = (double)rand() / (double)RAND_MAX; // 0->1
			// ADD A MOLECULE
			if (ranf2 < 0.5) {
                system.checkpoint("doing molecule add move.");
				addMolecule(system, model);			
                system.checkpoint("done with molecule add move.");
			} // end add
 
			else { // REMOVE MOLECULE
                system.checkpoint("doing molecule delete move.");
				removeMolecule(system, model);
                system.checkpoint("done with molecule delete move.");
			} // end add vs. remove		
		return; // we did the add or remove so exit MC step.
		} // end doing an add/remove. 
	} // end if uVT		

	
	// DISPLACE / ROTATE (for all: NPT, uVT, NVT, NVE); NVE has special BoltzFact tho.
	// make sure it's a movable molecule
	system.checkpoint("NOT volume/add/remove :: Starting displace or rotate..");
    displaceMolecule(system, model);


    	system.checkpoint("done with displace/rotate");
    return; // done with move, so exit MC step	
}
