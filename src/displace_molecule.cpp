#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>
// rotate function already included in include tree

void displaceMolecule(System &system, string model) {
string movable="notyet";
    int randm = -1;
    while (movable != "M") {
            randm = (rand() % (int)(system.stats.count_movables)) + (int)system.stats.count_frozen_molecules;
            //printf("randm: %i\n",randm);
            movable = system.molecules[randm].MF;
    }           
	system.checkpoint("Got the random molecule.");
	
	double old_V=0.0; double new_V=0.0;

	// first calculate the system's current potential energy
		double* oldpotentials = getTotalPotential(system,model); //new double[2];
		old_V = (oldpotentials[0]+oldpotentials[1]+oldpotentials[2]+oldpotentials[3]);	// keep in K
        //printf("0 1 2 3: %f %f %f %f\n", oldpotentials[0], oldpotentials[1], oldpotentials[2], oldpotentials[3]);
        
    // save a temporary copy of molecule to go back if needed
    Molecule tmp_molecule = system.molecules[randm];

	// pick rotation or displacement
	int rot_disp_flag = -1;  // THIS WILL BE 0 FOR ROTATION, 1 FOR TRANSLATION
	double randRotTrans = (double)rand()/(double)RAND_MAX;
	if (randRotTrans < system.constants.rotate_prob && system.constants.rotate_option == "on") { // do a rotation
		
        // ROTATION
        system.checkpoint("doing a rotation move.");
		rot_disp_flag = 0; system.stats.rotate_attempts++;
		double randangle; string plane; 
		randangle = system.constants.rotate_angle_factor*(double)rand()/(double)RAND_MAX; // angle of rotation from 0 -> rotate_angle_factor
		double randxyz = (double)rand()/(double)RAND_MAX;
		// 1/3 change for a given plane
		if (randxyz < 0.33333) {
			plane="x";
		} else if (randxyz > 0.33333 && randxyz < 0.66667) {
			plane="y";
		} else { plane="z";}

		for (int i=0; i<system.molecules[randm].atoms.size(); i++) {
				double* rotated = rotatePoint(system, system.molecules[randm].atoms[i].pos[0], system.molecules[randm].atoms[i].pos[1], system.molecules[randm].atoms[i].pos[2], plane, randangle);
				system.molecules[randm].atoms[i].pos[0] = rotated[0];
				system.molecules[randm].atoms[i].pos[1] = rotated[1];
				system.molecules[randm].atoms[i].pos[2] = rotated[2];
		} // end for i
        //printf("rotate success.\n");
	} // end rotation option 
	else {
		// DISPLACEMENT
		rot_disp_flag = 1; system.stats.displace_attempts++;
        system.checkpoint("doing displacement move.");
		double randx,randy,randz;
        	randx = system.constants.displace_factor * (((double)rand() / (double)RAND_MAX)*2-1);
        	randy = system.constants.displace_factor * (((double)rand() / (double)RAND_MAX)*2-1);
        	randz = system.constants.displace_factor * (((double)rand() / (double)RAND_MAX)*2-1);

			for (int i=0; i<system.molecules[randm].atoms.size(); i++) {
					system.molecules[randm].atoms[i].pos[0] += randx;
        			system.molecules[randm].atoms[i].pos[1] += randy;
        			system.molecules[randm].atoms[i].pos[2] += randz;
			}
    } // end rotation/displacement if/else

	// check P.B.C. (move the molecule back in the box if needed)
    checkInTheBox(system, randm);

		double* newpotentials = getTotalPotential(system,model); //new double[2];
                new_V = (newpotentials[0]+newpotentials[1]+newpotentials[2]+newpotentials[3]); // keep in K

	double potential_diff = new_V - old_V;

	// now accept or reject the move based on Boltzmann probability
	double boltzmann_factor;
	if (system.constants.ensemble != "nve") {
	        boltzmann_factor = exp(-potential_diff/(system.constants.temp)); // no kb because energy is already in K
	   
	} else if (system.constants.ensemble == "nve") { // for nve only
		boltzmann_factor = pow((system.constants.total_energy - new_V),3.0*system.stats.count_movables/2.0) / pow((system.constants.total_energy - old_V),(3.0*system.stats.count_movables/2.0)); // again, no kb because energy already in K

	}    
            if (boltzmann_factor != INFINITY) {
            if (rot_disp_flag == 0)
                system.stats.rotate_bf_sum += boltzmann_factor;
            else if (rot_disp_flag == 1)
                system.stats.displace_bf_sum += boltzmann_factor;
            } // end BF not inf

	// make ranf for probability pick
	double ranf = (double)rand() / (double)RAND_MAX; // a value between 0 and 1
	//printf("ranf: %10f\n",ranf);

	// apply selection Frenkel Smit p. 30 
	// accept move
	if (ranf < boltzmann_factor) {
		if (rot_disp_flag == 0)
			system.stats.rotate_accepts++;
		else if (rot_disp_flag == 1)
			system.stats.displace_accepts++;
		//printf("accepted\n");
	}
	// reject move (by moving it back)
	else {
		// for whole molecule
		for (int i=0; i<system.molecules[randm].atoms.size(); i++) {
			system.molecules[randm].atoms[i].pos[0] = tmp_molecule.atoms[i].pos[0];
            system.molecules[randm].atoms[i].pos[1] = tmp_molecule.atoms[i].pos[1];
            system.molecules[randm].atoms[i].pos[2] = tmp_molecule.atoms[i].pos[2];
		}

            // check P.B.C. (move the molecule back in the box if needed)
            checkInTheBox(system, randm);
	
	} // end displace/rotate (all ensembles)
    return;
}
