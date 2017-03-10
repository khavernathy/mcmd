#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>
#include <rotatepoint.cpp>

void addMolecule(System &system, string model) {

    system.checkpoint("starting addMolecule");
	system.stats.insert_attempts++;

	// get current energy.
	double* old_potentials = getTotalPotential(system, model);
	double old_potential = old_potentials[0] + old_potentials[1] + old_potentials[2] + old_potentials[3];
	
	system.molecules.push_back(system.proto);
	system.stats.count_movables += 1;	

	//id in vector INDEX. .ID is PDB ID
	int last_molecule_id = (int)system.molecules.size()-1;
	system.molecules[last_molecule_id].PDBID = system.molecules[last_molecule_id -1].PDBID + 1; // .PDBID is the last one +1
	int last_molecule_PDBID = system.molecules[last_molecule_id].PDBID;
	//printf("The last (added) molecule id is %i\n", last_molecule_id);
	//printf("And its .ID is %i\n", last_molecule_ID);
    
    // set mol_id variable on individual atoms
    for (int i=0; i<system.molecules[last_molecule_id].atoms.size(); i++) {
        system.molecules[last_molecule_id].atoms[i].mol_PDBID = last_molecule_PDBID;
        system.molecules[last_molecule_id].atoms[i].PDBID = system.constants.total_atoms + 1; 
        system.constants.total_atoms += 1;
        // make sure charge is there
        //printf("the added molecule atom %i has charge = %f\n", i, system.molecules[last_molecule_id].atoms[i].C);
    }

	// for random placement
	double randx = 0.5*system.pbc.x_length*(((double)rand()/(double)RAND_MAX)*2-1); // (-1 -> 1) * 1/2length
	double randy = 0.5*system.pbc.y_length*(((double)rand()/(double)RAND_MAX)*2-1);
	double randz = 0.5*system.pbc.z_length*(((double)rand()/(double)RAND_MAX)*2-1);
	double deltax = randx - system.proto.atoms[0].pos[0]; // fixed random distances from the prototype molecule.
	double deltay = randy - system.proto.atoms[0].pos[1];
	double deltaz = randz - system.proto.atoms[0].pos[2];

	// translate the new molecule's atoms to random place.
	for (int i=0; i<system.molecules[last_molecule_id].atoms.size(); i++) {
		system.molecules[last_molecule_id].atoms[i].pos[0] += deltax;
		system.molecules[last_molecule_id].atoms[i].pos[1] += deltay;
		system.molecules[last_molecule_id].atoms[i].pos[2] += deltaz;	
	}
	
	// rotate the molecule here by random amount.
	// get random plane of rotation
		string plane;
		double randp = (double)rand()/(double)RAND_MAX;
		if (randp < 0.33333) plane="x";	
		else if (randp < 0.66667) plane="y";
		else plane="z";
		// get random angle from 0 -> 360
		double randangle = system.constants.rotate_angle_factor * (double)rand()/(double)RAND_MAX;
	// rotate atoms in molecule.
	for (int i=0; i<system.molecules[last_molecule_id].atoms.size(); i++) {
					
		double* rotated_points = rotatePoint(system, system.molecules[last_molecule_id].atoms[i].pos[0], system.molecules[last_molecule_id].atoms[i].pos[1], system.molecules[last_molecule_id].atoms[i].pos[2], plane, randangle);
	
		system.molecules[last_molecule_id].atoms[i].pos[0] = rotated_points[0];
		system.molecules[last_molecule_id].atoms[i].pos[1] = rotated_points[1];
		system.molecules[last_molecule_id].atoms[i].pos[2] = rotated_points[2];	
	}

    // **IMPORTANT: MAKE SURE THE MOLECULE IS IN THE BOX** 
    // maybe not needed actually. displaces move the molecules into box methinks 
	
	// FULLY DONE ADDING MOLECULE TO SYSTEM IN PLACE. NOW GET NEW ENERGY		
	double* new_potentials = getTotalPotential(system, model);
	double new_potential = new_potentials[0] + new_potentials[1] + new_potentials[2] + new_potentials[3];
	double energy_delta = (new_potential - old_potential);
	
	// BOLTZMANN ACCEPT OR REJECT
        double boltz_factor = system.pbc.volume * system.constants.pres * system.constants.ATM2REDUCED/(system.constants.temp * (double)(system.stats.count_movables)) *
                            exp(-energy_delta/system.constants.temp);

	system.stats.insert_bf_sum += boltz_factor;				

	double ranf = (double)rand()/(double)RAND_MAX;
	if (ranf < boltz_factor) {
		//printf("accepted add.\n");
		system.stats.insert_accepts++; //accept (keeps new molecule)
		//system.stats.count_movables++;
	} else {
		// remove the new molecule.
		//printf("rejected add.\n");
		system.molecules.pop_back(); 
		system.constants.total_atoms -= (int)system.proto.atoms.size();
		system.stats.count_movables--;
	//	printf("new size of sys.molecules after erasing: %i\n", (int)system.molecules.size());
	}
    system.checkpoint("done with addMolecule");
return;
}
