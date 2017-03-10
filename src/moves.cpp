#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

/* ADD A MOLECULE */
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
    double boltz_factor = get_boltzmann_factor(system, old_potential, new_potential, "add");     

	double ranf = (double)rand()/(double)RAND_MAX;
	if (ranf < boltz_factor) {
		system.stats.insert_accepts++; //accept (keeps new molecule)
	} else {
		// remove the new molecule.
		system.molecules.pop_back(); 
		system.constants.total_atoms -= (int)system.proto.atoms.size();
		system.stats.count_movables--;
	}
    system.checkpoint("done with addMolecule");
return;
}


/* REMOVE A MOLECULE */
void removeMolecule(System &system, string model) {
    system.checkpoint("starting removeMolecule");
    system.stats.remove_attempts++;
    
    if ((int)system.stats.count_movables == 1) // IMPORTANT: CANCEL THE DELETE IF ONLY 1 MOVABLE MOLECULE LEFT
        return;
   
    // get original energy.
    double* old_potentials = getTotalPotential(system, model);
    double old_potential = old_potentials[0] + old_potentials[1] + old_potentials[2];

    system.checkpoint("getting random movable.");
    // select random movable molecule
    string movable="notyet";
    int randm = -1;
    while (movable != "M") {
            randm = (rand() % (int)(system.stats.count_movables)) + (int)system.stats.count_frozen_molecules;
          //  printf("randm: %i\n",randm);
            movable = system.molecules[randm].MF;
    }
    //printf("The molecule id to be deleted is %i\n",randm);
    system.checkpoint("random movable selected.");

    // save a copy of this moleucule.
    Molecule tmp_molecule = system.molecules[randm];

    // delete the molecule
    system.molecules.erase(system.molecules.begin() + randm);
    system.stats.count_movables--;
    system.constants.total_atoms -= (int)system.proto.atoms.size();

    // get new energy
    double* new_potentials = getTotalPotential(system, model);
    double new_potential = new_potentials[0] + new_potentials[1] + new_potentials[2];
    double energy_delta = (new_potential - old_potential);

    //printf("doing boltzmann -- ");
    // calculate BOLTZMANN FACTOR
    double boltz_factor = get_boltzmann_factor(system, old_potential, new_potential, "remove");

    // accept or reject
    double ranf = (double)rand()/(double)RAND_MAX;
    if (ranf < boltz_factor) {
	    //printf("accepted remove.\n");
	    system.stats.remove_accepts++;
    } else {
	    //printf("rejected remove.\n");
	    // put the molecule back.
	    system.molecules.push_back(tmp_molecule);
        system.constants.total_atoms += (int)system.proto.atoms.size();
	    system.stats.count_movables++;
	}	 // end boltz accept/reject
return;
}


/* DISPLACE (TRANSLATE OR ROTATE) */
void displaceMolecule(System &system, string model) {
    system.stats.displace_attempts++;
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
        
    // save a temporary copy of molecule to go back if needed
    Molecule tmp_molecule = system.molecules[randm];

	// pick rotation or displacement

    if (system.molecules[randm].atoms.size() > 1) { // try rotation
	double randRotTrans = (double)rand()/(double)RAND_MAX;
	if (randRotTrans < system.constants.rotate_prob && system.constants.rotate_option == "on") { // do a rotation
		
        // ROTATION
        system.checkpoint("doing a rotation move.");
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
	} // end rotation option
    } // end if atoms-in-mol > 1 
	else {
		// DISPLACEMENT
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

	// now accept or reject the move based on Boltzmann probability
	double boltzmann_factor = get_boltzmann_factor(system, old_V, new_V, "displace");

	// make ranf for probability pick
	double ranf = (double)rand() / (double)RAND_MAX; // a value between 0 and 1

	// apply selection Frenkel Smit p. 30 
	// accept move
	if (ranf < boltzmann_factor) {
			system.stats.displace_accepts++;
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
