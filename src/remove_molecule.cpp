#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>


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
    double boltz_factor = system.constants.temp * ((double)(system.stats.count_movables) + 1.0)/(system.constants.volume*system.constants.pres*system.constants.ATM2REDUCED) *
                            exp(-energy_delta/system.constants.temp); 


    //printf("boltzmann success -- ");
    system.stats.remove_bf_sum += boltz_factor;
    //printf("removing bf: %f\n",boltz_factor);

    //printf("Make sure molecules[0] is a thing: %s\n",system.molecules[0].name.c_str());
    //printf("Make sure tmp_atoms[atomids[0]] is a thing: %s \n",tmp_atoms[(int)atomids[0]].name.c_str());
    //printf("Make sure tmp_atoms[0] is a thing: %s \n", tmp_atoms[0].name.c_str());
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
        //printf("put molecule id=%i back.\n",randm);
	    //printf("The re-added molecule has .ID = %i\n",system.molecules[(int)(system.molecules.size()-1)].ID);
	}	 // end boltz accept/reject
    //system.checkpoint("done with removeMolecule");
return;
}
