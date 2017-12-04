/* Douglas M Franz
 * Space group, USF, 2017
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

// function to determine UFF atom-type based on
// name of element and number of bonds
string getUFFlabel(System &system, string name, int num_bonds) {
    return "hey";
}

// function to find all bonds for all atoms.
void findBonds(System &system) {
    int i,j,l;
    double r;
    int local_bonds=0;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            local_bonds = 0;
            // for each atom, we find its bonded neighbors by a distance search
            // (by only searching atoms on this same molecule)
            for (l=0; l<system.molecules[i].atoms.size(); l++) {
               if (j==l) continue; // don't do self-atom
               double* distances = getDistanceXYZ(system, i,j,i,l);
               r = distances[3];
               if (r < system.constants.bondlength) {
                    local_bonds++;
                    system.molecules[i].atoms[j].bonds.insert(std::pair<int,double>(l,r));
               } // end if r < bond-length      
            } // end pair (i,j) -- (i,l)  
         
            // based on the total number of bonds to this atom, 
            // determine the atom-type from UFF.
            system.molecules[i].atoms[j].UFFlabel = getUFFlabel(system, system.molecules[i].atoms[j].name, system.molecules[i].atoms[j].bonds.size()); 

        } // end j
    } // end i
}

/*
 * Essentially everything below comes from the parameters and
 * functional forms prescribed in UFF, via
 * J. Am. Chem. Soc. Vol. 114, No. 25, 1992
 * */

// get the total potential from bond stretches
// via the Morse potential
double stretch_energy(System &system) {
    return 0;
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_energy(System &system) {
    return 0;
}

// get the total potential from torsions
// via simple Fourier small cosine expansion
double torsions_energy(System &system) {
    return 0;
}




