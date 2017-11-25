#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

/*
 * Functions for single-point energy of a given molecule
*/

void singlePointEnergy(System &system) {

    int i; // atom id
    int n; // dimensions
    double multipole0 = 0; // charge sum
    double multipole1[3] = {0,0,0}; // dipole moment
    double dipole_magnitude = 0;

    system.molecules[0].calc_center_of_mass();

    // multiple 0 (total charge)
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        multipole0 += system.molecules[0].atoms[i].C/system.constants.E2REDUCED;
        
        for (n=0; n<3; n++)
            multipole1[n] += system.molecules[0].atoms[i].C/system.constants.E2REDUCED * (system.molecules[0].atoms[i].pos[n] - system.molecules[0].com[n]) * system.constants.eA2D;

    }   

    for (n=0; n<3; n++)
        dipole_magnitude += multipole1[n]*multipole1[n];
    dipole_magnitude = sqrt(dipole_magnitude);

    printf("Multipole 0 (total charge)  = { %f } e\n", multipole0);
    printf("Multipole 1 (dipole moment) = { %f %f %f } Debye\n", multipole1[0], multipole1[1], multipole1[2]);
    printf("                  magnitude = %f Debye\n", dipole_magnitude);

    


    return;
}
