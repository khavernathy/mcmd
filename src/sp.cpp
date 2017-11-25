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
    double q, x, y, z, comx, comy, comz, ax,ay,az; // temporary charge holder, pos. vector
    double multipole0 = 0; // charge sum
    double multipole1[3] = {0,0,0}; // dipole moment
    double dipole_magnitude = 0;
    double multipole2[6] = {0,0,0,0,0,0}; // quadrupole moment
        // 0: xx; 1: yy: 2: zz; 3: xy; 4: xz; 5: yz

    printf("Input atoms: %s\n=============================================\n\n", system.constants.atom_file.c_str());

    system.molecules[0].calc_center_of_mass();

    // multiple 0 (total charge)
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        q = system.molecules[0].atoms[i].C/system.constants.E2REDUCED;
        multipole0 += q;
        x = system.molecules[0].atoms[i].pos[0];
        y = system.molecules[0].atoms[i].pos[1];
        z = system.molecules[0].atoms[i].pos[2];

        comx = system.molecules[0].com[0];
        comy = system.molecules[0].com[1];
        comz = system.molecules[0].com[2];

        ax = x-comx; ay = y-comy; az=z-comz;

        multipole1[0] += q*ax;
        multipole1[1] += q*ay;
        multipole1[2] += q*az;

        multipole2[0] += q*(1.5*(ax*ax) - 0.5*(ax*ax + ay*ay + az*az));
        multipole2[1] += q*(1.5*(ay*ay) - 0.5*(ax*ax + ay*ay + az*az));
        multipole2[2] += q*(1.5*(az*az) - 0.5*(ax*ax + ay*ay + az*az));
        multipole2[3] += q*1.5*ax*ay;
        multipole2[4] += q*1.5*ax*az;
        multipole2[5] += q*1.5*ay*az;
    }   

    for (n=0; n<3; n++)
        multipole1[n] *= system.constants.eA2D;

    for (n=0; n<3; n++)
        dipole_magnitude += multipole1[n]*multipole1[n];
    dipole_magnitude = sqrt(dipole_magnitude);

    for (n=0; n<6; n++)
        multipole2[n] *= system.constants.eA2D;

    printf("Multipole 0 (total charge)      = { %5.5f } e\n", multipole0);
    printf("Multipole 1 (dipole moment)     = { %5.5f %5.5f %5.5f } Debye\n", multipole1[0], multipole1[1], multipole1[2]);
    printf("            (dipole magnitude)  = { %5.5f } Debye\n", dipole_magnitude);
    printf("Multipole 2 (quadrupole moment) = { %5.5f %5.5f %5.5f } Debye A\n", multipole2[0], multipole2[3], multipole2[4]);
    printf("                                  { %5.5f %5.5f %5.5f } Debye A\n", multipole2[3], multipole2[1], multipole2[5]);
    printf("                                  { %5.5f %5.5f %5.5f } Debye A\n", multipole2[4], multipole2[5], multipole2[2]);


    return;
}
