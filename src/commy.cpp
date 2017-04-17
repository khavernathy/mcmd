#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>


double commy(System &system) {
    // the communist potential from 1961 van der Waals paper
    double total_pot=0; 
    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    int i,j,k,l;
    double r,r7;
    double numerator= system.constants.HBARC * -23.0;
    double FourPi = 12.566370614359172;
    double polar1, polar2;    

    for (i = 0; i < system.molecules.size(); i++) {
    for (j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (k = i+1; k < system.molecules.size(); k++) { // so if one frozen molecule, frozen-frozen is ignored.
    for (l =0; l < system.molecules[k].atoms.size(); l++) {

        polar1 = system.molecules[i].atoms[j].polar;
        polar2 = system.molecules[k].atoms[l].polar;

        if (polar1 == 0. || polar2 == 0.) continue; // skip 0-energy

        // calculate distance between atoms
        double* distances = getDistanceXYZ(system, i, j, k, l);
        r = distances[3]; //printf("r %f\n", r);
        r7= r*r; // r2
        r7 *= r7 * r7; // r6
        r7 *= r; // r7

    

        if ((r <= cutoff)) {
    
            //printf("num %f denom %f polar1 %f polar2 %f\n", numerator, FourPi*r7, system.molecules[i].atoms[j].polar, system.molecules[k].atoms[l].polar);
            total_pot += (numerator/(FourPi*r7)) * polar1 * polar2;
        }
        
    }  // loop l
    } // loop k 
    } //loop j
    } // loop i


//    printf("total commy: %e \n", total_pot);
    return total_pot;

}



