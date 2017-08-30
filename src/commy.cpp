#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

/*
    note: this should not be used for publishable research.
    It's just a toy-potential for testing different 
    mathematical forms 
*/

double commy(System &system) {
    // the communist potential from 1961 van der Waals paper
    // http://iopscience.iop.org/article/10.1070/PU1961v004n02ABEH003330/meta;jsessionid=F9941A012802331F1E4FF488A0F004A9.c4.iopscience.cld.iop.org  
    double total_pot=0; 
    const double cutoff = system.pbc.cutoff;
    int i,j,k,l;
    double r,sr6,r7;
    const double numerator= system.constants.HBARC * -23.0;
    const double FourPi = 12.566370614359172;
    double polar1, polar2;    
    double attractive, repulsive; // energies
    double eps,sig;

    for (i = 0; i < system.molecules.size(); i++) {
    for (j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (k = i+1; k < system.molecules.size(); k++) { // so if one frozen molecule, frozen-frozen is ignored.
    for (l =0; l < system.molecules[k].atoms.size(); l++) {

        attractive=0; repulsive=0;

        // do mixing rules
        eps = system.molecules[i].atoms[j].eps,sig=system.molecules[i].atoms[j].sig;
        if (eps != system.molecules[k].atoms[l].eps)
            eps = sqrt(eps * system.molecules[k].atoms[l].eps);
        if (sig != system.molecules[k].atoms[l].sig)
            sig = 0.5 * (sig + system.molecules[k].atoms[l].sig);

        polar1 = system.molecules[i].atoms[j].polar;
        polar2 = system.molecules[k].atoms[l].polar;

        // calculate distance between atoms
        double* distances = getDistanceXYZ(system, i, j, k, l);
        r = distances[3]; //printf("r %f\n", r);
        if (sig != 0 && eps != 0) {
            sr6 = sig/r;
            sr6 *= sr6;
            sr6 *= sr6*sr6;
            repulsive = 4.0*eps*(sr6*sr6);
        }

        if (polar1 != 0 && polar2 != 0) {
            r7= r*r; // r2
            r7 *= r7 * r7; // r6
            r7 *= r; // r7
            attractive = ((numerator/(FourPi*r7)) * polar1 * polar2) / 50; // I'm doing the division here to scale down
        }

        if ((r <= cutoff)) {
            // communist potential is attractive contribution; using normal LJ repulsion contribution
            total_pot += attractive + repulsive;
        }
        
    }  // loop l
    } // loop k 
    } //loop j
    } // loop i

    return total_pot;

}



