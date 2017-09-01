#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>


/* the Tang-Toennies potential for the entire system */
double tt(System &system) {
    double potential = 0;
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            for (int k=i+1; k<system.molecules.size(); k++) {
                for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                    
                    double* distances = getDistanceXYZ(system, i,j,k,l);

                }
            }
        }
    }    


    return potential;
}
