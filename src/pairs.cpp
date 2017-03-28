#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

void make_pairs(System &system) {
    // we need to call this any time N changes..

    int i,j,k,l,p;
    int count=0;
    int currentsize = (int)system.pairs.size();
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (k=0; k<system.molecules.size(); k++) {
                for (l=0; l<system.molecules[k].atoms.size(); l++) {

                    if (count < currentsize) {
                       // the pair already exists.
                        double* distances = getDistanceXYZ(system, i,j,k,l);     
                        system.pairs[count].r = distances[3];
                        for (p=0; p<3; p++) system.pairs[count].d[p] = distances[p];

                        if (system.pairs[count].r == system.pairs[count].prev_r) 
                            system.pairs[count].recalculate=0;
                        else { 
                            system.pairs[count].recalculate=1;
                            system.pairs[count].prev_r = system.pairs[count].r;
                            for (p=0; p<3; p++) system.pairs[count].prev_d[p] = system.pairs[count].d[p];
                        }

                    } else {
                        // the pair doesn't exist yet.
                        Pair newpair;
                        newpair.id_set = {i,j,k,l};
                    
                        double* distances = getDistanceXYZ(system, i,j,k,l);
                        newpair.r = distances[3];
                        for (p=0; p<3; p++) newpair.d[p] = distances[p]; 

                        // since default for recalculate = 1, no need to set for new pair.

                        system.pairs.push_back(newpair);
                    }
                    
                   // printf("pair recalculate = %i\n",system.pairs[count].recalculate);           
 
                    count++;
                }
            }
        }
    }
    // trim excess pairs if we deleted some atoms in last MC move.
    if (count < currentsize) {
        while (currentsize != (count)) {
            system.pairs.pop_back();
            currentsize--; 
        } 
    }  
    //printf("current size: %i; pairs now: %i\n", (int)system.pairs.size(), count);
}

