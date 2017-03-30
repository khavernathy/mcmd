#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

void make_pairs(System &system) {
    int i,j,k,l,p,q,r;
    int molsize = (int)system.molecules.size();
    int frozenatoms = system.stats.count_frozens;

    if (system.stats.MCstep == 0 || system.constants.ensemble != ENSEMBLE_NVT) {
        // (re-)set up sizes for 4D array
        system.pairs.resize(molsize);
        for (i=0; i<molsize; ++i) {
            system.pairs[i].resize(frozenatoms);
            for (j=0; j< frozenatoms; ++j) {
                system.pairs[i][j].resize(molsize);
                for (k=0; k<molsize; ++k) {
                    system.pairs[i][j][k].resize(frozenatoms);
                }
            }
        }
    }
    // done resizing distarray.

    // gets all the needed pairwise distances for the entire system in one shot.
    for (i=0; i<molsize; i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (k=i; k<molsize; k++) {
                for (l=0; l<system.molecules[k].atoms.size(); l++) {
                    if (i==k && j==l) continue; // always skip self-atom distance
                    if (system.molecules[i].frozen && system.molecules[k].frozen) continue; // I don't think I ever use frozen-pair distances.    
                    double *distances = getDistanceXYZ(system, i,j,k,l);
                    double r = distances[3];
                    system.pairs[i][j][k][l].r = r;
                    for (p=0; p<3; p++) system.pairs[i][j][k][l].d[p] = distances[p];                    

                    if (r != system.pairs[i][j][k][l].prev_r) {
                        // mark for recalculate
                        system.pairs[i][j][k][l].recalculate = 1;
                    } else {
                        system.pairs[i][j][k][l].recalculate = 0;
                        system.pairs[i][j][k][l].prev_r = r;
                        for (p=0;p<3;p++) system.pairs[i][j][k][l].prev_d[p] = distances[p];
                    }


                    //printf("system.dist.distarray[%i][%i][%i][%i] = %f\n", i,j,k,l, system.dist.distarray[i][j][k][l]);
                }
            }
        }
    }

}

