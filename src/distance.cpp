#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double getDistance(System &system, int i, int j, int k, int l) {
   // calculate distance between atoms
        double d[3],r,rsq;

        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];

        rsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        r = sqrt(rsq); 
        return r;
}

double * getDistanceXYZ(System &system, int i, int j, int k, int l) {
    // calculate distance between atoms
        double d[3],r,rsq;

        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];

        rsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        r = sqrt(rsq);

        static double output[4];

        for (int n=0; n<3; n++) output[n] = d[n];
        output[3] = r;        

        return output;
}
