#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double getDistance(System &system, int i, int j, int k, int l) {
   // calculate distance between atoms
        double d[3],r,rsq;

        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];

        // apply 1/2 box cutoff:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
        /*
        if (dx > system.constants.x_max) { dx = dx - system.constants.x_length; }
        if (dx < system.constants.x_min) { dx = dx + system.constants.x_length; }
        if (dy > system.constants.y_max) { dy = dy - system.constants.y_length; }
        if (dy < system.constants.y_min) { dy = dy + system.constants.y_length; }
        if (dz > system.constants.z_max) { dz = dz - system.constants.z_length; }
        if (dz < system.constants.z_min) { dz = dz + system.constants.z_length; }
        */    

        rsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        r = sqrt(rsq); 
        return r;
}

double * getDistanceXYZ(System &system, int i, int j, int k, int l) {
    // calculate distance between atoms
        double d[3],r,rsq;

        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];

        // apply 1/2 box cutoff:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
       /*
         if (dx > system.constants.x_max) { dx = dx - system.constants.x_length; }
        if (dx < system.constants.x_min) { dx = dx + system.constants.x_length; }
        if (dy > system.constants.y_max) { dy = dy - system.constants.y_length; }
        if (dy < system.constants.y_min) { dy = dy + system.constants.y_length; }
        if (dz > system.constants.z_max) { dz = dz - system.constants.z_length; }
        if (dz < system.constants.z_min) { dz = dz + system.constants.z_length; }
*/
        rsq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
        r = sqrt(rsq);

        static double output[4];

        for (int n=0; n<3; n++) output[n] = d[n];
        output[3] = r;        

        return output;
}
