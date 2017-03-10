#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double getDistance(System &system, int i, int j, int k, int l) {
   // calculate distance between atoms
        double dx,dy,dz,r,rsq;

        dx = (system.molecules[i].atoms[j].pos[0] - system.molecules[k].atoms[l].pos[0]);
        dy = (system.molecules[i].atoms[j].pos[1] - system.molecules[k].atoms[l].pos[1]);
        dz = (system.molecules[i].atoms[j].pos[2] - system.molecules[k].atoms[l].pos[2]);

        // apply 1/2 box cutoff:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
        /*
        if (dx > system.constants.x_max) { dx = dx - system.constants.x_length; }
        if (dx < system.constants.x_min) { dx = dx + system.constants.x_length; }
        if (dy > system.constants.y_max) { dy = dy - system.constants.y_length; }
        if (dy < system.constants.y_min) { dy = dy + system.constants.y_length; }
        if (dz > system.constants.z_max) { dz = dz - system.constants.z_length; }
        if (dz < system.constants.z_min) { dz = dz + system.constants.z_length; }
        */    

        rsq = dx*dx + dy*dy + dz*dz;
        r = sqrt(rsq); 
        return r;
}

double * getDistanceXYZ(System &system, int i, int j, int k, int l) {
    // calculate distance between atoms
        double dx,dy,dz,r,rsq;

        dx = (system.molecules[i].atoms[j].pos[0] - system.molecules[k].atoms[l].pos[0]);
        dy = (system.molecules[i].atoms[j].pos[1] - system.molecules[k].atoms[l].pos[1]);
        dz = (system.molecules[i].atoms[j].pos[2] - system.molecules[k].atoms[l].pos[2]);

        // apply 1/2 box cutoff:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
       /*
         if (dx > system.constants.x_max) { dx = dx - system.constants.x_length; }
        if (dx < system.constants.x_min) { dx = dx + system.constants.x_length; }
        if (dy > system.constants.y_max) { dy = dy - system.constants.y_length; }
        if (dy < system.constants.y_min) { dy = dy + system.constants.y_length; }
        if (dz > system.constants.z_max) { dz = dz - system.constants.z_length; }
        if (dz < system.constants.z_min) { dz = dz + system.constants.z_length; }
*/
        rsq = dx*dx + dy*dy + dz*dz;
        r = sqrt(rsq);

        static double output[4];

        output[0] = dx;
        output[1] = dy;
        output[2] = dz;
        output[3] = r;        

        return output;
}
