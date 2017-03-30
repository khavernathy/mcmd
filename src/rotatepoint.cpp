#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

double * rotatePoint(System &system, double x, double y, double z, int plane, double angle) {

	double finalx, finaly, finalz;// printf("%f %f %f\n",x,y,z);
    // the function takes an ANGLE 0->360 so need to preconvert from rads if neededo
	angle = angle*M_PI/180.0;
	if (plane == 0) {
		finalx = x;
		finaly = y*cos(angle) -z*sin(angle);
		finalz = y*sin(angle) + z*cos(angle);
	} else if (plane == 1) {
		finalx = x*cos(angle) + z*sin(angle);
		finaly = y;
		finalz = -x*sin(angle) + z*cos(angle); 
	} else if (plane == 2) {
		finalx = x*cos(angle) - y*sin(angle);
		finaly = x*sin(angle) + y*cos(angle);
		finalz = z;
	}

	static double output[3]; //double* output[3]; //vector<double> output; //printf("%f %f %f\n",finalx,finaly,finalz);
        output[0] = finalx;
        output[1] = finaly;
	    output[2] = finalz;
        return output;
}

double * rotatePointRadians(System &system, double x, double y, double z, int plane, double angle) {

	double finalx, finaly, finalz;// printf("%f %f %f\n",x,y,z);
    // the function takes an ANGLE 0->360 so need to preconvert from rads if neededo
	
	if (plane == 0) {
		finalx = x;
		finaly = y*cos(angle) -z*sin(angle);
		finalz = y*sin(angle) + z*cos(angle);
	} else if (plane == 1) {
		finalx = x*cos(angle) + z*sin(angle);
		finaly = y;
		finalz = -x*sin(angle) + z*cos(angle); 
	} else if (plane == 2) {
		finalx = x*cos(angle) - y*sin(angle);
		finaly = x*sin(angle) + y*cos(angle);
		finalz = z;
	}

	static double output[3]; //double* output[3]; //vector<double> output; //printf("%f %f %f\n",finalx,finaly,finalz);
        output[0] = finalx;
        output[1] = finaly;
	    output[2] = finalz;
        return output;
}
// -------------- end rotatePoint func ---------------------

/*
void rotateMolecule(System &system, Molecule mol, string plane, double angle) {
    angle = angle*M_PI/180.0;

    // for each atom in molecule, rotate it.
    for (int i=0; i<mol.atoms.size(); i++) {
        
        double* rotated = rotatePoint(system, mol.atoms[i].pos[0], mol.atoms[i].pos[1], mol.atoms[i].pos[2], plane, angle);

        for (int n=0; n<3; n++) mol.atoms[i].pos[n] = rotated[n];
         
    }
}*/
