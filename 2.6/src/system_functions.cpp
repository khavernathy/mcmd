#include <iostream>
#include <string>
#include <strings.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
using namespace std;


/* MOVE ALL ATOMS SUCH THAT THEY ARE CENTERED ABOUT 0,0,0 */
void centerCoordinates(System &system) {
	int size = system.constants.total_atoms; //system.atoms.size();
		
	double xtemp[size];
	double ytemp[size];
	double ztemp[size];
	int temp_index=0;

	for (int j=0; j<system.molecules.size(); j++) {
	for (int i=0; i<system.molecules[j].atoms.size(); i++) {
		        xtemp[temp_index] = system.molecules[j].atoms[i].pos[0];
                ytemp[temp_index] = system.molecules[j].atoms[i].pos[1];
                ztemp[temp_index] = system.molecules[j].atoms[i].pos[2];
	    temp_index++;
    }
	}

	double xmax= *std::max_element(xtemp, xtemp+size);
	double ymax= *std::max_element(ytemp, ytemp+size);
	double zmax= *std::max_element(ztemp, ztemp+size);

	double xmin= *std::min_element(xtemp, xtemp+size);
	double ymin= *std::min_element(ytemp, ytemp+size);
	double zmin= *std::min_element(ztemp, ztemp+size);

	//cout << xmax << ymax << zmax << xmin << ymin << zmin;
	for (int i=0; i<system.molecules.size(); i++) {
	for (int j=0; j<system.molecules[i].atoms.size(); j++) {
		system.molecules[i].atoms[j].pos[0] = system.molecules[i].atoms[j].pos[0] - (xmin + (xmax - xmin)/2.0);
		system.molecules[i].atoms[j].pos[1] = system.molecules[i].atoms[j].pos[1] - (ymin + (ymax - ymin)/2.0);
		system.molecules[i].atoms[j].pos[2] = system.molecules[i].atoms[j].pos[2] - (zmin + (zmax - zmin)/2.0);
	}
	}
}


/* CALCULATE CENTER OF MASS OF THE SYSTEM */
double * centerOfMass(System &system) {
   //     long unsigned int size = system.atoms.size();

        double x_mass_sum=0.0; double y_mass_sum=0.0; double z_mass_sum=0.0; double mass_sum=0.0;

	for (int j=0; j<system.molecules.size(); j++) {
        for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                double atom_mass = system.molecules[j].atoms[i].m;
                mass_sum += atom_mass;

                x_mass_sum += system.molecules[j].atoms[i].pos[0]*atom_mass;
                y_mass_sum += system.molecules[j].atoms[i].pos[1]*atom_mass;
                z_mass_sum += system.molecules[j].atoms[i].pos[2]*atom_mass;
        }
	}

        double comx = x_mass_sum/mass_sum;
        double comy = y_mass_sum/mass_sum;
        double comz = z_mass_sum/mass_sum;

        double* output = new double[3];
        output[0] = comx; 
        output[1] = comy; 
        output[2] = comz;
        return output;
}


