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

double getrand() {
    return (double)rand()/(double)RAND_MAX; // a value between 0 and 1.
}

/* GET BOX LIMIT COORDINATE FOR PBC */
double getBoxLimit(System &system, int planeid, double x, double y, double z) {

    // 0: -x; 1: +x; 2: -y; 3: +y; 4: -z; 5: +z;
    if (planeid == 0) return (-system.pbc.D[0] - system.pbc.B[0]*y - system.pbc.C[0]*z)/system.pbc.A[0]; // -x
    if (planeid == 1) return (-system.pbc.D[1] - system.pbc.B[1]*y - system.pbc.C[1]*z)/system.pbc.A[1]; // +x
    if (planeid == 2) return (-system.pbc.D[2] - system.pbc.A[2]*x - system.pbc.C[2]*z)/system.pbc.B[2]; // -y
    if (planeid == 3) return (-system.pbc.D[3] - system.pbc.A[3]*x - system.pbc.C[3]*z)/system.pbc.B[3]; // +y
    if (planeid == 4) return (-system.pbc.D[4] - system.pbc.A[4]*x - system.pbc.B[4]*y)/system.pbc.C[4]; // -z
    if (planeid == 5) return (-system.pbc.D[5] - system.pbc.A[5]*x - system.pbc.B[5]*y)/system.pbc.C[5]; // +z
}


/* MOVE ALL ATOMS SUCH THAT THEY ARE CENTERED ABOUT 0,0,0 */
void centerCoordinates(System &system) {
	int size = system.constants.total_atoms;	
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


/* CHECK IF MOLECULE IS IN BOX AND MOVE BACK IN IF NOT */
void checkInTheBox(System &system, int i) { // i is molecule id passed in function call.

// first the easy systems, alpha=beta=gamma=90
// saves some time because box limits need not be recomputed
if (system.pbc.alpha == 90 && system.pbc.beta == 90 && system.pbc.gamma == 90) {
    for (int j=0; j<system.molecules[i].atoms.size(); j++) {
        if (system.molecules[i].atoms[j].pos[0]  > system.pbc.x_max) {
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[0] -= system.pbc.x_length;
                }
        }
        else if (system.molecules[i].atoms[j].pos[0] < system.pbc.x_min) {
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
	                system.molecules[i].atoms[k].pos[0] += system.pbc.x_length;
                }
        }

        if (system.molecules[i].atoms[j].pos[1] > system.pbc.y_max) {
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] -= system.pbc.y_length;
                }
        }
        else if (system.molecules[i].atoms[j].pos[1] < system.pbc.y_min) {
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] += system.pbc.y_length;
                }
        }

        if (system.molecules[i].atoms[j].pos[2]  > system.pbc.z_max) {
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] -= system.pbc.z_length;
                }
        }
        else if (system.molecules[i].atoms[j].pos[2] < system.pbc.z_min) {
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] += system.pbc.z_length;
                }
        }	
    } // end atom loop
} // end if alpha=beta=gamma=90

// the universal treatment, alpha != beta ?= gamma
else {
    double comv[3]; // center of mass of molecules
    double box_limit[6]; //-x, +x, -y, +y, -z, +z limits, which are functions of atom position in other dims.
    double newlimit[6]; // an array to hold temporary new planes for movements.

    for (int n=0; n<3; n++) comv[n] = system.molecules[i].com[n];

  
    // find appropriate values for box limits based on atom coordinates.
    // check in this order: (-x, +x, -y, +y, -z, +z)
    for (int n=0; n<6; n++)
        box_limit[n] = getBoxLimit(system, n, comv[0], comv[1], comv[2]);

    //printf("box_limit values ::\n");
    //for (int n=0; n<6; n++) printf("%i : %.5f\n",n,box_limit[n]);

    // if outside box checks are correct but moves are not..
    // NOTE: the r_c cutoff is not actually used in MD

        if (system.molecules[i].com[0] < box_limit[0]) { // left of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);

                double ydiff = system.molecules[i].atoms[k].pos[1] - newlimit[2];
                double zdiff = system.molecules[i].atoms[k].pos[2] - newlimit[4];

                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0] + system.pbc.x_length, system.molecules[i].atoms[k].pos[1] , system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0] + system.pbc.x_length, system.molecules[i].atoms[k].pos[1] , system.molecules[i].atoms[k].pos[2]);

                system.molecules[i].atoms[k].pos[0] += system.pbc.x_length;
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydiff;
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdiff;
             }
        }


        if (system.molecules[i].com[0] > box_limit[1]) { // right of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);

                double ydiff = system.molecules[i].atoms[k].pos[1] - newlimit[2];
                double zdiff = system.molecules[i].atoms[k].pos[2] - newlimit[4];

                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0] - system.pbc.x_length, system.molecules[i].atoms[k].pos[1] , system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0] - system.pbc.x_length, system.molecules[i].atoms[k].pos[1] , system.molecules[i].atoms[k].pos[2]);

                system.molecules[i].atoms[k].pos[0] -= system.pbc.x_length;
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydiff;
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdiff;
            }

        }

        if (system.molecules[i].com[1] < box_limit[2]) { // below box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);

                double xdiff = system.molecules[i].atoms[k].pos[0] - newlimit[0];
                double zdiff = system.molecules[i].atoms[k].pos[2] - newlimit[4];

                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1] + system.pbc.y_length, system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1] + system.pbc.y_length, system.molecules[i].atoms[k].pos[2]);
    
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdiff;
                system.molecules[i].atoms[k].pos[1] += system.pbc.y_length;
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdiff;

            }

        }

        if (system.molecules[i].com[1] > box_limit[3]) { // above box

            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);

                double xdiff = system.molecules[i].atoms[k].pos[0] - newlimit[0];
                double zdiff = system.molecules[i].atoms[k].pos[2] - newlimit[4];

                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1] - system.pbc.y_length, system.molecules[i].atoms[k].pos[2]);
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1] - system.pbc.y_length, system.molecules[i].atoms[k].pos[2]);

    
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdiff;
                system.molecules[i].atoms[k].pos[1] -= system.pbc.y_length;
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdiff;
        
            }
        }

        if (system.molecules[i].com[2] < box_limit[4]) { // behind box
        
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
        
                double xdiff = system.molecules[i].atoms[k].pos[0] - newlimit[0];
                double ydiff = system.molecules[i].atoms[k].pos[1] - newlimit[2];

                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2] + system.pbc.z_length);
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2] + system.pbc.z_length);

                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdiff;
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydiff; 
                system.molecules[i].atoms[k].pos[2] += system.pbc.z_length;

            }
        }

        if (system.molecules[i].com[2] > box_limit[5]) { // in front of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);

                double xdiff = system.molecules[i].atoms[k].pos[0] - newlimit[0];
                double ydiff = system.molecules[i].atoms[k].pos[1] - newlimit[2];

                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2] - system.pbc.z_length);
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2] - system.pbc.z_length);

                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdiff;
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydiff; 
                system.molecules[i].atoms[k].pos[2] -= system.pbc.z_length;

            }
        }


} // end if non-90/90/90 
} // end pbc function


void addAtomToProto(System &system, string name, string molname, string MF, double x, double y, double z, double mass, double charge, double polarizability, double epsilon, double sigma) {
    // initialize
    Atom atom;
    // PDBIDS
    int last_mol_index = system.molecules.size() - 1;
    int last_mol_pdbid = system.molecules[last_mol_index].PDBID;
    int last_atom_pdbid = system.molecules[last_mol_index].atoms[system.molecules[last_mol_index].atoms.size() - 1].PDBID;
    atom.PDBID = (int)(last_atom_pdbid + 1);
    atom.mol_PDBID = (int)(last_mol_pdbid + 1);

    atom.name = name;
    atom.mol_name = molname;
    atom.MF = MF;
    atom.pos[0] = x;
    atom.pos[1] = y;
    atom.pos[2] = z;
    atom.m = mass * system.constants.cM;
    atom.C = charge;
    atom.polar = polarizability;
    atom.eps = epsilon;
    atom.sig = sigma;

    system.proto.mass += atom.m;
    // send it over
    system.proto.atoms.push_back(atom);
    
}

void moleculePrintout(System &system) {
    // CONFIRM ATOMS AND MOLECULES PRINTOUT

    // an excessive full list of atoms/molecules.
/*    
    for (int b=0; b<system.molecules.size(); b++) {
        printf("Molecule PDBID = %i: %s has %i atoms and is %s; The first atom has PDBID = %i\n", system.molecules[b].PDBID, system.molecules[b].name.c_str(), (int)system.molecules[b].atoms.size(), system.molecules[b].MF.c_str(), system.molecules[b].atoms[0].PDBID);
        for (int c=0; c<system.molecules[b].atoms.size(); c++) {
            system.molecules[b].atoms[c].printAll();
            printf("\n");
        }
        
    } 
*/

    if (system.constants.mode == "mc") { // prototype is only used for MC.
    // CHANGE THE PROTOTYPE IF USER SUPPLIED A KEYWORD IN INPUT
    // THIS WILL OVERWRITE ANY PROTOTYPE IN THE INPUT ATOMS FILE if user put it there, e.g. whatever.pdb
        if (system.constants.sorbate_name != "") {
            string sorbmodel = system.constants.sorbate_name;
            
            // clear the prototype (even if there isn't one).
            system.proto.reInitialize(); 
            int last_mol_index = system.molecules.size() -1;
            int last_mol_pdbid = system.molecules[last_mol_index].PDBID;
            system.proto.PDBID = last_mol_pdbid+1;


            //std::cout << "THE SORB MODEL WAS SUPPLIED: " << sorbmodel.c_str(); printf("\n");
            // each call takes 12 arguments
            // HYDROGEN H2
            if (sorbmodel == "h2_buch") {
                addAtomToProto(system, "H2G", "H2", "M", 0.0, 0.0, 0.0, 2.016, 0.0, 0.0, 34.2, 2.96);
                system.proto.name = "H2";
            }
            else if (sorbmodel == "h2_bss") {
                addAtomToProto(system, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.74640, 0.0, 8.85160, 3.2293);
                addAtomToProto(system, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, "H2N", "H2", "M", 0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                addAtomToProto(system, "H2N", "H2", "M", -0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                system.proto.name = "H2";
            }
            else if (sorbmodel == "h2_dl") {
                addAtomToProto(system, "H2G", "H2", "M", 0.00, 0.0, 0.0, -0.93600, 0.0, 36.7, 2.958, 0.0);
                addAtomToProto(system, "H2E", "H2", "M", -0.370, 0.0, 0.0, 1.008, 0.468, 0.0, 0.0, 0.0);
                addAtomToProto(system, "H2E", "H2", "M", 0.370, 0.0, 0.0, 1.008, 0.468, 0.0, 0.0, 0.0);
                system.proto.name = "H2";
            }
            else if (sorbmodel == "h2_bssp") {
                addAtomToProto(system, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.7464, 0.6938, 12.76532, 3.15528);
                addAtomToProto(system, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.3732, 0.00044, 0.0, 0.0);
                addAtomToProto(system, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.00044, 0.0, 0.0);
                addAtomToProto(system, "H2N", "H2", "M", 0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                addAtomToProto(system, "H2N", "H2", "M", -0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                system.proto.name = "H2";
            }
            // CARBON DIOXIDE CO2
            else if (sorbmodel == "co2_phast") {
                addAtomToProto(system, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77106, 0.0, 8.52238, 3.05549);
                addAtomToProto(system, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, "CON", "CO2", "M", 1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                addAtomToProto(system, "CON", "CO2", "M", -1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                system.proto.name = "CO2";
            }
            else if (sorbmodel == "co2_phast*") {
                addAtomToProto(system, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77134, 1.22810, 19.61757, 3.03366);
                addAtomToProto(system, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, "CON", "CO2", "M", 1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                addAtomToProto(system, "CON", "CO2", "M", -1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                system.proto.name = "CO2";
            }
            else if (sorbmodel == "co2_trappe") {
                addAtomToProto(system, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.7, 0.0, 27.0, 2.80);
                addAtomToProto(system, "COE", "CO2", "M", 1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05);
                addAtomToProto(system, "COE", "CO2", "M", -1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05); 
                system.proto.name = "CO2";
            }
            // NITROGEN N2
            else if (sorbmodel == "n2_mcquarrie") {
                addAtomToProto(system, "N2G", "N2", "M", 0.0, 0.0, 0.0, 28.01344, 0.0, 0.0, 95.1, 3.7);
                system.proto.name = "N2";
            }
            // METHANE CH4
            else if (sorbmodel == "ch4_trappe") {
                addAtomToProto(system, "CHG", "CH4", "M", 0.0, 0.0, 0.0, 16.0426, 0.0, 0.0, 148.0, 3.73);
                system.proto.name = "CH4";
            }
            else if (sorbmodel == "ch4_9site") {
                addAtomToProto(system, "CHG", "CH4", "M", 0.0, 0.0, 0.0, 12.011, -0.5868, 0.0, 58.53869, 2.22416);
                addAtomToProto(system, "CHE", "CH4", "M", 0.0, 0.0, 1.099, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system, "CHE", "CH4", "M", 1.036, 0.0, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system, "CHE", "CH4", "M", -0.518, -0.897, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system, "CHE", "CH4", "M", -0.518, 0.897, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system, "MOV", "CH4", "M", 0.0, 0.0, 0.816, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system, "MOV", "CH4", "M", 0.769, 0.0, -0.271, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system, "MOV", "CH4", "M", -0.385, -0.668, -0.271, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system, "MOV", "CH4", "M", -0.385,   0.668,  -0.271,  0.00000,  0.00000,  0.00000, 16.85422,  2.96286);
                system.proto.name = "CH4";
            }
            else if (sorbmodel == "ch4_9site*") {  
                addAtomToProto(system, "CHG",  "CH4", "M", 0.000,  -0.000,  -0.000, 12.01100, -0.58680,  1.09870, 45.09730,  2.16247); //  0.00000  0.00000
                addAtomToProto(system, "CHE", "CH4", "M", 0.000,  0.000,1.099,  1.00790,  0.14670,  0.42460,  0.00000,  0.0000);
                addAtomToProto(system, "CHE",  "CH4", "M",  1.036,  -0.000,  -0.366,  1.00790,  0.14670,  0.42460, 0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system, "CHE",  "CH4", "M", -0.518, -0.897, -0.366,  1.00790,  0.14670,  0.42460,  0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system, "CHE",  "CH4", "M", -0.518,   0.897,  -0.366,  1.00790,  0.14670,  0.42460,  0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system, "MOV", "CH4", "M", 0.000,  -0.000,   0.814,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                addAtomToProto(system, "MOV", "CH4", "M", 0.768,  -0.000,  -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787);//  0.00000  0.00000
                addAtomToProto(system, "MOV",  "CH4", "M", -0.383,  -0.666, -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                addAtomToProto(system, "MOV",  "CH4", "M", -0.383,   0.666,  -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                system.proto.name = "CH4";
            }
            // ACETYLENE C2H2
            else if (sorbmodel == "c2h2") {
                addAtomToProto(system, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                addAtomToProto(system, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                system.proto.name = "ACE";
            }
            else if (sorbmodel == "c2h2*") {
                addAtomToProto(system, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                addAtomToProto(system, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                system.proto.name = "ACE";
            }
            // ETHLYENE C2H4
            else if (sorbmodel == "c2h4") {
                addAtomToProto(system, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                system.proto.name = "ETH";
            }
            else if (sorbmodel == "c2h4*") {
                addAtomToProto(system, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                system.proto.name = "ETH";
            }
            // ETHANE C2H6
            else if (sorbmodel == "c2h6") {
                addAtomToProto(system, "C2G", "ETH", "M", -0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897); //  0.00000  0.00000
                addAtomToProto(system, "C2G", "ETH", "M", 0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                system.proto.name = "ETH";
            }
            else if (sorbmodel == "c2h6*") {
                addAtomToProto(system, "C2G", "ETH", "M", -0.762,0.000, 0.000,12.01100, -0.04722,0.6967,98.63326,3.37151);
                addAtomToProto(system, "C2G", "ETH", "M", 0.762,0.000,0.000,12.01100,-0.04722,0.6967,98.63326,  3.37151);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                system.proto.name = "ETH";
            }
            // WATER H20
            else if (sorbmodel == "h2o") {
                addAtomToProto(system, "OXY", "H2O", "M", 0.0, 0.0, 0.0, 16.0, -0.834, 0.0, 76.42, 3.151);
                addAtomToProto(system, "HYD", "H2O", "M", -0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                addAtomToProto(system, "HYD", "H2O", "M", 0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                system.proto.name = "H2O";
            }

            // USER SORBATE MODEL NOT FOUND; ERROR OUT
            else {
                std::cout << "ERROR: The sorbate model name you supplied, " << sorbmodel.c_str() << ", was not found in the database. Check your spelling or use a manual model in your input atoms file."; printf("\n");
                std::exit(0);
            }
        } // end if sorbate name != ""    
       
        // finally, zero the prototype's coordinates
        system.proto.calc_center_of_mass();
        for (int i=0; i<system.proto.atoms.size(); i++) {
            for (int k=0; k<3; k++) system.proto.atoms[i].pos[k] -= system.proto.com[k];
        }
        system.proto.calc_center_of_mass(); 

        printf("Prototype molecule has PDBID %i ( name %s ) and has %i atoms\n", system.proto.PDBID, system.proto.name.c_str(), (int)system.proto.atoms.size());
        system.proto.printAll();
        //for (int i=0; i<system.proto.atoms.size(); i++) {
         //   system.proto.atoms[i].printAll();
        //}    


    } // end if MC
    
    // END MOLECULE PRINTOUT
}

void setupBox(System &system) {
    if (system.pbc.alpha == 90 && system.pbc.beta == 90 && system.pbc.gamma == 90) {
	system.pbc.x_max = system.pbc.x_length/2.0;
	system.pbc.x_min = -system.pbc.x_max;
	system.pbc.y_max = system.pbc.y_length/2.0;
	system.pbc.y_min = -system.pbc.y_max;
	system.pbc.z_max = system.pbc.z_length/2.0;
	system.pbc.z_min = -system.pbc.z_max;
    } 
    system.pbc.calcVolume();
    system.pbc.calcRecip();
    system.pbc.calcCutoff();
    system.constants.ewald_alpha = 3.5/system.pbc.cutoff;
    system.pbc.calcBoxVertices();
    system.pbc.calcPlanes();
    system.pbc.printBasis();
}

void setCheckpoint(System &system) {
    // saves variables of interest to temporary storage JIC

    // VALUES
        // ENERGY
    system.last.rd = system.stats.rd.value;
        system.last.lj_lrc = system.stats.lj_lrc.value;
        system.last.lj_self_lrc = system.stats.lj_self_lrc.value;
        system.last.lj = system.stats.lj.value;
    system.last.es = system.stats.es.value;
        system.last.es_self = system.stats.es_self.value;
        system.last.es_real = system.stats.es_real.value;
        system.last.es_recip = system.stats.es_recip.value;
    system.last.polar = system.stats.polar.value;
    system.last.potential = system.stats.potential.value;
        // VOLUME
    system.last.volume = system.stats.volume.value;
        // DENSITY
    system.last.density = system.stats.density.value;
        // CHEMICAL POTENTIAL
    system.last.chempot = system.stats.chempot.value;
        // Z
    system.last.z = system.stats.z.value;
        // QST
    system.last.qst = system.stats.qst.value;
        // Nmov
    system.last.Nmov = system.stats.Nmov.value;

}

void revertToCheckpoint(System &system) {
    // reverts variables of interest if needed
    
    // VALUES
        // ENERGY
    system.stats.rd.value = system.last.rd;
        system.stats.lj_lrc.value = system.last.lj_lrc;
        system.stats.lj_self_lrc.value = system.last.lj_self_lrc;
        system.stats.lj.value = system.last.lj;
    system.stats.es.value = system.last.es;
        system.stats.es_self.value = system.last.es_self;
        system.stats.es_real.value = system.last.es_real;
        system.stats.es_recip.value = system.last.es_recip;
    system.stats.polar.value = system.last.polar;
    system.stats.potential.value = system.last.potential;    
        // VOLUME
    system.stats.volume.value = system.last.volume;
        // DENSITY
    system.stats.density.value = system.last.density;
        // CHEMICAL POTENTIAL
    system.stats.chempot.value = system.last.chempot;
        // Z
    system.stats.z.value = system.last.z;
        // QST
    system.stats.qst.value = system.last.qst;    
        // Nmov
    system.stats.Nmov.value = system.last.Nmov;
    
}

void initialize(System &system) {
        system.stats.Nsq.name = "Nsq";
        system.stats.NU.name = "NU";
        system.stats.qst.name = "qst";
        system.stats.rd.name = "rd";
        system.stats.es.name = "es";
        system.stats.polar.name = "polar";
        system.stats.potential.name = "potential"; 
        system.stats.density.name = "density";
        system.stats.volume.name = "volume";
        system.stats.z.name = "z";
        system.stats.Nmov.name = "Nmov";
        system.stats.wtp.name = "wtp";
        system.stats.wtpME.name = "wtpME";
        system.stats.lj_lrc.name = "lj_lrc";
        system.stats.lj_self_lrc.name = "lj_self_lrc";
        system.stats.lj.name = "lj";
        system.stats.es_self.name = "es_self";
        system.stats.es_real.name = "es_real";
        system.stats.es_recip.name = "es_recip";
        system.stats.chempot.name = "chempot";
        system.stats.totalmass.name = "totalmass";
        system.stats.frozenmass.name = "frozenmass";
        system.stats.movablemass.name = "moveablemass";
        system.stats.pressure.name = "pressure";
        system.stats.temperature.name= "temperature";
}
