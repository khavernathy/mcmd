#include <iostream>
#include <string>
#ifdef WINDOWS
       #include <string.h>
#else
       #include <strings.h>
#endif
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
using namespace std;

double getrand() {
#ifndef WINDOWS
    return drand48();
#else
    return (std::rand()/(RAND_MAX+1));
#endif
}

/* GET BOX LIMIT COORDINATE FOR PBC */
double getBoxLimit(System &system, int planeid, double x, double y, double z) {

    // 0: -x; 1: +x; 2: -y; 3: +y; 4: -z; 5: +z;
    if (planeid == 0) return (-system.pbc.D[0] - system.pbc.B[0]*y - system.pbc.C[0]*z)/system.pbc.A[0]; // -x
    else if (planeid == 1) return (-system.pbc.D[1] - system.pbc.B[1]*y - system.pbc.C[1]*z)/system.pbc.A[1]; // +x
    else if (planeid == 2) return (-system.pbc.D[2] - system.pbc.A[2]*x - system.pbc.C[2]*z)/system.pbc.B[2]; // -y
    else if (planeid == 3) return (-system.pbc.D[3] - system.pbc.A[3]*x - system.pbc.C[3]*z)/system.pbc.B[3]; // +y
    else if (planeid == 4) return (-system.pbc.D[4] - system.pbc.A[4]*x - system.pbc.B[4]*y)/system.pbc.C[4]; // -z
    else if (planeid == 5) return (-system.pbc.D[5] - system.pbc.A[5]*x - system.pbc.B[5]*y)/system.pbc.C[5]; // +z
    else {
        printf("error: planeid is not valid.\n");
        std::exit(0);
    }
}


/* MOVE ALL ATOMS SUCH THAT THEY ARE CENTERED ABOUT 0,0,0 */
void centerCoordinates(System &system) {
	printf("Centering all coordinates...\n");
    	int size = system.constants.total_atoms;
	#ifndef WINDOWS
	double xtemp[size];
	double ytemp[size];
	double ztemp[size];
	#else
	double* xtemp = new double[size];
	double* ytemp = new double[size];
	double* ztemp = new double[size];
	#endif
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

    printf("Pre-centered values:\n");
    printf("----> xmax: %.5f ymax: %.5f zmax: %.5f\n", xmax, ymax, zmax);
    printf("----> xmin: %.5f ymin: %.5f zmin: %.5f\n", xmin, ymin, zmin);
    printf("----> xlen: %.5f ylen: %.5f zlen: %.5f\n", xmax-xmin, ymax-ymin, zmax-zmin);

	for (int i=0; i<system.molecules.size(); i++) {
	for (int j=0; j<system.molecules[i].atoms.size(); j++) {
		system.molecules[i].atoms[j].pos[0] = system.molecules[i].atoms[j].pos[0] - (xmin + (xmax - xmin)/2.0);
		system.molecules[i].atoms[j].pos[1] = system.molecules[i].atoms[j].pos[1] - (ymin + (ymax - ymin)/2.0);
		system.molecules[i].atoms[j].pos[2] = system.molecules[i].atoms[j].pos[2] - (zmin + (zmax - zmin)/2.0);
	}
	}
    printf("Done centering coordinates.\n\n");
}


/* CALCULATE CENTER OF MASS OF THE SYSTEM */
double * centerOfMass(System &system) {

    double x_mass_sum=0.0; double y_mass_sum=0.0; double z_mass_sum=0.0; double mass_sum=0.0;

	for (int j=0; j<system.molecules.size(); j++) {
        for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                double atom_mass = system.molecules[j].atoms[i].mass;
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
        
        system.molecules[i].calc_center_of_mass();

        if (system.molecules[i].com[0] > system.pbc.x_max) { // right of box
            system.molecules[i].diffusion_corr[0] += system.pbc.x_length;
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[0] -= system.pbc.x_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
        else if (system.molecules[i].com[0] < system.pbc.x_min) {
            system.molecules[i].diffusion_corr[0] -= system.pbc.x_length;
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
	                system.molecules[i].atoms[k].pos[0] += system.pbc.x_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
        if (system.molecules[i].com[1] > system.pbc.y_max) {
            system.molecules[i].diffusion_corr[1] += system.pbc.y_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] -= system.pbc.y_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
        else if (system.molecules[i].com[1] < system.pbc.y_min) {
            system.molecules[i].diffusion_corr[1] -= system.pbc.y_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] += system.pbc.y_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
        if (system.molecules[i].com[2]  > system.pbc.z_max) {
            system.molecules[i].diffusion_corr[2] += system.pbc.z_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] -= system.pbc.z_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
        else if (system.molecules[i].com[2] < system.pbc.z_min) {
            system.molecules[i].diffusion_corr[2] -= system.pbc.z_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] += system.pbc.z_length;
                }
                system.molecules[i].calc_center_of_mass();
        }
} // end if alpha=beta=gamma=90

// the universal treatment, alpha != beta ?= gamma
else {
    double box_limit[6]; //-x, +x, -y, +y, -z, +z limits, which are functions of atom position in other dims.
    double newlimit[6]; // an array to hold temporary new planes for movements.
    double tmp_com[3]; // temporary COM to track diffusion corrections

    // first refresh COM
    system.molecules[i].calc_center_of_mass();

    // find appropriate values for box limits based on atom coordinates.
    // check in this order: (-x, +x, -y, +y, -z, +z)
    for (int n=0; n<6; n++)
        box_limit[n] = getBoxLimit(system, n, system.molecules[i].com[0], system.molecules[i].com[1], system.molecules[i].com[2]);

        if (system.molecules[i].com[0] < box_limit[0]) { // left of box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
                system.molecules[i].calc_center_of_mass();
                for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n]; 
        }
        else if (system.molecules[i].com[0] > box_limit[1]) { // right of box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
            system.molecules[i].calc_center_of_mass();
            for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n];
        }
        if (system.molecules[i].com[1] < box_limit[2]) { // below box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
            system.molecules[i].calc_center_of_mass();
            for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n];
        }
        else if (system.molecules[i].com[1] > box_limit[3]) { // above box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
            system.molecules[i].calc_center_of_mass();
            for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n];

        }
        if (system.molecules[i].com[2] < box_limit[4]) { // behind box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
            system.molecules[i].calc_center_of_mass();
            for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n];
        }
        else if (system.molecules[i].com[2] > box_limit[5]) { // in front of box
            for (int n=0;n<3;n++) tmp_com[n] = system.molecules[i].com[n];
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
            system.molecules[i].calc_center_of_mass();
            for (int n=0;n<3;n++) system.molecules[i].diffusion_corr[n] -= system.molecules[i].com[n] - tmp_com[n];
        }
} // end if non-90/90/90
} // end pbc function


void addAtomToProto(System &system, int protoid, string name, string molname, string MF, double x, double y, double z, double mass, double charge, double polarizability, double epsilon, double sigma) {
    // initialize
    Atom atom;
    // PDBIDS
    int last_mol_index, last_mol_pdbid, last_atom_pdbid;
    if (system.molecules.size() > 0) {
        last_mol_index = system.molecules.size() - 1;
        last_mol_pdbid = system.molecules[last_mol_index].PDBID;
        last_atom_pdbid = system.molecules[last_mol_index].atoms[system.molecules[last_mol_index].atoms.size() - 1].PDBID;
        atom.PDBID = (int)(last_atom_pdbid + 1);
        atom.mol_PDBID = (int)(last_mol_pdbid + 1);
    } else {
        atom.PDBID = 1;
        atom.mol_PDBID = 1;
    }


    atom.name = name;
    atom.mol_name = molname;
    if (MF == "M") atom.frozen = 0;
    else if (MF == "F") atom.frozen = 1;
    atom.pos[0] = x;
    atom.pos[1] = y;
    atom.pos[2] = z;
    atom.mass = mass;// * system.constants.cM;
    atom.C = charge * system.constants.E2REDUCED;
    atom.polar = polarizability;
    atom.eps = epsilon;
    atom.sig = sigma;

    system.proto[protoid].mass += atom.mass;
    // push the atom into the protoype molecule
    system.proto[protoid].atoms.push_back(atom);

}

void moleculePrintout(System &system) {
    // CONFIRM ATOMS AND MOLECULES PRINTOUT

    // an excessive full list of atoms/molecules.
    /*
    for (int b=0; b<system.molecules.size(); b++) {
        printf("Molecule PDBID = %i: %s has %i atoms and frozen=%i; The first atom has PDBID = %i\n", system.molecules[b].PDBID, system.molecules[b].name.c_str(), (int)system.molecules[b].atoms.size(), system.molecules[b].frozen, system.molecules[b].atoms[0].PDBID);
        for (int c=0; c<system.molecules[b].atoms.size(); c++) {
            system.molecules[b].atoms[c].printAll();
            printf("\n");
        }

    }
    */

    if (system.constants.mode == "mc" || system.constants.mode == "md" || system.constants.mode == "sp") { // prototypes will be considered for MC and MD (for multisorb stuff)  
    // CHANGE THE PROTOTYPE IF USER SUPPLIED A KEYWORD IN INPUT
    // THIS WILL OVERWRITE ANY PROTOTYPE IN THE INPUT ATOMS FILE if user put it there, e.g. whatever.pdb
        if (system.constants.sorbate_name.size() > 0) {
            // clear the initial prototype
            if (system.proto.size() > 0) {
                system.proto[0].reInitialize();
            }
            else if (system.proto.size() == 0) {
                Molecule firstie;
                system.proto.push_back(firstie);
            }
        for (int i=0; i<system.constants.sorbate_name.size(); i++) {
            Molecule noobie;
            if (i>0) system.proto.push_back(noobie);
            int last_mol_index;
            int last_mol_pdbid=0;

            string sorbmodel = system.constants.sorbate_name[i];
            if (system.molecules.size() > 0) {
                last_mol_index = system.molecules.size() -1;
                last_mol_pdbid = system.molecules[last_mol_index].PDBID;
            }
            system.proto[i].PDBID = last_mol_pdbid+1;
            system.proto[i].frozen = 0;

            // each call to addAtomToProto takes 12 arguments which correspond to PDB-style input
            // proto id, atom name, molecule name, x, y, z, m, q, a, eps, sig
            // HYDROGEN H2
            if (sorbmodel == "h2_buch") {
                addAtomToProto(system, i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 2.016, 0.0, 0.0, 34.2, 2.96);
                system.proto[i].name = "H2";
                system.proto[i].dof = 3;
            }
            else if (sorbmodel == "h2_bss") {
                addAtomToProto(system, i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.74640, 0.0, 8.85160, 3.2293);
                addAtomToProto(system, i, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, i, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"H2N", "H2", "M", 0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                addAtomToProto(system, i,"H2N", "H2", "M", -0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                system.proto[i].name = "H2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "h2_dl") {
                addAtomToProto(system, i,"H2G", "H2", "M", 0.00, 0.0, 0.0,   0.0, -0.93600, 0.0, 36.7, 2.958);
                addAtomToProto(system, i,"H2E", "H2", "M", -0.370, 0.0, 0.0, 1.008, 0.468, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"H2E", "H2", "M", 0.370, 0.0, 0.0,  1.008, 0.468, 0.0, 0.0, 0.0);
                system.proto[i].name = "H2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "h2_bssp") {
                addAtomToProto(system,i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.7464, 0.6938, 12.76532, 3.15528);
                addAtomToProto(system,i, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.3732, 0.00044, 0.0, 0.0);
                addAtomToProto(system,i, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.00044, 0.0, 0.0);
                addAtomToProto(system,i, "H2N", "H2", "M", 0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                addAtomToProto(system,i, "H2N", "H2", "M", -0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                system.proto[i].name = "H2";
                system.proto[i].dof = 5;
            }
            // HELIUM He (useful for theor. calc. of pore vol.). See http://pubs.acs.org/doi/abs/10.1021/jp050948l
            // and http://onlinelibrary.wiley.com/doi/10.1002/aic.690470521/abstract
            // and our paper ESI http://pubs.rsc.org/en/Content/ArticleLanding/2014/CC/c4cc03070b#!divAbstract
            else if (sorbmodel == "he" || sorbmodel == "He") {
                addAtomToProto(system, i,"He", "He", "M", 0.0, 0.0, 0.0, 4.002602, 0.0, 0.0, 10.220, 2.280);
                system.proto[i].name = "He";
                system.proto[i].dof = 3;
            }
            else if (sorbmodel == "he_hogan") {
                addAtomToProto(system, i,"He", "He", "M", 0.0, 0.0, 0.0, 4.0026, 0.0, 0.2049407, 9.071224, 2.653089);
                system.proto[i].name = "He";
                system.proto[i].dof = 3;
            }
            // OTHER NOBLE GASES
            else if (sorbmodel == "ne_hogan") {
                addAtomToProto(system, i, "Ne", "Ne", "M", 0.0, 0.0, 0.0, 20.1797, 0.0, 0.3913212, 36.824138, 2.785823); 
                system.proto[i].name = "Ne";
                system.proto[i].dof = 3;
            }
            else if (sorbmodel == "ar_hogan") {
                addAtomToProto(system, i, "Ar", "Ar", "M", 0.0, 0.0, 0.0, 39.948, 0.0, 1.6392212, 128.326802, 3.371914);
                system.proto[i].name = "Ar";
                system.proto[i].dof = 3;
            } else if (sorbmodel == "kr_hogan") {
                addAtomToProto(system, i, "Kr", "Kr", "M", 0.0, 0.0, 0.0, 83.798, 0.0, 2.5004096, 183.795833, 3.601271);
               system.proto[i].name = "Kr"; 
                system.proto[i].dof = 3;
            } else if (sorbmodel == "xe_hogan") {
                addAtomToProto(system, i, "Xe", "Xe", "M", 0.0, 0.0, 0.0, 131.293, 0.0, 4.0232578, 237.985247, 3.956802);
                system.proto[i].name = "Xe"; 
                system.proto[i].dof = 3;
            }
                 // NITRIC OXIDE

            else if (sorbmodel == "no_3site") {
                addAtomToProto(system, i,"N", "NOG", "M", 0.0000, 0.0000, -0.6154, 14.007, -0.250, 0.000, 80.4839, 3.57);
                addAtomToProto(system, i,"O", "NOG", "M", 0.0000, 0.0000,  0.5356, 15.999, -0.345, 0.000, 100.6049, 3.66);
                addAtomToProto(system, i,"CoM", "NOG", "M", 0.0001, 0.0017, -0.0309, 0.0000,  0.595, 0.000, 0.00000, 0.000);
                system.proto[i].name = "NO";
                system.proto[i].dof = 5;
            }


            // CARBON MONOXIDE
             else if (sorbmodel == "co_3site") {
                addAtomToProto(system, i,"C", "COG", "M", -1.8980, -0.2855, 0.0000, 12.011, -0.750, 0.000, 13.1709, 3.830);
                addAtomToProto(system, i,"O", "COG", "M", -0.7735, -0.1963, 0.0000, 15.999, -0.850, 0.000, 80.0312, 3.120);
                addAtomToProto(system, i,"CoM", "COG", "M", -1.3488, -0.2402, 0.0000, 0.0000,  1.600, 0.000, 0.00000, 0.000);
                system.proto[i].name = "CO";
                system.proto[i].dof = 5;
            }


            // NIMESULIDE  DRUG
                                     
             else if (sorbmodel == "nimesulide2") {
                addAtomToProto(system, i,"O", "NIM", "M", -0.3020, -1.1081, -0.8817, 15.9998, -0.4756, 0.000, 30.18100, 3.12500);
                addAtomToProto(system, i,"C", "NIM", "M", -0.1556,  0.2232, -0.6141, 12.0110,  0.3424, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -2.2735, -2.3917, -1.2636, 12.0110, -0.3846, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -1.4206, -1.7658, -0.3693, 12.0110,  0.5541, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -1.2189,  1.0684, -0.3402, 12.0110, -0.4551, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -1.6209, -1.8231,  1.0040, 12.0110, -0.3595, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -2.7179, -2.5264,  1.4888, 12.0110, -0.0758, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -3.5896, -3.1646,  0.6083, 12.0110, -0.2133, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -3.3641, -3.0992, -0.7633, 12.0110, -0.0621, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M",  0.3462,  2.9371, -0.1611, 12.0110, -0.3012, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M",  1.3924,  2.0790, -0.4658, 12.0110, -0.2068, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"C", "NIM", "M", -0.9324,  2.4063, -0.1109, 12.0110,  0.1183, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"O", "NIM", "M", -1.7996,  4.4474,  0.5211, 15.9998, -0.4937, 0.000, 30.18100, 3.12500);
                addAtomToProto(system, i,"N", "NIM", "M", -2.0540,  3.3052,  0.1932, 14.0070,  0.8429, 0.000, 34.70800, 3.26700);
                addAtomToProto(system, i,"O", "NIM", "M", -3.1843,  2.8666,  0.1045, 15.9998, -0.4923, 0.000, 30.18100, 3.12500);
                addAtomToProto(system, i,"C", "NIM", "M",  1.1512,  0.7273, -0.6812, 12.0110,  0.2595, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"N", "NIM", "M",  2.2284, -0.1545, -0.9867, 14.0070, -0.6642, 0.000, 34.70800, 3.26700);
                addAtomToProto(system, i,"S", "NIM", "M",  3.1420, -0.6438,  0.3455, 32.0100,  1.2724, 0.000, 137.8290, 3.60000);
                addAtomToProto(system, i,"O", "NIM", "M",  2.3908, -1.5798,  1.1817, 15.9998, -0.6042, 0.000, 30.18100, 3.12500);
                addAtomToProto(system, i,"O", "NIM", "M",  3.6188,  0.5728,  0.9929, 15.9998, -0.5976, 0.000, 30.18100, 3.12500);
                addAtomToProto(system, i,"C", "NIM", "M",  4.4642, -1.4989, -0.4495, 12.0110, -0.5414, 0.000, 52.81700, 3.43700);
                addAtomToProto(system, i,"H", "NIM", "M", -2.0795, -2.3224, -2.3272, 1.00790,  0.2027, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M", -2.2337,  0.6970, -0.2920, 1.00790,  0.2613, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M", -0.9281, -1.3253,  1.6737, 1.00790,  0.2025, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M", -2.8876, -2.5796,  2.5576, 1.00790,  0.1527, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M", -4.4421, -3.7123,  0.9918, 1.00790,  0.1609, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M", -4.0394, -3.5944, -1.4510, 1.00790,  0.1503, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  0.5192,  3.9891,  0.0161, 1.00790,  0.2301, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  2.4049,  2.4530, -0.5428, 1.00790,  0.2301, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  1.9010, -0.9883, -1.4793, 1.00790,  0.3876, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  4.9884, -0.7965, -1.0944, 1.00790,  0.1698, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  5.1126, -1.8646,  0.3476, 1.00790,  0.1805, 0.000, 22.13300, 2.57600);
                addAtomToProto(system, i,"H", "NIM", "M",  4.0410, -2.3311, -1.0122, 1.00790,  0.2090, 0.000, 22.13300, 2.57600);
                system.proto[i].name = "NIM";
            }
                       // ASPIRIN DRUG 
             else if (sorbmodel == "aspirin") {
                addAtomToProto(system, i,"O", "ASA", "M", -1.1325, -0.7005, -0.7980, 15.9998, -0.75018, 0.00, 85.51420, 2.996986);
                addAtomToProto(system, i,"O", "ASA", "M",  1.0182,  2.8790, -0.2302, 15.9998, -0.83274, 0.00, 105.6352, 2.95836);
                addAtomToProto(system, i,"O", "ASA", "M", -1.0687,  1.9501,  0.1565, 15.9998, -0.71355, 0.00, 105.6352, 2.95836);
                addAtomToProto(system, i,"O", "ASA", "M", -1.8361, -0.9178,  1.4272, 15.9998, -0.74120, 0.00, 105.6352, 2.95836);
                addAtomToProto(system, i,"C", "ASA", "M",  0.2396, -0.7175, -0.3803, 12.0110,  0.60184, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  0.8969,  0.4933, -0.0524, 12.0110, -0.32781, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  0.9287, -1.9418, -0.4127, 12.0110, -0.30031, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  2.2765,  0.4467,  0.2817, 12.0110, -0.00665, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  2.2991, -1.9769, -0.0736, 12.0110, -0.00867, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  2.9725, -0.7797,  0.2675, 12.0110, -0.20353, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M",  0.1598,  1.7860, -0.0275, 12.0110,  1.08094, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M", -2.1231, -0.7140,  0.2361, 12.0110,  1.16628, 0.00, 43.26000, 3.39624);
                addAtomToProto(system, i,"C", "ASA", "M", -3.4795, -0.4451, -0.3531, 12.0110, -0.47488, 0.00, 55.03000, 3.39624);
                addAtomToProto(system, i,"H", "ASA", "M",  0.5029,  3.7182, -0.1631, 1.00790,  0.51145, 0.00, 0.00000, 0.0000);
                addAtomToProto(system, i,"H", "ASA", "M",  0.3848, -2.8482, -0.6817, 1.00790,  0.17586, 0.00, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "ASA", "M",  2.7903,  1.3740,  0.5330, 1.00790,  0.14167, 0.00, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "ASA", "M",  2.8373, -2.9256, -0.0887, 1.00790,  0.12072, 0.00, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "ASA", "M",  4.0309, -0.8034,  0.5307, 1.00790,  0.13166, 0.00, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "ASA", "M", -3.5558,  0.6255, -0.5873, 1.00790,  0.18078, 0.00, 7.897000, 2.6433);
                addAtomToProto(system, i,"H", "ASA", "M", -4.2472, -0.7180,  0.3785, 1.00790,  0.11701, 0.00, 7.897000, 2.6433);
                addAtomToProto(system, i,"H", "ASA", "M", -3.6127, -1.0160, -1.2797, 1.00790,  0.13130, 0.00, 7.897000, 2.6433);
                   system.proto[i].name = "ASA";
            }
                       // NIMESULIDE  DRUG
             else if (sorbmodel == "naproxen") {
                addAtomToProto(system, i,"C", "NAP", "M",  1.3197, -1.7554,  0.2851, 12.0110, -0.22975, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  2.7059, -1.7144,  0.1949, 12.0110, -0.25888, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  3.3717, -0.4640, -0.0275, 12.0110,  0.55480, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  2.6561,  0.7310, -0.0970, 12.0110, -0.45645, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  1.2173,  0.7020, -0.0386, 12.0110,  0.13412, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  0.5401, -0.5442,  0.2414, 12.0110,  0.19791, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M", -0.8908, -0.5670,  0.3020, 12.0110, -0.38815, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M",  0.4362,  1.9042, -0.0960, 12.0110, -0.15942, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M", -0.9573,  1.8566, -0.0105, 12.0110, -0.24515, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M", -1.6379,  0.6138,  0.2056, 12.0110,  0.16861, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M", -3.1670,  0.6057,  0.3607, 12.0110, -0.02976, 0.0, 55.0300, 1.9080);
                addAtomToProto(system, i,"C", "NAP", "M", -3.7627, -0.4506, -0.5493, 12.0110,  1.00663, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"O", "NAP", "M", -3.9194, -1.6631, -0.2784, 15.9998, -0.71289, 0.0, 105.635, 1.6612);
                addAtomToProto(system, i,"O", "NAP", "M", -4.1112,  0.0830, -1.7988, 15.9998, -0.90661, 0.0, 105.635, 1.6612);
                addAtomToProto(system, i,"C", "NAP", "M", -3.6055,  0.3855,  1.8322, 12.0110, -0.08785, 0.0, 55.0300, 1.9080);
                addAtomToProto(system, i,"O", "NAP", "M",  4.7781, -0.5684, -0.1012, 15.9998, -0.59226, 0.0, 85.5142, 1.6837);
                addAtomToProto(system, i,"C", "NAP", "M",  5.5451,  0.6790, -0.3050, 12.0110,  0.25781, 0.0, 43.2600, 1.9080);
                addAtomToProto(system, i,"H", "NAP", "M",  0.8109, -2.7086,  0.4534, 1.00790,  0.16859, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M",  3.3110, -2.6214,  0.2343, 1.00790,  0.15784, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M",  3.1497,  1.6875, -0.2807, 1.00790,  0.19471, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M", -1.3953, -1.5245,  0.4707, 1.00790,  0.19664, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M",  0.9396,  2.8597, -0.2730, 1.00790,  0.12683, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M", -1.5392,  2.7798, -0.0836, 1.00790,  0.14411, 0.0, 7.89700, 1.4590);
                addAtomToProto(system, i,"H", "NAP", "M", -3.5453,  1.5734, -0.0025, 1.00790,  0.07602, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M", -4.4424, -0.6308, -2.3958, 1.00790,  0.54182, 0.0, 0.00000, 0.0000);
                addAtomToProto(system, i,"H", "NAP", "M", -3.2677, -0.5977,  2.1864, 1.00790,  0.03295, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M", -3.1598,  1.1667,  2.4636, 1.00790,  0.01686, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M", -4.7007,  0.4347,  1.9231, 1.00790,  0.03004, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M",  6.5868,  0.3522, -0.3233, 1.00790,  0.04872, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M",  5.2726,  1.1481, -1.2604, 1.00790,  0.00460, 0.0, 7.89700, 1.4850);
                addAtomToProto(system, i,"H", "NAP", "M",  5.3743,  1.3748,  0.5278, 1.00790,  0.00755, 0.0, 7.89700, 1.4850);
                system.proto[i].name = "NAP";
            }
                         else if (sorbmodel == "floubiprofe") {
                addAtomToProto(system, i,"F", "FBP", "M", -0.9202,  2.2358, -0.0412, 18.9984, -0.282541, 0.0, 25.16000, 2.99700);
                addAtomToProto(system, i,"O", "FBP", "M",  4.5296,  0.3530,  1.3520, 15.9998, -0.640210, 0.0, 108.1503, 3.13300);
                addAtomToProto(system, i,"O", "FBP", "M",  3.7195, -1.7101,  1.0074, 15.9998, -0.820002, 0.0, 108.1503, 3.13300);
                addAtomToProto(system, i,"C", "FBP", "M",  3.3078, -0.1236, -0.6936, 12.0110, -0.098228, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.7894, -0.1064, -0.5095, 12.0110,  0.176176, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -1.0026, -0.1273, -0.1726, 12.0110, -0.233496, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.0863,  1.0733, -0.3585, 12.0110, -0.432004, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.0873, -1.3040, -0.4754, 12.0110, -0.313453, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -0.2784,  1.0456, -0.1928, 12.0110,  0.487284, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -0.2819, -1.3099, -0.3104, 12.0110, -0.025577, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  3.8910,  1.1772, -1.2560, 12.0110, -0.270875, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -2.4801, -0.1412, -0.0146, 12.0110,  0.174684, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  3.9286, -0.4166,  0.6521, 12.0110,  0.969497, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -3.1128,  0.6643,  0.9265, 12.0110, -0.203481, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -3.2539, -0.9881, -0.7996, 12.0110, -0.184674, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -4.4840,  0.6072,  1.0886, 12.0110, -0.095761, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -4.6266, -1.0418, -0.6380, 12.0110, -0.115052, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -5.2455, -0.2469,  0.3096, 12.0110, -0.156147, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"H", "FBP", "M",  3.8032,  1.9787, -0.5366, 1.00790,  0.089195, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  4.9448,  1.0464, -1.4686, 1.00790,  0.087171, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  3.3859,  1.4543, -2.1719, 1.00790,  0.078532, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  4.0618, -1.9042,  1.8925, 1.00790,  0.506691, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -6.3090, -0.2877,  0.4362, 1.00790,  0.134234, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -5.2093, -1.6967, -1.2549, 1.00790,  0.131042, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -2.5316,  1.3362,  1.5212, 1.00790,  0.156844, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -2.7827, -1.5922, -1.5492, 1.00790,  0.132650, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -4.9568,  1.2297,  1.8220, 1.00790,  0.126193, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  1.5698,  2.0256, -0.3702, 1.00790,  0.225729, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  3.5509, -0.9554, -1.3412, 1.00790,  0.089026, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -0.8062, -2.2429, -0.2685, 1.00790,  0.119426, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  1.6210, -2.2286, -0.5494, 1.00790,  0.187127, 0.0, 7.545000, 2.50802);
                   system.proto[i].name = "FBP";
            }
                         else if (sorbmodel == "floubiprofen") {
                addAtomToProto(system, i,"F", "FBP", "M", -0.9202,  2.2358, -0.0412, 18.9984, -0.282541, 0.0, 25.16000, 2.99700);
                addAtomToProto(system, i,"O", "FBP", "M",  4.5296,  0.3530,  1.3520, 15.9998, -0.640210, 0.0, 108.1503, 3.13300);
                addAtomToProto(system, i,"O", "FBP", "M",  3.7195, -1.7101,  1.0074, 15.9998, -0.820002, 0.0, 108.1503, 3.13300);
                addAtomToProto(system, i,"C", "FBP", "M",  3.3078, -0.1236, -0.6936, 12.0110, -0.098228, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.7894, -0.1064, -0.5095, 12.0110,  0.176176, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -1.0026, -0.1273, -0.1726, 12.0110, -0.233496, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.0863,  1.0733, -0.3585, 12.0110, -0.432004, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  1.0873, -1.3040, -0.4754, 12.0110, -0.313453, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -0.2784,  1.0456, -0.1928, 12.0110,  0.487284, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -0.2819, -1.3099, -0.3104, 12.0110, -0.025577, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  3.8910,  1.1772, -1.2560, 12.0110, -0.270875, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -2.4801, -0.1412, -0.0146, 12.0110,  0.174684, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M",  3.9286, -0.4166,  0.6521, 12.0110,  0.969497, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -3.1128,  0.6643,  0.9265, 12.0110, -0.203481, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -3.2539, -0.9881, -0.7996, 12.0110, -0.184674, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -4.4840,  0.6072,  1.0886, 12.0110, -0.095761, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -4.6266, -1.0418, -0.6380, 12.0110, -0.115052, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"C", "FBP", "M", -5.2455, -0.2469,  0.3096, 12.0110, -0.156147, 0.0, 47.78700, 3.46000);
                addAtomToProto(system, i,"H", "FBP", "M",  3.8032,  1.9787, -0.5366, 1.00790,  0.089195, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  4.9448,  1.0464, -1.4686, 1.00790,  0.087171, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  3.3859,  1.4543, -2.1719, 1.00790,  0.078532, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  4.0618, -1.9042,  1.8925, 1.00790,  0.506691, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -6.3090, -0.2877,  0.4362, 1.00790,  0.134234, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -5.2093, -1.6967, -1.2549, 1.00790,  0.131042, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -2.5316,  1.3362,  1.5212, 1.00790,  0.156844, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -2.7827, -1.5922, -1.5492, 1.00790,  0.132650, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -4.9568,  1.2297,  1.8220, 1.00790,  0.126193, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  1.5698,  2.0256, -0.3702, 1.00790,  0.225729, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  3.5509, -0.9554, -1.3412, 1.00790,  0.089026, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M", -0.8062, -2.2429, -0.2685, 1.00790,  0.119426, 0.0, 7.545000, 2.50802);
                addAtomToProto(system, i,"H", "FBP", "M",  1.6210, -2.2286, -0.5494, 1.00790,  0.187127, 0.0, 7.545000, 2.50802);
                   system.proto[i].name = "FBP";
            }
            // DMA cations   6-31+g
            else if (sorbmodel == "dma") {
                addAtomToProto(system, i,"C", "DMA", "M",  1.2542, -0.2761,  0.0000, 12.0110,  0.565846, 0.000, 47.85600, 3.470000);
                addAtomToProto(system, i,"H", "DMA", "M",  2.0958,  0.4201, -0.0000, 1.00790, -0.183826, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"H", "DMA", "M",  1.2678, -0.8919, -0.8992, 1.00790, -0.141727, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"H", "DMA", "M",  1.2678, -0.8919,  0.8992, 1.00790, -0.141884, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"C", "DMA", "M", -1.2542, -0.2761,  0.0000, 12.0110,  0.566114, 0.000, 47.85600, 3.470000);
                addAtomToProto(system, i,"H", "DMA", "M", -2.0958,  0.4201, -0.0000, 1.00790, -0.183894, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"H", "DMA", "M", -1.2678, -0.8919, -0.8992, 1.00790, -0.141876, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"H", "DMA", "M", -1.2678, -0.8919,  0.8992, 1.00790, -0.141898, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"N", "DMA", "M",  0.0000,  0.5293, -0.0000, 14.0070,  0.587023, 0.000, 38.93400, 3.269000);
                addAtomToProto(system, i,"H", "DMA", "M",  0.0000,  1.1676,  0.8330, 1.00790, -0.391833, 0.000, 7.649000, 2.850000);
                addAtomToProto(system, i,"H", "DMA", "M",  0.0000,  1.1676, -0.8330, 1.00790, -0.392046, 0.000, 7.649000, 2.850000);
                system.proto[i].name = "DMA";
                system.proto[i].dof = 6;
            }
              // DMA cations   aug-ccpvdz
            else if (sorbmodel == "dma2") {
                addAtomToProto(system, i,"C", "DMA", "M",  1.2542, -0.2761,  0.0000, 12.0110,  0.895648, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "DMA", "M",  2.0958,  0.4201, -0.0000, 1.00790, -0.314477, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "DMA", "M",  1.2678, -0.8919, -0.8992, 1.00790, -0.294736, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "DMA", "M",  1.2678, -0.8919,  0.8992, 1.00790, -0.294736, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"C", "DMA", "M", -1.2542, -0.2761,  0.0000, 12.0110,  0.895648, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "DMA", "M", -2.0958,  0.4201, -0.0000, 1.00790, -0.314477, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "DMA", "M", -1.2678, -0.8919, -0.8992, 1.00790, -0.294736, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "DMA", "M", -1.2678, -0.8919,  0.8992, 1.00790, -0.294736, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"N", "DMA", "M",  0.0000,  0.5293, -0.0000, 14.0070,  1.476108, 0.000, 85.51, 1.875000);
                addAtomToProto(system, i,"H", "DMA", "M",  0.0000,  1.1676,  0.8330, 1.00790, -0.729753, 0.000, 7.897, 0.600000);
                addAtomToProto(system, i,"H", "DMA", "M",  0.0000,  1.1676, -0.8330, 1.00790, -0.729753, 0.000, 7.897, 0.600000);
                system.proto[i].name = "DMA";
            }
            // TMA cations
            else if (sorbmodel == "tma") {
                addAtomToProto(system, i,"C", "TMA", "M", -1.2307, -0.8566, -0.0083, 12.0110,  0.28207, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "TMA", "M", -2.1026, -0.2021, -0.0032, 1.00790, -0.19171, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M", -1.2161, -1.4845,  0.8828, 1.00790, -0.19252, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M", -1.2149, -1.4688, -0.9103, 1.00790, -0.19342, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"C", "TMA", "M",  1.2178, -0.8747, -0.0069, 12.0110,  0.31334, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "TMA", "M",  2.0994, -0.2333, -0.0006, 1.00790, -0.20002, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M",  1.1928, -1.5024,  0.8842, 1.00790, -0.20064, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M",  1.1940, -1.4867, -0.9088, 1.00790, -0.20325, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"N", "TMA", "M",  0.0000,  0.0000,  0.0000, 14.0070,  1.15473, 0.000, 85.51, 1.875000);
                addAtomToProto(system, i,"C", "TMA", "M",  0.0056,  0.8549,  1.2318, 12.0110,  0.33304, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "TMA", "M",  0.9068,  1.4684,  1.2180, 1.00790, -0.20982, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M", -0.8863,  1.4817,  1.2169, 1.00790, -0.20663, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M",  0.0002,  0.1993,  2.1029, 1.00790, -0.20425, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"C", "TMA", "M",  0.0072,  0.8764, -1.2167, 12.0110,  0.35447, 0.000, 55.03, 1.908000);
                addAtomToProto(system, i,"H", "TMA", "M",  0.9084,  1.4895, -1.1909, 1.00790, -0.21088, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M",  0.0030,  0.2361, -2.0991, 1.00790, -0.21084, 0.000, 7.897, 1.100000);
                addAtomToProto(system, i,"H", "TMA", "M", -0.8847,  1.5028, -1.1919, 1.00790, -0.21368, 0.000, 7.897, 1.100000);
                system.proto[i].name = "TMA";
            }
            else if (sorbmodel == "co2_phast*" || sorbmodel == "co2_phastp") {
                addAtomToProto(system, i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77134, 1.22810, 19.61757, 3.03366);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "CON", "CO2", "M", 1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                system.proto[i].name = "CO2" ;
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "co2_phastq*" || sorbmodel == "co2_phastqp") {
                addAtomToProto(system, i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.66134, 1.22810, 26.89402, 3.18054);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.33067, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.33067, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "CON", "CO2", "M", 1.187, 0.0, 0.0, 0.0, 0.0, 0.0, 70.24356, 2.75458);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.187, 0.0, 0.0, 0.0, 0.0, 0.0, 70.24356, 2.75458);
                system.proto[i].name = "CO2" ;
                system.proto[i].dof = 5;


            }
            else if (sorbmodel == "co2_trappe") {
                addAtomToProto(system,i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.7, 0.0, 27.0, 2.80);
                addAtomToProto(system,i, "COE", "CO2", "M", 1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05);
                addAtomToProto(system,i, "COE", "CO2", "M", -1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05);
                system.proto[i].name = "CO2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "co2_becker") {
                // A polar model adapted from TraPPE. Becker et al. 10.1021/acs.jpcc.6b12052
                addAtomToProto(system,i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.7, 0.916, 23.4, 2.8);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.16, 0.0, 0.0, 16.0, -0.35, 0.575, 73.08, 3.05);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.16, 0.0, 0.0, 16.0, -0.35, 0.575, 73.08, 3.05);
                system.proto[i].name = "CO2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "co2_epm2") {
                // Harris and Yung, 1995, J. Phys. Chem.
                // https://pubs.acs.org/doi/pdf/10.1021/j100031a034
                addAtomToProto(system,i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.6512, 0.00, 28.129, 2.757);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.149, 0.0, 0.0, 16.0, -0.3256, 0.00, 80.507, 3.033);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.149, 0.0, 0.0, 16.0, -0.3256, 0.00, 80.507, 3.033);
                system.proto[i].name = "CO2";
                system.proto[i].dof = 5;
            }

            // NITROGEN N2
            else if (sorbmodel == "n2_mcquarrie") {
                addAtomToProto(system,i, "N2G", "N2", "M", 0.0, 0.0, 0.0, 28.01344, 0.0, 0.0, 95.1, 3.7);
                system.proto[i].name = "N2";
                system.proto[i].dof = 3;
            } else if (sorbmodel == "n2_trappe") {
                addAtomToProto(system,i, "N2G", "N2", "M", 0.0, 0.0, 0.0, 0.0, 0.964, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "N2E", "N2", "M", 0.550, 0.0, 0.0, 14.0067, -0.482, 0.0, 36.0, 3.31);
                addAtomToProto(system,i, "N2E", "N2", "M", -0.550, 0.0, 0.0, 14.0067, -0.482, 0.0, 36.0, 3.31);
                system.proto[i].name = "N2";
                system.proto[i].dof = 5;
            } else if (sorbmodel == "n2_nonpolar") {
                addAtomToProto(system,i,"N2G", "N2", "M", 0.000,   0.000,   0.000,  0.00000,  1.04742,  0.00000, 17.60293,  3.44522);
                addAtomToProto(system,i,"N2E", "N2", "M", 0.549,   0.000,   0.000, 14.00670, -0.52371,  0.00000,  0.00000,  0.00000);
                addAtomToProto(system,i,"N2E", "N2", "M", -0.549,   0.000,   0.000, 14.00670, -0.52371,  0.00000,  0.00000,  0.00000);
                addAtomToProto(system,i,"N2N", "N2", "M", 0.738,   0.000,   0.000,  0.00000,  0.00000,  0.00000, 18.12772,  3.15125);
                addAtomToProto(system,i,"N2N", "N2", "M", -0.738,   0.000,   0.000,  0.00000,  0.00000,  0.00000, 18.12772,  3.15125);
                system.proto[i].name = "N2";
                system.proto[i].dof = 5;
            } else if (sorbmodel == "n2_polar") {
                addAtomToProto(system,i,"N2G", "N2", "M", 0.192, 1.574, -4.563,  0.00000, 1.04742, 1.45590, 20.63650,  3.42344);
                addAtomToProto(system,i,"N2E", "N2", "M", -0.346, 1.624, -4.662, 14.00670, -0.52371,  0.51380,  0.00000,  0.00000);
                addAtomToProto(system,i,"N2E", "N2", "M", 0.730, 1.524,  -4.463, 14.00670, -0.52371,  0.51380,  0.00000,  0.00000);
                addAtomToProto(system,i,"N2N", "N2", "M", -0.537, 1.642, -4.698,  0.00000, 0.00000,  0.00000, 16.14200,  3.16141);
                addAtomToProto(system,i,"N2N", "N2", "M", 0.921, 1.506,  -4.428,  0.00000,  0.00000,  0.00000, 16.14200,  3.16141);
                system.proto[i].name = "N2";
                system.proto[i].dof = 5;
            }

            // METHANE CH4
            else if (sorbmodel == "ch4_trappe") {
                addAtomToProto(system,i, "CHG", "CH4", "M", 0.0, 0.0, 0.0, 16.0426, 0.0, 0.0, 148.0, 3.73);
                system.proto[i].name = "CH4";
                system.proto[i].dof = 3;
            }
            else if (sorbmodel == "ch4_9site") {
                addAtomToProto(system,i, "CHG", "CH4", "M", 0.0, 0.0, 0.0, 12.011, -0.5868, 0.0, 58.53869, 2.22416);
                addAtomToProto(system,i, "CHE", "CH4", "M", 0.0, 0.0, 1.099, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "CHE", "CH4", "M", 1.036, 0.0, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "CHE", "CH4", "M", -0.518, -0.897, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "CHE", "CH4", "M", -0.518, 0.897, -0.366, 1.0079, 0.14670, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "MOV", "CH4", "M", 0.0, 0.0, 0.816, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system,i, "MOV", "CH4", "M", 0.769, 0.0, -0.271, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system,i, "MOV", "CH4", "M", -0.385, -0.668, -0.271, 0.0, 0.0, 0.0, 16.85422, 2.96286);
                addAtomToProto(system,i, "MOV", "CH4", "M", -0.385,   0.668,  -0.271,  0.00000,  0.00000,  0.00000, 16.85422,  2.96286);
                system.proto[i].name = "CH4";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "ch4_9site*") {
                addAtomToProto(system,i, "CHG",  "CH4", "M", 0.000,  -0.000,  -0.000, 12.01100, -0.58680,  1.09870, 45.09730,  2.16247); //  0.00000  0.00000
                addAtomToProto(system,i, "CHE", "CH4", "M", 0.000,  0.000,1.099,  1.00790,  0.14670,  0.42460,  0.00000,  0.0000);
                addAtomToProto(system,i, "CHE",  "CH4", "M",  1.036,  -0.000,  -0.366,  1.00790,  0.14670,  0.42460, 0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system,i, "CHE",  "CH4", "M", -0.518, -0.897, -0.366,  1.00790,  0.14670,  0.42460,  0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system,i, "CHE",  "CH4", "M", -0.518,   0.897,  -0.366,  1.00790,  0.14670,  0.42460,  0.00000,  0.00000); //  0.00000  0.00000
                addAtomToProto(system,i, "MOV", "CH4", "M", 0.000,  -0.000,   0.814,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                addAtomToProto(system,i, "MOV", "CH4", "M", 0.768,  -0.000,  -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787);//  0.00000  0.00000
                addAtomToProto(system,i, "MOV",  "CH4", "M", -0.383,  -0.666, -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                addAtomToProto(system,i, "MOV",  "CH4", "M", -0.383,   0.666,  -0.270,  0.00000,  0.00000,  0.00000, 18.57167,  2.94787); //  0.00000  0.00000
                system.proto[i].name = "CH4";
                system.proto[i].dof = 6;
            }
            // ACETYLENE C2H2
            else if (sorbmodel == "c2h2" || sorbmodel == "acetylene") {
                addAtomToProto(system,i, "CoM", "ACE", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system,i, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system,i, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                addAtomToProto(system,i, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                system.proto[i].name = "ACE";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "c2h2*" || sorbmodel == "acetylene*") {
                addAtomToProto(system,i, "CoM", "ACE", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system,i, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system,i, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                addAtomToProto(system,i, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                system.proto[i].name = "ACE";
                system.proto[i].dof = 5;
            }
            // ETHLYENE C2H4
            else if (sorbmodel == "c2h4" || sorbmodel == "ethlyene" || sorbmodel == "ethene") {
                addAtomToProto(system,i, "CoM", "ETH", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "c2h4*" || sorbmodel == "ethlyene*" || sorbmodel == "ethene*") {
                addAtomToProto(system,i, "CoM", "ETH", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 7;
            }
            else if (sorbmodel == "c2h4_trappe") {
                addAtomToProto(system,i, "CH2",  "ETH",  "M", -0.665,   0.000,   0.000, 14.02200, -0.000,  0.00000, 85.000,  3.675);
                addAtomToProto(system,i, "CH2",  "ETH",  "M", 0.665,   0.000,   0.000, 14.02200, -0.000,  0.00000, 85.000,  3.675);
                system.proto[i].dof = 5;
            }

            // ETHANE C2H6
            else if (sorbmodel == "c2h6" || sorbmodel == "ethane") {
                addAtomToProto(system,i, "CoM", "ETH", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897); //  0.00000  0.00000
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "c2h6*" || sorbmodel == "ethane*") {
                addAtomToProto(system,i, "CoM", "ETH", "M", 0,0,0, 0, 0, 0, 0, 0);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.762,0.000, 0.000,12.01100, -0.04722,0.6967,98.63326,3.37151);
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.762,0.000,0.000,12.01100,-0.04722,0.6967,98.63326,  3.37151);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "c2h6_trappe") {
                addAtomToProto(system, i, "C2G",  "ETH",  "M", 0.770,   0.000,   0.000,  15.035, 0.0000,  0.00000, 98.00,  3.750);
                addAtomToProto(system, i, "C2G",  "ETH",  "M", -0.770,   0.000,   0.000,  15.035, 0.0000,  0.00000, 98.00,  3.750);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 5;
            }
            // PROPENE C3H6 -- Franz made w copt dft aug-cc-pvtz chelpg charge fit to ESP. 
            else if (sorbmodel == "c3h6" || sorbmodel == "propene" || sorbmodel == "propylene") {
                addAtomToProto(system,i,"C","PRO","M", -0.316,   1.089,   0.071, 12.01100, -0.63057, 1.288599, 52.837986, 3.430851);
                addAtomToProto(system,i,"C","PRO","M",0.332,  -0.079,   0.154, 12.01100, 0.16400, 1.288599, 52.837986, 3.430851);
                addAtomToProto(system,i,"H","PRO","M",0.213 ,  2.028 ,  0.200, 1.00790, 0.20492, 0.413835, 22.141632, 2.571134);
                addAtomToProto(system,i,"H","PRO","M",-1.382 ,  1.137,  -0.127, 1.00790, 0.22965, 0.413835, 22.141632, 2.571134);
                addAtomToProto(system,i,"C","PRO","M",-0.319,  -1.413,  -0.008, 12.01100, -0.40470, 1.288599, 52.837986, 3.430851);
                addAtomToProto(system,i,"H","PRO","M",1.402,  -0.078,   0.353, 1.00790, 0.06821, 0.413835, 22.141632, 2.571134);
                addAtomToProto(system,i,"H","PRO","M",-1.402,  -1.330,  -0.145, 1.00790, 0.13160, 0.413835, 22.141632 ,2.571134);
                addAtomToProto(system,i,"H","PRO","M",-0.135 , -2.028  , 0.878, 1.00790, 0.11707, 0.413835, 22.141632, 2.571134);
                addAtomToProto(system,i,"H","PRO","M",0.099 , -1.929,  -0.878, 1.00790, 0.11983, 0.413835, 22.141632, 2.571134);
                system.proto[i].name = "PRO";
                system.proto[i].dof = 6;
            }
            // WATER H2O (TIP3P)
            else if (sorbmodel == "h2o" || sorbmodel == "water" || sorbmodel == "tip3p" || sorbmodel == "h2o_tip3p") {
                addAtomToProto(system,i, "OXY", "H2O", "M", 0.0, 0.0, 0.0, 16.0, -0.834, 0.0, 76.42, 3.151);
                addAtomToProto(system,i, "HYD", "H2O", "M", -0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "HYD", "H2O", "M", 0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                system.proto[i].name = "H2O";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "tip4p" || sorbmodel == "water_tip4p" || sorbmodel == "h2o_tip4p") {
                addAtomToProto(system,i, "OXY", "H2O", "M", 0.0, 0.0, 0.0, 16.0, 0.0, 0.0, 78.0, 3.154);
                addAtomToProto(system,i, "HYD", "H2O","M",  0.58588, 0.75695, 0.0, 1.008, 0.52, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "HYD", "H2O","M",  0.58588, -0.75695, 0.0, 1.008, 0.52, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "M", "H2O", "M", 0.15, 0.0, 0.0, 0.0, -1.04, 0.0, 0.0, 0.0);
                system.proto[i].name = "H2O";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "h2o_roney") { // Ben Roney's polarizable H2O model

                addAtomToProto(system,i,"OXY",  "H2O", "M", 0.000 ,  0.000 ,  0.000, 15.99900, -0.66900,  0.83700, 78.22000,  3.16600); 
                addAtomToProto(system,i,"HYD","H2O","M",-0.761 , -0.588,   0.000,  1.00790,  0.33450,  0.49500,  0.00000,  0.00000);  
                addAtomToProto(system,i,"HYD",  "H2O", "M", 0.761 , -0.588 ,  0.000 , 1.00790,  0.33450,  0.49500,  0.00000,  0.00000);
                system.proto[i].name = "H2O";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "h2o_pol3") {
                addAtomToProto(system,i,"OXY",  "H2O", "M", 0.000,   0.000,   0.000, 15.99900, -0.73000,  0.52800, 78.50225,  3.59600); 
                addAtomToProto(system,i,"HYD","H2O", "M",-0.816 , -0.577,   0.000,  1.00790,  0.36500,  0.17000,  0.00000,  0.00000);
                addAtomToProto(system,i,"HYD","H2O", "M",   0.816,  -0.577,   0.000,  1.00790,  0.36500,  0.17000,  0.00000,  0.00000);
                system.proto[i].name = "H2O";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "h2o_franz") {
                addAtomToProto(system,i,"OXY",  "H2O", "M", 0.000,   0.000,   0.000, 15.99900, -0.73000,  0.52800, 78.50225,  3.59600); 
                addAtomToProto(system,i,"HYD","H2O", "M",-0.816 , -0.577,   0.000,  1.00790,  0.36500,  0.17000,  22.14, 2.571);
                addAtomToProto(system,i,"HYD","H2O", "M",   0.816,  -0.577,   0.000,  1.00790,  0.36500,  0.17000,  22.14, 2.571);
                system.proto[i].name = "H2O";
                system.proto[i].dof = 6;
            }

            // SULFUR DIOXIDE -- Peng et al. http://pubs.acs.org/doi/pdf/10.1021/acs.jpcc.7b01925 (Table 1)
            else if (sorbmodel == "so2") {
                addAtomToProto(system,i,"SUL", "SO2", "M", 0.00000, 0.00000, 0.0, 32.060, 0.4700, 2.4744476, 154.4, 3.585);
                addAtomToProto(system,i,"OXY", "SO2", "M", -1.4321, 0.00000, 0.0, 15.999, -0.235, 0.852, 62.3, 2.993);
                addAtomToProto(system,i,"OXY", "SO2", "M", 0.70520, 1.24636, 0.0, 15.999, -0.235, 0.852, 62.3, 2.993);
                system.proto[i].name = "SO2";
                system.proto[i].dof = 6;
            }
            // METHANOL CH3OH -- my model: UFF sig/eps; charges from dft aug-cc-pvtz on CCSDT=FULL aug-cc-pvtz geometry (CCCBDB); van Deuynen polariz's
            else if (sorbmodel == "methanol" || sorbmodel == "ch3oh" || sorbmodel == "ch4o") {
                addAtomToProto(system,i, "CME", "MET", "M", -0.046520, 0.662923, 0.0, 12.0107, 0.074830, 1.2886, 52.84, 3.431);
                addAtomToProto(system,i, "OME", "MET", "M", -0.046520, -0.754865, 0.0, 15.9994, -0.577001, 0.852, 30.19, 3.118);
                addAtomToProto(system,i, "HME", "MET", "M", -1.086272, 0.976083, 0.0, 1.00794, 0.085798, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HME", "MET", "M", 0.437877, 1.070546, 0.888953, 1.00794, 0.015921, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HME", "MET", "M", 0.437877, 1.070546, -0.888953, 1.00794, 0.016558, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HME", "MET", "M", 0.861794, -1.055796, 0.0, 1.00794, 0.383895, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "CoM", "MET", "M", -0.020179, -0.063588, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                system.proto[i].name = "MET";
                system.proto[i].dof = 6;
            }
            // ETHANOL C2H5OH -- my model: UFF sig/eps; charges from dft aug-cc-pvtz pbe0 on CCSDT aug-cc-pvtz geomtry (CCCBDB); van Deuynen polariz's

            else if (sorbmodel == "ethanol" || sorbmodel == "c2h5oh" || sorbmodel == "c2h6o") {
                addAtomToProto(system,i, "CET", "ETH", "M", 1.161583, -0.406755, 0.0, 12.0107, -0.352167, 1.2886, 52.84, 3.431);
                addAtomToProto(system,i, "CET", "ETH", "M", 0.0, 0.552718, 0.0, 12.0107, 0.362477, 1.2886, 52.84, 3.431);
                addAtomToProto(system,i, "OET", "ETH", "M", -1.187114, -0.212860, 0.0, 15.9994, -0.619551, 0.852, 30.19, 3.118);
                addAtomToProto(system,i, "HET", "ETH", "M", -1.932434, 0.383817, 0.0, 1.00794, 0.367986, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HET", "ETH", "M", 2.102860, 0.135840, 0.0, 1.00794, 0.072973, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HET", "ETH", "M", 1.122347, -1.039829, 0.881134, 1.00794, 0.115255, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HET", "ETH", "M", 1.122347, -1.039829, -0.881134, 1.00794, 0.114702, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HET", "ETH", "M", 0.056147, 1.193553, 0.880896, 1.00794, -0.031564, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "HET", "ETH", "M", 0.056147, 1.193553, -0.880896, 1.00794, -0.030110, 0.41380, 22.14, 2.571);
                addAtomToProto(system,i, "CoM", "ETH", "M", -0.054141, -0.017774, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                system.proto[i].name = "ETH";
                system.proto[i].dof = 6;
            }

            // BENZENE
            else if (sorbmodel == "benzene" || sorbmodel == "c6h6") {
                double wp = 1.4*sqrt(3)/2.0;
                double s=3.695, e=50.5, m=13.01833;
                addAtomToProto(system,i,"CH","BNZ","M", 1.4, 0, 0, m, 0, 0,e,s);
                addAtomToProto(system,i,"CH", "BNZ", "M", 0.7,wp,0, m, 0,0,e,s);
                addAtomToProto(system,i,"CH", "BNZ", "M", -0.7, wp, 0, m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",-1.4,0,0,m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",-0.7,-wp,0,m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",0.7,-wp,0,m,0,0,e,s); 
                system.proto[i].name = "BNZ";
                system.proto[i].dof = 6;
            }
            // TraPPE benzene 9-site
            else if (sorbmodel == "c6h6_9site" || sorbmodel == "benzene_9site") {
                double wp = 1.4*sqrt(3)/2.0;
                double s=3.74, e=48.0, m=13.01833;
                addAtomToProto(system,i,"CH","BNZ","M", 1.4, 0, 0, m, 0, 0,e,s);
                addAtomToProto(system,i,"CH", "BNZ", "M", 0.7,wp,0, m, 0,0,e,s);
                addAtomToProto(system,i,"CH", "BNZ", "M", -0.7, wp, 0, m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",-1.4,0,0,m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",-0.7,-wp,0,m,0,0,e,s);
                addAtomToProto(system,i,"CH","BNZ","M",0.7,-wp,0,m,0,0,e,s); 
                addAtomToProto(system,i,"CoM","BNZ", "M", 0,0,0,0,2.42,0,0,0);
                addAtomToProto(system,i,"PI","BNZ","M",0.785,0,0,0,-1.21,0,0,0);
                addAtomToProto(system,i,"PI","BNZ","M",-0.785,0,0,0,-1.21,0,0,0);
                system.proto[i].name = "BNZ";
                system.proto[i].dof = 6;
            }
            // CYCLOHEXANE: hybrid of TraPPE model and CCCBDB
            // coordinates are from CCD cc-pvDZ geometry CCCBDB, rest is TraPPE params.
            else if (sorbmodel == "cyclohexane" || sorbmodel == "c6h12") {
                double s=3.91, e=52.5, m=84.16/6;
                addAtomToProto(system,i,"CH2","CHX","M",-1.2221530,  -0.6603440,  0.3930520, m, 0, 0,e,s);
                addAtomToProto(system,i,"CH2","CHX","M",0.0001370 ,  1.5281720 ,  -0.0001270, m, 0, 0, e,s);
                addAtomToProto(system,i,"CH2","CHX","M",1.2221530,   -0.6604920 , -0.3929090, m, 0,0,e,s);
                addAtomToProto(system,i,"CH2","CHX","M",1.2221530  , 0.6603440 ,  0.3930520, m, 0,0,e,s);
                addAtomToProto(system,i,"CH2","CHX","M", -0.0001370,  -1.5281720,  -0.0001270,m,0,0,e,s);
                addAtomToProto(system,i,"CH2","CHX","M", -1.2221530 , 0.6604920 ,  -0.3929090,m,0,0,e,s);
                system.proto[i].name = "CHX"; 
                system.proto[i].dof = 6;
            }
        
            // OXYGEN : TraPPE : using p-table mass
            else if (sorbmodel == "o2" || sorbmodel == "oxygen") {
                double s=3.02, e=49.0, m=15.9997;
                addAtomToProto(system,i,"O","O2", "M", -0.605,0,0,m,-0.113,0,e,s);
                addAtomToProto(system,i,"CoM","O2","M", 0,0,0,0,0.226,0,0,0);
                addAtomToProto(system,i,"O","O2", "M", 0.605, 0,0,m,-0.113,0,e,s);
                system.proto[i].name = "O2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "o2_1" || sorbmodel == "oxygen") {
                double s=3.51, e=108.5303, m=15.9997;
                addAtomToProto(system,i,"O","O2", "M", -0.605,0,0,m,-0.113,0,e,s);
                addAtomToProto(system,i,"CoM","O2","M", 0,0,0,0,0.226,0,0,0);
                addAtomToProto(system,i,"O","O2", "M", 0.605, 0,0,m,-0.113,0,e,s);
                system.proto[i].name = "O2";
                system.proto[i].dof = 5;
            }
    
            // AMMONIA : TraPPE (hybrid: p-table mass, using CCCBDB experimental geometry for NH3: 1966 Herzberg)
            else if (sorbmodel == "nh3" || sorbmodel == "ammonia") {
                addAtomToProto(system,i,"M","NH3","M", 0,0,0.08, 0, -1.23, 0,0,0);
                addAtomToProto(system,i,"N","NH3","M", 0,0,0,14.0067, 0,0,185.0,3.42);
                addAtomToProto(system,i,"H","NH3","M", 0.0000,  -0.9377, -0.3816, 1.00794, 0.41, 0, 0,0);
                addAtomToProto(system,i,"H","NH3","M", 0.8121,  0.4689,  -0.3816, 1.00794, 0.41, 0,0,0);
                addAtomToProto(system,i,"H","NH3","M", -0.8121, 0.4689,  -0.3816, 1.00794, 0.41, 0,0,0);
                system.proto[i].name = "AMM";
                system.proto[i].dof = 6;
            }

            // OCTANE
            else if (sorbmodel == "octane" || sorbmodel == "c8h18") {
                // only 8 sites with params. H are dummies for viz.
                // MP2 6-31G* geom from cccbdb
                double q=0., a=0., se=3.75, sm=3.95, ee=98., em=46., ch3=15.0345, ch2=14.02658;
                addAtomToProto(system,i,"C1","OCT","M", -2.7969990,  3.5306670,   0.0000000,ch3,q,a,ee,se);
                addAtomToProto(system,i,"C2","OCT","M", 0.0002360,   0.7653030,   0.0000000,ch2 ,q, a,em,sm);
                addAtomToProto(system,i,"C3","OCT","M", -0.0002360,  -0.7653030 , 0.0000000,ch2,q,a,em,sm);
                addAtomToProto(system,i,"C4","OCT","M",-1.3965500,   1.3904300,   0.0000000   ,ch2,q,a,em,sm);
                addAtomToProto(system,i,"C5","OCT","M",1.3965500,    -1.3904300,  0.0000000   ,ch2,q,a,em,sm);
                addAtomToProto(system,i,"C6","OCT","M",-1.3965500,   2.9196300,   0.0000000  ,ch2,q,a,em,sm);
                addAtomToProto(system,i,"C7","OCT","M",1.3965500,    -2.9196300,  0.0000000  ,ch2,q,a,em,sm);
                addAtomToProto(system,i,"C8","OCT","M",2.7969990 ,   -3.5306670,  0.0000000,ch3,q,a,ee,se);
                addAtomToProto(system,i,"H","OCT","M",0.5528930,    1.1233910,  0.8793430        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",0.5528930,    1.1233910,  -0.8793430        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-0.5528930,   -1.1233910, 0.8793430        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-0.5528930,   -1.1233910, -0.8793430        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-1.9500370,   1.0332110,  -0.8794270        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-1.9500370,   1.0332110,  0.8794270        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",1.9500370,    -1.0332110, -0.8794270        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",1.9500370,    -1.0332110, 0.8794270        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-0.8416580,   3.2747210,  0.8781290        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-0.8416580,   3.2747210,  -0.8781290        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",0.8416580,    -3.2747210, 0.8781290        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",0.8416580,    -3.2747210, -0.8781290        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-2.7682840,   4.6251440,  0.0000000        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-3.3616990,   3.2150460,  -0.8828340        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",-3.3616990,   3.2150460,  0.8828340        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",2.7682840,    -4.6251440, 0.0000000        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",3.3616990,    -3.2150460, -0.8828340        ,0.,0.,0.,0.,0.);
                addAtomToProto(system,i,"H","OCT","M",3.3616990,    -3.2150460, 0.8828340        ,0.,0.,0.,0.,0.);
                system.proto[i].name = "OCT";
                system.proto[i].dof = 6;
            }
            else if (sorbmodel == "bf4" || sorbmodel=="bf4-") {
                // boron tetrafluoride anion (-1)
                // polarizability B alpha= 0.6634 from http://pubs.rsc.org/-/content/articlehtml/2017/cp/c7cp00924k#fn1
                double qB = 0.785727, qF = -0.446432;
                double sB = 4.083*system.constants.uff2mpmc, sF = 3.364*system.constants.uff2mpmc;
                double eB = 0.180/system.constants.kbk, eF = 0.05/system.constants.kbk;
                double aB = 0.6634, aF = 3.0013;
                double mB = 10.81, mF = 18.998;
                addAtomToProto(system,i,"B","BF4","M", 4.0557951459    ,   10.9438141063,       10.6815655407, mB, qB, aB, eB, sB);
                addAtomToProto(system,i,"F","BF4","M", 4.4811413817    ,   12.1935648431,       11.2443468872, mF, qF, aF, eF, sF);
 
                addAtomToProto(system,i,"F","BF4","M",  2.9022370967   ,    11.1553422360   ,     9.8652950248, mF, qF, aF, eF, sF);
                addAtomToProto(system,i,"F","BF4","M",  3.7305284247     ,  10.0328608252     ,  11.7431569709, mF, qF, aF, eF, sF);
                addAtomToProto(system,i,"F","BF4","M", 5.1075017722    ,   10.3969420588   ,     9.8785197493, mF, qF, aF, eF, sF);
                system.proto[i].name = "BF4";
                system.proto[i].dof = 6; 
            }
            // USER SORBATE MODEL NOT FOUND; ERROR OUT
            else {
                std::cout << "ERROR: The sorbate model name you supplied, " << sorbmodel.c_str() << ", was not found in the database. Check your spelling or use a manual model in your input atoms file."; printf("\n");
                std::exit(0);
            }

            // add proto to system molecule list if there's nothing there yet (i.e. no input file is OK)
            if (system.molecules.size() == 0) {
                system.molecules.push_back(system.proto[i]);
                system.stats.count_movables++;
                system.constants.total_atoms += system.proto[i].atoms.size();
            }

        } // end for all sorbate names
        } // end if sorbate name vector<string> has values (user inputs)

        // finally, zero the prototype coordinates and set fugacities (from user input)

        if (system.proto.size() > 0 && system.constants.sorbate_fugacity.size() > 0) { // this is needed to avoid seg fault
            // error out if the sorbates do not all have fugacities
            if (system.proto.size() != system.constants.sorbate_fugacity.size()) {
                printf("ERROR: The sorbate molecule count does not match the assigned fugacities count. Check your input file to see if all sorbates have a fugacity assigned.\n");
                std::exit(0);
            }
            for (int i=0; i<system.proto.size(); i++) {
            system.proto[i].fugacity = system.constants.sorbate_fugacity[i];
            system.proto[i].calc_center_of_mass();
            for (int j=0; j<system.proto[i].atoms.size(); j++) {
                for (int k=0; k<3; k++) system.proto[i].atoms[j].pos[k] -= system.proto[i].com[k];
            }
            system.proto[i].calc_center_of_mass();

            } // end protos loop
        } // end if we have protos yet

        // check for Manually-entered DOFs
        if (system.proto.size() > 0 && system.constants.sorbate_dof.size() > 0) {
            // use manual sorbate Degrees of Freedom
            if (system.proto.size() != system.constants.sorbate_dof.size()) {
                printf("ERROR: The number of manual sorbate DOFs you supplied (%i) does not match the number of sorbates in the system (%i)\n", (int)system.constants.sorbate_dof.size(), (int)system.proto.size()); 
                std::exit(0);
            }

            // apply the manual DOFs
            for (int z=0; z<system.proto.size(); z++)
                system.proto[z].dof = system.constants.sorbate_dof[z];
        }

        // finally, show the current proto molecules
        printf("\n::: PROTOTYPE (SORBATE) MOLECULES :::\n");
        for (int i=0; i<system.proto.size(); i++) {
        printf("\n:: %i :: Prototype molecule %i has PDBID %i ( name %s ) and has %i atoms\n",i,i, system.proto[i].PDBID, system.proto[i].name.c_str(), (int)system.proto[i].atoms.size());
            //system.proto.printAll();
            for (int j=0; j<system.proto[i].atoms.size(); j++) {
                system.proto[i].atoms[j].printAll();
            }
        }
        printf("\n");

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
    //system.pbc.printBasis();
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
    system.last.potential_sq = system.stats.potential_sq.value;
        // VOLUME
    system.last.volume = system.stats.volume.value;
        // CHEMICAL POTENTIAL
    system.last.chempot = system.stats.chempot.value;
        // Z
    system.last.z = system.stats.z.value;
        // QST
    system.last.qst = system.stats.qst.value;
    system.last.qst_nvt = system.stats.qst_nvt.value;

    for (int i=0; i<system.proto.size(); i++) {
        // DENSITY // Nmov
        system.last.density[i] = system.stats.density[i].value;
        system.last.Nmov[i] = system.stats.Nmov[i].value;
    }

        // N
    system.last.total_atoms = system.constants.total_atoms;
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
    system.stats.potential_sq.value = system.last.potential_sq;
        // VOLUME
    system.stats.volume.value = system.last.volume;
            // CHEMICAL POTENTIAL
    system.stats.chempot.value = system.last.chempot;
        // Z
    system.stats.z.value = system.last.z;
        // QST
    system.stats.qst.value = system.last.qst;
    system.stats.qst_nvt.value = system.last.qst_nvt;

    for (int i=0; i<system.proto.size(); i++) {
        // DENSITY // Nmov
        system.stats.density[i].value = system.last.density[i];
        system.stats.Nmov[i].value = system.last.Nmov[i];
    }

        // N
    system.constants.total_atoms = system.last.total_atoms;

}

void initialize(System &system) {

        // i only needed these for testing, I think they're useless right now,
        // so I'm not dealing with the problem of adding a name to each vector
        // instance of Nmov

        system.stats.Nsq.name = "Nsq";
        system.stats.NU.name = "NU";
        system.stats.qst.name = "qst";
        system.stats.rd.name = "rd";
        system.stats.es.name = "es";
        system.stats.polar.name = "polar";
        system.stats.potential.name = "potential";
        system.stats.volume.name = "volume";
        system.stats.z.name = "z";
        for (int i=0; i<system.proto.size(); i++) {
            system.stats.Nmov[i].name="Nmov";
            system.stats.wtp[i].name="wtp";
            system.stats.wtpME[i].name="wtpME";
            system.stats.movablemass[i].name="movablemass";
            system.stats.density[i].name="density";
        }
        system.stats.lj_lrc.name = "lj_lrc";
        system.stats.lj_self_lrc.name = "lj_self_lrc";
        system.stats.lj.name = "lj";
        system.stats.es_self.name = "es_self";
        system.stats.es_real.name = "es_real";
        system.stats.es_recip.name = "es_recip";
        system.stats.chempot.name = "chempot";
        system.stats.totalmass.name = "totalmass";
        system.stats.frozenmass.name = "frozenmass";
        system.stats.pressure.name = "pressure";
        system.stats.temperature.name= "temperature";
        system.stats.fdotr_sum.name="fdotr_sum";
        
        system.stats.dist_within.name="dist_within";
        system.stats.heat_capacity.name = "heat_capacity";
}

void setupFugacity(System &system) {
    // called in main()
    if (system.constants.fugacity_single == 1) {
        // we don't do anything unless user put a fugacity sorbate in input.
        if (system.constants.fugacity_single_sorbate == "h2") system.proto[0].fugacity = h2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "n2") system.proto[0].fugacity = n2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "co2") system.proto[0].fugacity = co2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "ch4") system.proto[0].fugacity = ch4_fugacity(system, system.constants.temp, system.constants.pres); 
        else if (system.constants.fugacity_single_sorbate == "off") system.proto[0].fugacity = system.constants.pres;     
    
    } else if (system.proto.size() == 1) {
        // presumably this means we want to use use I.G. approx for the sorbate
        system.proto[0].fugacity = system.constants.pres;
    }
    
    // (over)write CO2 fugacity with Luci's fits if user desires
    if (system.constants.co2_fit_fugacity) {
        printf("Using Laratelli CO2 fugacity fits with pressure = %f atm.\n",system.constants.pres);
        for (int i=0; i<system.constants.sorbate_name.size();i++) {
            if (system.constants.sorbate_name[i] == "co2_phast*" || system.constants.sorbate_name[i]=="co2_phastp") {
                system.proto[i].fugacity = phast_star(system.constants.pres);
            } else if (system.constants.sorbate_name[i] == "co2_phastq*" || system.constants.sorbate_name[i]=="co2_phastqp") {
                system.proto[i].fugacity = phast_q_star(system.constants.pres);
            } else if (system.constants.sorbate_name[i] == "co2_trappe") {
                system.proto[i].fugacity = trappe(system.constants.pres);
            }
        }
    }

    for (int i=0; i<system.proto.size(); i++) {
        printf("Fugacity for prototype %i ( %s ) = %f atm\n", i, system.proto[i].name.c_str(), system.proto[i].fugacity);
    }

    return;


}

void setupNBias(System &system) {
    // lots of messy math here but all I'm doing is converting
    // various units to/from N_molecules. 
    
    // frozen atoms total mass
    for (int c=0; c<system.molecules.size();c++) {
        for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].mass*system.constants.amu2kg/system.constants.NA;
            if (system.molecules[c].frozen) system.stats.frozenmass.value += thismass;
        }
    }
    

    string unit = system.constants.bias_uptake_unit;
    double x = system.constants.bias_uptake;
    double thevalue=0; // converted to N for safekeeping
    if (unit == "n" || unit == "") {
        thevalue = x;
    } else if (unit == "wt%") { 
        thevalue = (-x*system.stats.frozenmass.value/1000.) / system.proto[0].mass*system.constants.amu2kg / (x/100. -1) / 100.;
    } else if (unit == "wt%ME") {
        thevalue = x * (system.stats.frozenmass.value*1000)*system.constants.NA/1000/(system.proto[0].mass*system.constants.amu2kg*1000*system.constants.NA)/1000/10*100;
    } else if (unit == "cm^3/g") {
        thevalue = x * system.stats.frozenmass.value * system.constants.NA / 1000. / 22.4;
    } else if (unit == "mmol/g") {
        thevalue = x*system.stats.frozenmass.value*system.constants.NA/1000;
    } else if (unit == "mg/g") {
        thevalue = x * (system.stats.frozenmass.value*1000)*system.constants.NA/1000/(system.proto[0].mass*system.constants.amu2kg*1000*system.constants.NA)/1000;
    } else if (unit == "g/mL" || unit == "g/cm^3") {
        thevalue = x * 1e6 / 1e30 * (system.pbc.volume) * system.constants.NA / (system.proto[0].mass*system.constants.amu2kg*1000*system.constants.NA);
    }

    system.constants.bias_uptake = thevalue;
    printf("Calculated bias-uptake goal-N = %f molecules (in the box) = %f %s\n", thevalue, x, unit.c_str());
}


void fragmentMaker(System &system) {
    int numfrags = system.constants.numfrags; // number of base- frags to make
        // (will be multiplied by the number of atoms_per_frag's chosen
    
    int currentfrag = 0; // counter for fragments (local)
    int globalfrag = 0; // counter for fragments (global)
    int successfrags = 0; // counter for successful, non-duplicate, written fragments
    int currentatom = 0; // counter for atoms in a frag.
    double bondlength = system.constants.frag_bondlength; // Angstroms
    vector<vector<int>> fragment_atom_ids; // to check for duplicate fragments. 
    int duplicatefrags = 0; // counter for duplicates just to inform user

    // first, get the unique atom names.
    // because we want fragments which represent 
    // each atom in a buried environment
    printf("\n:: Making unique atoms vector for fragment creation.\n");
    vector<string> atomlabels;    
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) {
            for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                string theatom = system.molecules[i].atoms[j].name;
                if (std::find(atomlabels.begin(), atomlabels.end(), theatom) != atomlabels.end()) {
                    // the atom is already in our unique set "atomlabels"
                    continue;
                } else {
                    atomlabels.push_back(theatom);
                }
            }
        }
    }

    // confirm the unique atoms to user
    for (int x=0; x<atomlabels.size(); x++)
        printf(":: ---> Unique atom %i = %s\n", x+1, atomlabels[x].c_str());

    // loop through fragsizes (could be multiple depending on user input)
    for (int fs=0; fs<(int)system.constants.fragsize.size(); fs++) {
        // the size for this round of frags
        int fragsize = (int)system.constants.fragsize[fs];
        currentfrag=0; // reset
        printf("Making %i-atom fragments...\n", fragsize);
    // make fragment for each unique atom
    // cycles through uniques until numfrags is reached.
    for (int x=0; x<atomlabels.size(); x++) {
        string theatom = atomlabels[x];
        currentatom=0; // reset
        bondlength = system.constants.frag_bondlength; // reset

        vector<Atom> currentFrag; // the fragment being built
        vector<vector<double>> building_points; // holds the coordinates of the fractal'd focus-atoms
        vector<int> fragPDBIDs; // holds the PDBIDs of each atom in the frag to avoid duplicates  
            // get the most central occurrence of the atom
            int found=0;
            double probe_radius=0.0;
            while (!found) {
                probe_radius += 3.0; // expand by 3.0 A each time to find the atom
            
                for (int i=0;i<system.molecules.size();i++) {
                    if (!system.molecules[i].frozen) continue; // skip movers, we want frozens
                    for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                        if (system.molecules[i].atoms[j].name != theatom) continue; // only consider atoms of interest
                        // found-it checker
                        if (fabs(system.molecules[i].atoms[j].pos[0]) < probe_radius &&
                            fabs(system.molecules[i].atoms[j].pos[1]) < probe_radius &&
                            fabs(system.molecules[i].atoms[j].pos[2]) < probe_radius) {
                            found=1;
                            probe_radius=0;
                            vector<double> tempcoords = vector<double>(3);
                            for (int n=0;n<3;n++) tempcoords[n] = system.molecules[i].atoms[j].pos[n];
                            building_points.push_back(tempcoords);
                            currentFrag.push_back(system.molecules[i].atoms[j]); // add the atom to frag.
                            fragPDBIDs.push_back(system.molecules[i].atoms[j].PDBID); // and the PDBID to the vector for housekeeping
                            currentatom++;
                        } 
                    } // end j
                } // end i
            } // end while loop to find central atom
 
            // build it, son.
            double r; // pair distance
            int builders_size; // = (int)building_points.size();
            while (currentatom < fragsize) {   
                // find atoms (pseudo-)bonded to the builder atom, fractal style
                builders_size = (int)building_points.size();
                vector<vector<double>> temp_building_points;
                for (int bl=0; bl<builders_size; bl++) {
                for (int i=0;i<system.molecules.size();i++) {
                    if (!system.molecules[i].frozen) continue;
                    for (int j=0;j<system.molecules[i].atoms.size();j++) {
                        if (!system.molecules[i].atoms[j].frozen) continue;
                        double temp[3];
                        for (int n=0;n<3;n++) temp[n] = building_points[bl][n];
                        double * distances = getR(system, temp, system.molecules[i].atoms[j].pos, 0);
                        r = distances[3];
                        if (r <= bondlength) { // look for a bonded atom
                            // make sure it's not a duplicate atom
                            if (std::find(fragPDBIDs.begin(), fragPDBIDs.end(), system.molecules[i].atoms[j].PDBID) != fragPDBIDs.end()) {
                                continue; 
                            }
                            vector<double> tempvec = vector<double>(3);
                            for (int n=0;n<3;n++) tempvec[n] = system.molecules[i].atoms[j].pos[n];
                            temp_building_points.push_back(tempvec);
                            currentFrag.push_back(system.molecules[i].atoms[j]);
                            fragPDBIDs.push_back(system.molecules[i].atoms[j].PDBID);
                            currentatom++;
                        } // end if bonded-atom
                        if (currentatom >= fragsize) break; // make sure we kill the loop as soon as the max # of atoms hits (don't wait until this inner loop finishes)
                    } // end j
                } // end i (the pair-loop is over now.)
                } // end loop through builder atoms
                    //printf("made it here\n");
                    // if no additional bonders were detected, boost the bond-length so we can find some atoms
                    if (building_points == temp_building_points || temp_building_points.size()==0) {
                        bondlength += 0.1;
                    } else {
                        // reset the current (untapped) building points if new connections were found
                        building_points.clear();
                        building_points = temp_building_points;    
                    }
            } // end while loop adding atoms to frag



            // frag has been made.
            // write the frag and move to next
            // first check duplicate fragment
            int dupeFlag=0;
            sort( fragPDBIDs.begin(), fragPDBIDs.end() ); // sort so that the dupes are correctly found
            for (int fi=0; fi<fragment_atom_ids.size(); fi++) {
                if (fragPDBIDs == fragment_atom_ids[fi]) dupeFlag = 1;
            }

            // write the new fragment now, since it's not a duplicate.
            if (dupeFlag == 0) {
            string fragid = to_string(globalfrag);
            std::string fragtail = std::string(5 - fragid.length(), '0') + fragid;
            char suffix[7] = "-";
            strcat(suffix, fragtail.c_str());
            char filename[20] = "fragment";
            strcat(filename, suffix);
            char extension[5] = ".xyz";
            strcat(filename, extension);

            FILE *f = fopen(filename, "w");
            if (f == NULL) {
                printf("Error opening fragment file! Name = %s\n", filename);
                exit(EXIT_FAILURE);
            }
            fprintf(f, "%i\n", (int)currentFrag.size());
            double charge_sum=0.;
            for (int i=0; i<currentFrag.size(); i++) charge_sum += currentFrag[i].C / system.constants.E2REDUCED;
    
            fprintf(f, "Fragment # %i centered around %s with total charge = %.9f e produced by MCMD\n", globalfrag+1, theatom.c_str(), charge_sum);
            for (int i=0; i<currentFrag.size(); i++) {
                fprintf(f, "%s %.9f %.9f %.9f %.7f\n", currentFrag[i].name.c_str(), currentFrag[i].pos[0], currentFrag[i].pos[1], currentFrag[i].pos[2], currentFrag[i].C / system.constants.E2REDUCED);
            }
            fclose(f);

            successfrags++;
            printf("Built fragment %i with filename %s centered on %s\n", globalfrag+1, filename, theatom.c_str());
            currentfrag++;
            globalfrag++;
            fragment_atom_ids.push_back(fragPDBIDs); // to check duplicates
            } else {
                // the frag was a duplicate, so we skip it from writing
                duplicatefrags++;
                currentfrag++;
                globalfrag++;
            }

            if (x == ((int)atomlabels.size()-1) && currentfrag < numfrags) {
                x=-1; continue; // keep cycling the unique atoms if numfrags not reached yet
                // negative one bc it's about to do x++ again
            }
            if (currentfrag >= numfrags)  {
                break; // break out if numfrags reached
            }
            
        
    } // end x loop (unique atoms loop, which makes 1 frag per iter)
    } // end fs loop (frag-size loop, for different size frags)

    printf("%i duplicate fragments were detected and not written.\n", duplicatefrags);
    printf("%i fragments were successfully created.\n", successfrags);

}

void setupCrystalBuild(System &system) {
    // build out a periodic box according to user's desire
    int xdim,ydim,zdim; double xlen,ylen,zlen;
    int i,j;

        xdim = system.constants.crystalbuild_x;
        ydim = system.constants.crystalbuild_y;
        zdim = system.constants.crystalbuild_z;
        
        system.constants.free_volume *= (xdim*ydim*zdim); // adjusts free volume by multiplicity

        xlen= system.pbc.x_length;
        ylen = system.pbc.y_length;
        zlen = system.pbc.z_length;

        double origa = system.pbc.a;
        double origb = system.pbc.b;
        double origc = system.pbc.c;

        double orig_basis[3][3];
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                orig_basis[i][j] = system.pbc.basis[i][j];


        // assumes ONE frozen molecule!!
        int size = (int)system.molecules.size();
        int asize=0;
        for (int i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) {
                asize = system.molecules[i].atoms.size();
                break;
            }
        }

        printf("Building out crystal by %ix, %iy, %iz of the original.\n", xdim,ydim,zdim);
        printf(" --> using xlen = %f; ylen = %f; zlen = %f;\n", xlen,ylen,zlen);

        // DUPLICATE IN X
        if (fabs(xdim) > 1) {
            for (int iter=0; iter < fabs(xdim)-1; iter++) {
            system.pbc.a += origa;
            system.pbc.calcNormalBasis(); 
            setupBox(system);
            
            for (i =0; i <size; i++) {
                if (system.molecules[i].frozen ) {
                    for (j=0; j<asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        
                        for (int n=0;n<3;n++) {
                            if (xdim>0)
                                newatom.pos[n] += orig_basis[0][n]*(iter+1);
                            else
                                newatom.pos[n] -= orig_basis[0][n]*(iter+1);
                        }
                        
                        system.molecules[i].mass += newatom.mass;
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++; 
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        for (int n=0;n<3;n++) {
                            if (xdim>0) 
                                newmolecule.atoms[j].pos[n] += orig_basis[0][n]*(iter+1);
                            else
                                newmolecule.atoms[j].pos[n] -= orig_basis[0][n]*(iter+1);
                        }
                    }
                    system.stats.count_movables++;
                    system.molecules.push_back(newmolecule);
                }
            }
            } // end iterations x
        // reset the atom count
        size=(int)system.molecules.size(); 
        asize=0;
        for (int i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) {
                asize = system.molecules[i].atoms.size();
                break;
            }
        }


        }
        // DUPLICATE IN Y
        if (fabs(ydim) > 1) {
            for (int iter=0; iter < fabs(ydim)-1; iter++) {
            system.pbc.b += origb;
            system.pbc.calcNormalBasis();
            setupBox(system);

            for (i=0;i<size; i++) {
                if (system.molecules[i].frozen) {
                    for (j=0; j< asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        for (int n=0;n<3;n++) {
                            if (ydim>0)
                                newatom.pos[n] += orig_basis[1][n]*(iter+1);                    
                            else 
                                newatom.pos[n] -= orig_basis[1][n]*(iter+1);
                        }
                        
                        system.molecules[i].mass += newatom.mass;                     
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++;
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        for (int n=0;n<3;n++) {
                            if (ydim>0)
                                newmolecule.atoms[j].pos[n] += orig_basis[1][n]*(iter+1);
                            else
                                newmolecule.atoms[j].pos[n] -= orig_basis[1][n]*(iter+1);
                        }
                    }
                    system.stats.count_movables++;
                    system.molecules.push_back(newmolecule);
                }

            }
            } // end iterations y
        // reset the atom count.
        size=(int)system.molecules.size();
        asize=0;
        for (int i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) {
                asize = system.molecules[i].atoms.size();
                break;
            }
        }
        }
        // DUPLICATE IN Z
        if (fabs(zdim) > 1) {
        for (int iter=0; iter < fabs(zdim)-1; iter++) {
            system.pbc.c += origc;
            system.pbc.calcNormalBasis();
            setupBox(system);

            for (i=0;i<size; i++) {
                if (system.molecules[i].frozen) {
                    for (j=0; j<asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        for (int n=0;n<3;n++) {
                            if (zdim>0)
                                newatom.pos[n] += orig_basis[2][n]*(iter+1);                
                            else
                                newatom.pos[n] -= orig_basis[2][n]*(iter+1);        

                        }
                        system.molecules[i].mass += newatom.mass;
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++;
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        for (int n=0;n<3;n++) {
                            if (zdim>0)
                                newmolecule.atoms[j].pos[n] += orig_basis[2][n]*(iter+1);
                            else
                                newmolecule.atoms[j].pos[n] -= orig_basis[2][n]*(iter+1);
                        }
                    }
                    system.stats.count_movables++;
                    system.molecules.push_back(newmolecule);
                }
            }
            } // end iterations z
        // reset the atom count
        size=(int)system.molecules.size();
        asize=0;
        for (int i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) {
                asize = system.molecules[i].atoms.size();
                break;
            }
        }
        } 
    

    printf("Done building crystal.\n\n");
    if (system.constants.autocenter) centerCoordinates(system);
    setupBox(system); // final box setup JIC

    // this is a just in case measure. If user inputs an xyz file we need to re-write the PDBIDs before proceeding.
    // CONSOLIDATE ATOM AND MOLECULE PDBID's
            // quick loop through all atoms to make PDBID's pretty (1->N)
            if (system.molecules.size() > 0) {
            int molec_counter=1, atom_counter=1;
            for (int i=0; i<system.molecules.size(); i++) {
                system.molecules[i].PDBID = molec_counter;
                for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                    system.molecules[i].atoms[j].PDBID = atom_counter;
                    system.molecules[i].atoms[j].mol_PDBID = system.molecules[i].PDBID;
                    atom_counter++;
                } // end loop j
                molec_counter++;
            } // end loop i
            } // end if molecules exist
}

void scaleCharges(System &system) {
    for (int i=0; i < system.molecules.size(); i++) {
        for (int j=0; j < system.molecules[i].atoms.size(); j++) {
            if (system.molecules[i].atoms[j].frozen)
                system.molecules[i].atoms[j].C *= system.constants.scale_charges_factor;    
        }
    }
    printf("Finished scaling all atomic charges by %f.\n", system.constants.scale_charges_factor);
}

double getDOF(System &system) {
    // needs attention for more general definition
    int N = system.stats.count_movables;
    return (N>0) ? 3.0*N-3.0 : 3.0; 
}

int getNlocal(System &system, int protoid) {
    // get the number of molecules of a given sorbate-type
    int N=0;
    for (int i=0;i<system.molecules.size(); i++)
        if (system.molecules[i].name == system.proto[protoid].name)
            N++;

    return N;
}

void initialVelMD(System &system, int startupflag) {
        double randv;
        int userflag=0, i, n, z;
        double v_init, v_init_AVG;
        
        // user-defined velocities
        if (startupflag && system.constants.md_manual_init_vel) {
            userflag=1;
            for (i=0; i<system.molecules.size(); i++) {
                for (n=0; n<3; n++) {
                    randv = (getrand()*2 - 1) * sqrt(system.constants.md_init_vel * system.constants.md_init_vel/3.);
                    system.molecules[i].vel[n] = randv; // that is, +- 0->1 * user param
                }
            }
            printf("Assigned initial velocities via user-defined value: %f A/fs\n",system.constants.md_init_vel);
             
        // thermostat information (and apply init. vel if needed)
        }
        if (system.constants.thermostat_type == THERMOSTAT_ANDERSEN || system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) {
            if (system.constants.md_mode == MD_MOLECULAR) {
            v_init_AVG=0;
            for (z=0; z<system.proto.size(); z++) {   

                v_init = 1e-5*sqrt(system.constants.kb*system.constants.temp*system.proto[z].dof / (system.proto[z].mass*system.constants.amu2kg));
                system.proto[z].md_velx_goal = sqrt(v_init*v_init/3.);

                // apply thermostat velocities if user did not specify an init. vel.
                if (!userflag) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.proto[z].name == system.molecules[i].name) {
                        v_init_AVG += v_init * (1. / (double)system.stats.count_movables);
                        for (n=0;n<3;n++) {
                            double pm = (getrand() > 0.5) ? 1.0 : -1.0; 
                            system.molecules[i].vel[n] = pm * system.proto[z].md_velx_goal;
                        }
                    }
                } // end molecules loop i
                system.constants.md_init_vel = v_init_AVG;
                } // end if user-defined init. vel's
            } // end prototype loop z
            }
            // FLEXIBLE MODELING -- "FREE" BONDED ATOMS
            else if (system.constants.md_mode == MD_FLEXIBLE || system.constants.md_mode == MD_ATOMIC) {
                v_init_AVG=0;
                int N = (system.constants.flexible_frozen ? system.constants.total_atoms : system.constants.total_atoms - system.stats.count_frozens);
                unsigned int dof_total=3.*N - 3.0; // - (int)system.constants.uniqueBonds.size();
                for (unsigned int i=0;i<system.molecules.size();i++) {
                    if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
                for (unsigned int j=0;j<system.molecules[i].atoms.size();j++) {
                    // normalized atom velocity based on total DOF from bonds and such
                    v_init = 1e-5*sqrt(system.constants.kb*system.constants.temp*dof_total/N / (system.molecules[i].atoms[j].mass*system.constants.amu2kg));
                    v_init_AVG += v_init;
                    system.molecules[i].atoms[j].md_velx_goal = sqrt(v_init*v_init/3.);
                    // apply velocities to atoms
                    if (!userflag) {
                        for (n=0;n<3;n++) {
                            double pm = (getrand() > 0.5) ? 1.0 : -1.0;
                            system.molecules[i].atoms[j].vel[n] = pm * system.molecules[i].atoms[j].md_velx_goal;
                        }
                    }
                }
            }
                v_init_AVG /= (double)N;
                system.constants.md_init_vel = v_init_AVG;
            } // end if flexible MD
        } // end if thermostat

    // write initial velocities to the saved array for each molecule
    for (unsigned int i=0; i<system.molecules.size();i++) {
        for (unsigned int n=0;n<3;n++) {
            system.molecules[i].original_vel[n] = system.molecules[i].vel[n];
            if (system.molecules[i].vel[n] != 0) system.constants.zero_init_vel_flag = 0; // flag this so that we can change original_vel to make VACF nonzero on first step.
        }
    }
    

} // end initialVelMD()

void consolidatePDBIDs(System &system) {
// CONSOLIDATE ATOM AND MOLECULE PDBID's
            // quick loop through all atoms to make PDBID's pretty (1->N)
            int molec_counter=1, atom_counter=1;
            for (int i=0; i<system.molecules.size(); i++) {
                system.molecules[i].PDBID = molec_counter;
                for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                    system.molecules[i].atoms[j].PDBID = atom_counter;
                    system.molecules[i].atoms[j].mol_PDBID = system.molecules[i].PDBID;
                    atom_counter++;
                } // end loop j
                molec_counter++;
            } // end loop i            
}

int getNElectrons(System &system, int molid) {
    int i;
    int numelec=0;
    for (i=0; i < system.molecules[molid].atoms.size(); i++) {
         numelec += system.constants.elements[system.molecules[molid].atoms[i].name];
    }
    return numelec - system.constants.user_charge; // so if user says it's +2, we remove those.
}

int calcAtomPairs(System &system, int unique) {
    int pairs = 0;
    if (!unique) {
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            for (int k=0;k<system.molecules.size();k++) {
                for (int l=0;l<system.molecules[k].atoms.size();l++) {
                    pairs++;
                }
            }
        }
    }
    } else { // unique pairs ((N^2 - N) / 2)
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            for (int k=i+1;k<system.molecules.size();k++) {
                for (int l=0;l<system.molecules[k].atoms.size();l++) {
                    pairs++;
                }
            }
        }
    }

    }
    return pairs;
}

void updateMolecularDOFs(System &system) {
    unsigned int i,z;
    // updates the DOF of each molecule in the system
    // based on its prototype name info
    for (z=0;z<system.proto.size();z++) {
        for (i=0;i<system.molecules.size();i++) {
            if (system.proto[z].name == system.molecules[i].name) {
                system.molecules[i].dof = system.proto[z].dof;
            }
        }
    }
}
