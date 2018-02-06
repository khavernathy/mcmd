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
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
using namespace std;

double getrand() {
    return drand48();
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
    atom.m = mass * system.constants.cM;
    atom.C = charge * system.constants.E2REDUCED;
    atom.polar = polarizability;
    atom.eps = epsilon;
    atom.sig = sigma;

    system.proto[protoid].mass += atom.m;
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
            // CARBON DIOXIDE CO2
            else if (sorbmodel == "co2_phast") {
                addAtomToProto(system, i,"COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77106, 0.0, 8.52238, 3.05549);
                addAtomToProto(system, i,"COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"CON", "CO2", "M", 1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                system.proto[i].name = "CO2";
                system.proto[i].dof = 5;
            }
            else if (sorbmodel == "co2_phast*") {
                addAtomToProto(system, i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77134, 1.22810, 19.61757, 3.03366);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "CON", "CO2", "M", 1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
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
                double s=3.02, e=49.0, m=7.9997;
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

        system.stats.fdotrsum.name="fdotrsum";
        system.stats.dist_within.name="dist_within";
        system.stats.heat_capacity.name = "heat_capacity";
        system.stats.diffusion.name="diffusion";
}

void setupFugacity(System &system) {
    // called in main()
    if (system.constants.fugacity_single == 1) {
        // we don't do anything unless user put a fugacity sorbate in input.
        if (system.constants.fugacity_single_sorbate == "h2") system.proto[0].fugacity = h2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "n2") system.proto[0].fugacity = n2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "co2") system.proto[0].fugacity = co2_fugacity(system.constants.temp, system.constants.pres);
        else if (system.constants.fugacity_single_sorbate == "ch4") system.proto[0].fugacity = ch4_fugacity(system.constants.temp, system.constants.pres); 
        else if (system.constants.fugacity_single_sorbate == "off") system.proto[0].fugacity = system.constants.pres;     
    
    } else if (system.proto.size() == 1) {
        // presumably this means we want to use use I.G. approx for the sorbate
        system.proto[0].fugacity = system.constants.pres;
    }
    
    for (int i=0; i<system.proto.size(); i++) {
        printf("Calculated fugacity for prototype %i ( %s ) = %f atm\n", i, system.proto[i].name.c_str(), system.proto[i].fugacity);
    }

    return;


}

void setupNBias(System &system) {
    // lots of messy math here but all I'm doing is converting
    // various units to/from N_molecules. 
    
    // frozen atoms total mass
    for (int c=0; c<system.molecules.size();c++) {
        for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
            if (system.molecules[c].frozen) system.stats.frozenmass.value += thismass;
        }
    }
    

    string unit = system.constants.bias_uptake_unit;
    double x = system.constants.bias_uptake;
    double thevalue=0; // converted to N for safekeeping
    if (unit == "N" || unit == "") {
        thevalue = x;
    } else if (unit == "wt%") { 
        thevalue = (-x*system.stats.frozenmass.value/1000.) / system.proto[0].mass / (x/100. -1) / 100.;
    } else if (unit == "wt%ME") {
        thevalue = x * (system.stats.frozenmass.value*1000)*system.constants.NA/1000/(system.proto[0].mass*1000*system.constants.NA)/1000/10*100;
    } else if (unit == "cm^3/g") {
        thevalue = x * system.stats.frozenmass.value * system.constants.NA / 1000. / 22.4;
    } else if (unit == "mmol/g") {
        thevalue = x*system.stats.frozenmass.value*system.constants.NA/1000;
    } else if (unit == "mg/g") {
        thevalue = x * (system.stats.frozenmass.value*1000)*system.constants.NA/1000/(system.proto[0].mass*1000*system.constants.NA)/1000;
    } else if (unit == "g/mL" || unit == "g/cm^3") {
        thevalue = x * 1e6 / 1e30 * (system.pbc.volume) * system.constants.NA / (system.proto[0].mass*1000*system.constants.NA);
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
        if (xdim > 1) {
            for (int iter=0; iter < xdim-1; iter++) {
            system.pbc.a += origa;
            system.pbc.calcNormalBasis(); 
            setupBox(system);
            for (i =0; i <size; i++) {
                if (system.molecules[i].frozen ) {
                    for (j=0; j<asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        newatom.pos[0] += xlen*(iter+1);
                        system.molecules[i].mass += newatom.m;
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++; 
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        newmolecule.atoms[j].pos[0] += xlen*(iter+1);
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
        if (ydim > 1) {
            for (int iter=0; iter < ydim-1; iter++) {
            system.pbc.b += origb;
            system.pbc.calcNormalBasis();
            setupBox(system);

            for (i=0;i<size; i++) {
                if (system.molecules[i].frozen) {
                    for (j=0; j< asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        newatom.pos[1] += ylen*(iter+1);
                        if (system.pbc.gamma != 90.0) newatom.pos[0] -= (iter+1)*origb*sin((system.pbc.gamma-90.0)*M_PI/180.);
                        system.molecules[i].mass += newatom.m;                     
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++;
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        newmolecule.atoms[j].pos[1] += ylen*(iter+1);
                        if (system.pbc.gamma != 90.0) newmolecule.atoms[j].pos[0] -= (iter+1)*origb*sin((system.pbc.gamma-90.0)*M_PI/180.);
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
        if (zdim > 1) {
            for (int iter=0; iter < zdim-1; iter++) {
            system.pbc.c += origc;
            system.pbc.calcNormalBasis();
            setupBox(system);

            for (i=0;i<size; i++) {
                if (system.molecules[i].frozen) {
                    for (j=0; j<asize; j++) {
                        Atom newatom = system.molecules[i].atoms[j];
                        newatom.pos[2] += zlen*(iter+1);
                        if (system.pbc.beta != 90.0) newatom.pos[0] -= (iter+1)*origc*sin((system.pbc.beta-90.0)*M_PI/180.0);
                        system.molecules[i].mass += newatom.m;
                        system.molecules[i].atoms.push_back(newatom);
                        system.constants.total_atoms++;
                        system.stats.count_frozens++;
                    }
                } else if (system.constants.crystalbuild_includemovers) {
                    Molecule newmolecule = system.molecules[i];
                    for (j=0; j<newmolecule.atoms.size(); j++) {
                        system.constants.total_atoms++;
                        newmolecule.atoms[j].pos[2] += zlen*(iter+1);
                        if (system.pbc.beta != 90.0) newmolecule.atoms[j].pos[0] -= (iter+1)*origc*sin((system.pbc.beta-90.0)*M_PI/180.0);
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

void initialVelMD(System &system) {
        double randv;
        int userflag=0, i, n, z;
        double v_init, v_init_AVG;
        
        // user-defined velocities
        if (system.constants.md_manual_init_vel) {
            userflag=1;
            for (i=0; i<system.molecules.size(); i++) {
                for (n=0; n<3; n++) {
                    randv = (getrand()*2 - 1) * sqrt(system.constants.md_init_vel * system.constants.md_init_vel/3.);
                    system.molecules[i].vel[n] = randv; // that is, +- 0->1 * user param
                }
            }
            printf("Computed initial velocities via user-defined value: %f A/fs\n",system.constants.md_init_vel);
             
        // thermostat information (and apply init. vel if needed)
        }
        if (system.constants.thermostat_type == THERMOSTAT_ANDERSEN || system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) {
            v_init_AVG=0;
            for (z=0; z<system.proto.size(); z++) {   

                v_init = 1e-5*sqrt(system.constants.kb*system.constants.temp*system.proto[z].dof / system.proto[z].mass);
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
            
        } // end if thermostat
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


