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
    //for (int j=0; j<system.molecules[i].atoms.size(); j++) {
        if (system.molecules[i].com[0] > system.pbc.x_max) { // right of box
            system.molecules[i].diffusion_corr[0] += system.pbc.x_length;
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[0] -= system.pbc.x_length;
                }
        }
        else if (system.molecules[i].com[0] < system.pbc.x_min) {
            system.molecules[i].diffusion_corr[0] -= system.pbc.x_length;
                for (int k=0; k<system.molecules[i].atoms.size(); k++) {
	                system.molecules[i].atoms[k].pos[0] += system.pbc.x_length;
                }
        }
        if (system.molecules[i].com[1] > system.pbc.y_max) {
            system.molecules[i].diffusion_corr[1] += system.pbc.y_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] -= system.pbc.y_length;
                }
        }
        else if (system.molecules[i].com[1] < system.pbc.y_min) {
            system.molecules[i].diffusion_corr[1] -= system.pbc.y_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[1] += system.pbc.y_length;
                }
        }
        if (system.molecules[i].com[2]  > system.pbc.z_max) {
            system.molecules[i].diffusion_corr[2] += system.pbc.z_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] -= system.pbc.z_length;
                }
        }
        else if (system.molecules[i].com[2] < system.pbc.z_min) {
            system.molecules[i].diffusion_corr[2] -= system.pbc.z_length;
	            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                    system.molecules[i].atoms[k].pos[2] += system.pbc.z_length;
                }
        }
    //} // end atom loop
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
            system.molecules[i].diffusion_corr[0] -= system.pbc.x_length;
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
        else if (system.molecules[i].com[0] > box_limit[1]) { // right of box
            system.molecules[i].diffusion_corr[0] += system.pbc.x_length;
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
            system.molecules[i].diffusion_corr[1] -= system.pbc.y_length;
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
        else if (system.molecules[i].com[1] > box_limit[3]) { // above box
            system.molecules[i].diffusion_corr[1] += system.pbc.y_length;
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
            system.molecules[i].diffusion_corr[2] -= system.pbc.z_length;
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
        else if (system.molecules[i].com[2] > box_limit[5]) { // in front of box
            system.molecules[i].diffusion_corr[2] += system.pbc.z_length;
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
    // send it over
    system.proto[protoid].atoms.push_back(atom);

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

            std::cout << "THE SORB MODEL WAS SUPPLIED: " << sorbmodel.c_str(); printf("\n");
            // each call takes 12 arguments
            // HYDROGEN H2
            if (sorbmodel == "h2_buch") {
                addAtomToProto(system, i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 2.016, 0.0, 0.0, 34.2, 2.96);
                system.proto[i].name = "H2";
            }
            else if (sorbmodel == "h2_bss") {
                addAtomToProto(system, i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.74640, 0.0, 8.85160, 3.2293);
                addAtomToProto(system, i, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, i, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"H2N", "H2", "M", 0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                addAtomToProto(system, i,"H2N", "H2", "M", -0.329, 0.0, 0.0, 0.0, 0.0, 0.0, 4.06590, 2.3406);
                system.proto[i].name = "H2";
            }
            else if (sorbmodel == "h2_dl") {
                addAtomToProto(system, i,"H2G", "H2", "M", 0.00, 0.0, 0.0, -0.93600, 0.0, 36.7, 2.958, 0.0);
                addAtomToProto(system, i,"H2E", "H2", "M", -0.370, 0.0, 0.0, 1.008, 0.468, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"H2E", "H2", "M", 0.370, 0.0, 0.0, 1.008, 0.468, 0.0, 0.0, 0.0);
                system.proto[i].name = "H2";
            }
            else if (sorbmodel == "h2_bssp") {
                addAtomToProto(system,i, "H2G", "H2", "M", 0.0, 0.0, 0.0, 0.0, -0.7464, 0.6938, 12.76532, 3.15528);
                addAtomToProto(system,i, "H2E", "H2", "M", 0.371, 0.0, 0.0, 1.008, 0.3732, 0.00044, 0.0, 0.0);
                addAtomToProto(system,i, "H2E", "H2", "M", -0.371, 0.0, 0.0, 1.008, 0.37320, 0.00044, 0.0, 0.0);
                addAtomToProto(system,i, "H2N", "H2", "M", 0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                addAtomToProto(system,i, "H2N", "H2", "M", -0.363, 0.0, 0.0, 0.0, 0.0, 0.0, 2.16726, 2.37031);
                system.proto[i].name = "H2";
            }
            // HELIUM He (useful for theor. calc. of pore vol.). See http://pubs.acs.org/doi/abs/10.1021/jp050948l
            // and http://onlinelibrary.wiley.com/doi/10.1002/aic.690470521/abstract
            // and our paper ESI http://pubs.rsc.org/en/Content/ArticleLanding/2014/CC/c4cc03070b#!divAbstract
            else if (sorbmodel == "he" || sorbmodel == "He") {
                addAtomToProto(system, i,"He", "He", "M", 0.0, 0.0, 0.0, 4.002602, 0.0, 0.0, 10.220, 2.280);
            }
            // CARBON DIOXIDE CO2
            else if (sorbmodel == "co2_phast") {
                addAtomToProto(system, i,"COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77106, 0.0, 8.52238, 3.05549);
                addAtomToProto(system, i,"COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38553, 0.0, 0.0, 0.0);
                addAtomToProto(system, i,"CON", "CO2", "M", 1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.091, 0.0, 0.0, 0.0, 0.0, 0.0, 76.76607, 2.94473);
                system.proto[i].name = "CO2";
            }
            else if (sorbmodel == "co2_phast*") {
                addAtomToProto(system, i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.0107, 0.77134, 1.22810, 19.61757, 3.03366);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.162, 0.0, 0.0, 15.9994, -0.38567, 0.73950, 0.0, 0.0);
                addAtomToProto(system, i, "CON", "CO2", "M", 1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                addAtomToProto(system, i,"CON", "CO2", "M", -1.208, 0.0, 0.0, 0.0, 0.0, 0.0, 46.47457, 2.99429);
                system.proto[i].name = "CO2" ;
            }
            else if (sorbmodel == "co2_trappe") {
                addAtomToProto(system,i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.7, 0.0, 27.0, 2.80);
                addAtomToProto(system,i, "COE", "CO2", "M", 1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05);
                addAtomToProto(system,i, "COE", "CO2", "M", -1.160, 0.0, 0.0, 16.0, -0.35, 0.0, 79.0, 3.05);
                system.proto[i].name = "CO2";
            }
            else if (sorbmodel == "co2_becker") {
                // A polar model adapted from TraPPE. Becker et al. 10.1021/acs.jpcc.6b12052
                addAtomToProto(system,i, "COG", "CO2", "M", 0.0, 0.0, 0.0, 12.01, 0.7, 0.916, 23.4, 2.8);
                addAtomToProto(system, i, "COE", "CO2", "M", 1.16, 0.0, 0.0, 16.0, -0.35, 0.575, 73.08, 3.05);
                addAtomToProto(system, i, "COE", "CO2", "M", -1.16, 0.0, 0.0, 16.0, -0.35, 0.575, 73.08, 3.05);
            }
            // NITROGEN N2
            else if (sorbmodel == "n2_mcquarrie") {
                addAtomToProto(system,i, "N2G", "N2", "M", 0.0, 0.0, 0.0, 28.01344, 0.0, 0.0, 95.1, 3.7);
                system.proto[i].name = "N2";
            }
            // METHANE CH4
            else if (sorbmodel == "ch4_trappe") {
                addAtomToProto(system,i, "CHG", "CH4", "M", 0.0, 0.0, 0.0, 16.0426, 0.0, 0.0, 148.0, 3.73);
                system.proto[i].name = "CH4";
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
            }
            // ACETYLENE C2H2
            else if (sorbmodel == "c2h2" || sorbmodel == "acetylene") {
                addAtomToProto(system,i, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system,i, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 0.0, 81.35021, 3.40149);
                addAtomToProto(system,i, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                addAtomToProto(system,i, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.0, 0.00026, 4.77683);
                system.proto[i].name = "ACE";
            }
            else if (sorbmodel == "c2h2*" || sorbmodel == "acetylene*") {
                addAtomToProto(system,i, "C2G", "ACE", "M", 0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system,i, "C2G", "ACE", "M", -0.605, 0.0, 0.0, 12.011, -0.29121, 1.55140, 70.81797, 3.42964);
                addAtomToProto(system,i, "H2G", "ACE", "M", 1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                addAtomToProto(system,i, "H2G", "ACE", "M", -1.665, 0.0, 0.0, 1.008, 0.29121, 0.14480, 0.00026, 4.91793);
                system.proto[i].name = "ACE";
            }
            // ETHLYENE C2H4
            else if (sorbmodel == "c2h4" || sorbmodel == "ethlyene" || sorbmodel == "ethene") {
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 0.0, 69.08116, 3.51622);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.0, 3.169, 2.41504);
                system.proto[i].name = "ETH";
            }
            else if (sorbmodel == "c2h4*" || sorbmodel == "ethlyene*" || sorbmodel == "ethene*") {
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.666, 0.0, 0.0, 10.011, -0.34772, 1.6304, 52.22317, 3.58174);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, 0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.230, -0.921, 0.0, 1.0079, 0.17386, 0.19, 7.47472, 2.26449);
                system.proto[i].name = "ETH";
            }
            else if (sorbmodel == "c2h4_trappe") {
                addAtomToProto(system,i, "CH2",  "ETH",  "M", -0.665,   0.000,   0.000, 14.02200, -0.000,  0.00000, 85.000,  3.675);
                addAtomToProto(system,i, "CH2",  "ETH",  "M", 0.665,   0.000,   0.000, 14.02200, -0.000,  0.00000, 85.000,  3.675);
            }

            // ETHANE C2H6
            else if (sorbmodel == "c2h6" || sorbmodel == "ethane") {
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897); //  0.00000  0.00000
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.762,   0.000,   0.000, 12.01100, -0.04722,  0.00000, 141.80885,  3.28897);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.0, 0.62069, 2.88406);
                system.proto[i].name = "ETH";
            }
            else if (sorbmodel == "c2h6*" || sorbmodel == "ethane*") {
                addAtomToProto(system,i, "C2G", "ETH", "M", -0.762,0.000, 0.000,12.01100, -0.04722,0.6967,98.63326,3.37151);
                addAtomToProto(system,i, "C2G", "ETH", "M", 0.762,0.000,0.000,12.01100,-0.04722,0.6967,98.63326,  3.37151);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, 1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", -1.156, -0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, 0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, 0.508, -0.879, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                addAtomToProto(system,i, "H2G", "ETH", "M", 1.156, -1.015, 0.0, 1.0079, 0.01574, 0.4758, 2.60236, 2.57302);
                system.proto[i].name = "ETH";
            }
            else if (sorbmodel == "c2h6_trappe") {
                addAtomToProto(system, i, "C2G",  "ETH",  "M", 0.770,   0.000,   0.000,  15.035, 0.0000,  0.00000, 98.00,  3.750);
                addAtomToProto(system, i, "C2G",  "ETH",  "M", -0.770,   0.000,   0.000,  15.035, 0.0000,  0.00000, 98.00,  3.750);
            }
            // WATER H2O (TIP3P)
            else if (sorbmodel == "h2o" || sorbmodel == "water" || sorbmodel == "tip3p") {
                addAtomToProto(system,i, "OXY", "H2O", "M", 0.0, 0.0, 0.0, 16.0, -0.834, 0.0, 76.42, 3.151);
                addAtomToProto(system,i, "HYD", "H2O", "M", -0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                addAtomToProto(system,i, "HYD", "H2O", "M", 0.757, -0.586, 0.0, 1.008, 0.417, 0.0, 0.0, 0.0);
                system.proto[i].name = "H2O";
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


        // finally, show the current proto molecules
        for (int i=0; i<system.proto.size(); i++) {
        printf("Prototype molecule %i has PDBID %i ( name %s ) and has %i atoms\n",i, system.proto[i].PDBID, system.proto[i].name.c_str(), (int)system.proto[i].atoms.size());
            //system.proto.printAll();
            for (int j=0; j<system.proto[i].atoms.size(); j++) {
                system.proto[i].atoms[j].printAll();
            }
        }

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
        system.stats.csp.name = "csp";
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

    }
    
    printf("INPUT: Calculated fugacity for prototype %s = %f atm\n", system.constants.fugacity_single_sorbate.c_str(), system.proto[0].fugacity);

    return;


}
