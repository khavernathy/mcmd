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


/* GET BOX LIMIT COORDINATE FOR PBC */
double getBoxLimit(System &system, int planeid, double x, double y, double z) {

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
void checkInTheBox(System &system, int i) {

// first the easy systems, alpha=beta=gamma=90
// saves some time because box limits need not be recomputed
if (system.pbc.alpha == 90 && system.pbc.beta == 90 && system.pbc.gamma == 90) {
    if (system.constants.md_mode == "molecular") {
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
    } // end if molecular
    else if (system.constants.md_mode == "atomic") {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {        

        if (system.molecules[i].atoms[j].pos[0]  > system.pbc.x_max) {
                    system.molecules[i].atoms[j].pos[0] -= system.pbc.x_length;
        }
        else if (system.molecules[i].atoms[j].pos[0] < system.pbc.x_min) {
                    system.molecules[i].atoms[j].pos[0] += system.pbc.x_length;
        }

        if (system.molecules[i].atoms[j].pos[1] > system.pbc.y_max) {
                    system.molecules[i].atoms[j].pos[1] -= system.pbc.y_length;
        }
        else if (system.molecules[i].atoms[j].pos[1] < system.pbc.y_min) {
                    system.molecules[i].atoms[j].pos[1] += system.pbc.y_length;
        }

        if (system.molecules[i].atoms[j].pos[2]  > system.pbc.z_max) {
                    system.molecules[i].atoms[j].pos[2] -= system.pbc.z_length;
        }
        else if (system.molecules[i].atoms[j].pos[2] < system.pbc.z_min) {
                    system.molecules[i].atoms[j].pos[2] += system.pbc.z_length;
        } 
    } // end atom loop
    } // end if atomic
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
    if (system.constants.md_mode == "molecular") {

        if (system.molecules[i].com[0] < box_limit[0]) { // left of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new x coordinate
                system.molecules[i].atoms[k].pos[0] += system.pbc.x_length;
                // set new y based on x
                double ydist = system.molecules[i].atoms[k].pos[1] - box_limit[2];
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydist;
                // set new z based on new xy
                double zdist = system.molecules[i].atoms[k].pos[2] - box_limit[4];
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdist;
            }
        }
    
        if (system.molecules[i].com[0] > box_limit[1]) { // right of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new x coordinate
                system.molecules[i].atoms[k].pos[0] -= system.pbc.x_length;
                // set new y based on x
                double ydist = system.molecules[i].atoms[k].pos[1] - box_limit[2];
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydist;
                // set new z based on new xy
                double zdist = system.molecules[i].atoms[k].pos[2] - box_limit[4];
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdist;
            }

        }

        if (system.molecules[i].com[1] < box_limit[2]) { // below box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new y coordinate
                system.molecules[i].atoms[k].pos[1] += system.pbc.y_length;
                // set new x based on new y
                double xdist = system.molecules[i].atoms[k].pos[0] - box_limit[0];
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdist;
                // set new z based on new xy
                double zdist = system.molecules[i].atoms[k].pos[2] - box_limit[4];
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdist;
            }

        }

        if (system.molecules[i].com[1] > box_limit[3]) { // above box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new y coordinate
                system.molecules[i].atoms[k].pos[1] -= system.pbc.y_length;
                // set new x based on new y
                double xdist = system.molecules[i].atoms[k].pos[0] - box_limit[0];
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdist;
                // set new z based on new xy
                double zdist = system.molecules[i].atoms[k].pos[2] - box_limit[4];
                newlimit[4] = getBoxLimit(system, 4, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[2] = newlimit[4] + zdist;
            }
        }

        if (system.molecules[i].com[2] < box_limit[4]) { // behind box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new z coordinate
                system.molecules[i].atoms[k].pos[2] += system.pbc.z_length;
                // set new y based on new z
                double ydist = system.molecules[i].atoms[k].pos[1] - box_limit[2];
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydist;
                // set new x based on new yz
                double xdist = system.molecules[i].atoms[k].pos[0] - box_limit[0];
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdist;
            }

        }

        if (system.molecules[i].com[2] > box_limit[5]) { // in front of box
            for (int k=0; k<system.molecules[i].atoms.size(); k++) {
                // set new z coordinate
                system.molecules[i].atoms[k].pos[2] -= system.pbc.z_length;
                // set new y based on new z
                
                double ydist = system.molecules[i].atoms[k].pos[1] - box_limit[2];
                newlimit[2] = getBoxLimit(system, 2, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[1] = newlimit[2] + ydist;
                // set new x based on new yz
                double xdist = system.molecules[i].atoms[k].pos[0] - box_limit[0];
                newlimit[0] = getBoxLimit(system, 0, system.molecules[i].atoms[k].pos[0], system.molecules[i].atoms[k].pos[1], system.molecules[i].atoms[k].pos[2]);
                system.molecules[i].atoms[k].pos[0] = newlimit[0] + xdist;
            }
        }

    } else if (system.constants.md_mode == "atomic") {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
        int a=0;  // i'm too lazy to bother with an atomic algorithm because the prog really should just be molecules.
        }
    } // end if atomic     

} 
} // end pbc function


void moleculePrintout(System &system) {
    // CONFIRM ATOMS AND MOLECULES PRINTOUT
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
    printf("Prototype molecule has PDBID %i ( name %s ) and has %i atoms\n", system.proto.PDBID, system.proto.name.c_str(), (int)system.proto.atoms.size());
    /*for (int i=0; i<system.proto.atoms.size(); i++) {
        system.proto.atoms[i].printAll();
    } */    
    }
    
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
   
    // SUMS
        // ENERGY
    system.last.current_rd_sum = system.stats.current_rd_sum; 
        system.last.current_lj_lrc_sum = system.stats.current_lj_lrc_sum;
        system.last.current_lj_self_lrc_sum = system.stats.current_lj_self_lrc_sum;
        system.last.current_lj_sum = system.stats.current_lj_sum;
    system.last.current_es_sum = system.stats.current_es_sum;
        system.last.current_coulombic_self_sum = system.stats.current_coulombic_self_sum;
        system.last.current_coulombic_real_sum = system.stats.current_coulombic_real_sum;
        system.last.current_coulombic_reciprocal_sum = system.stats.current_coulombic_reciprocal_sum;
    system.last.current_polar_sum = system.stats.current_polar_sum;
    system.last.current_energy_sum = system.stats.current_energy_sum;
        // QST
    system.last.current_Nsq_sum = system.stats.current_Nsq_sum;
    system.last.current_NU_sum = system.stats.current_NU_sum;
    system.last.current_qst_sum = system.stats.current_qst_sum; 
        // VOLUME
    system.last.current_volume_sum = system.stats.current_volume_sum;    
        // DENSITY
    system.last.current_density_sum = system.stats.current_density_sum;
        // Z
    system.last.current_z_sum = system.stats.current_z_sum;    


    // VALUES
        // ENERGY
    system.last.rdU = system.stats.rdU;
        system.last.lj_lrc = system.stats.lj_lrc;
        system.last.lj_self_lrc = system.stats.lj_self_lrc;
        system.last.lj = system.stats.lj;
    system.last.esU = system.stats.esU;
        system.last.coulombic_self = system.stats.coulombic_self;
        system.last.coulombic_real = system.stats.coulombic_real;
        system.last.coulombic_reciprocal = system.stats.coulombic_reciprocal;
    system.last.polarU = system.stats.polarU;
    system.last.totalU = system.stats.totalU;
        // QST
    system.last.qst = system.stats.qst;    
        // VOLUME
    system.last.volume = system.stats.volume;
        // DENSITY
    system.last.density = system.stats.density;
        // CHEMICAL POTENTIAL
    system.last.chemical_potential = system.stats.chemical_potential;
        // Z
    system.last.Z = system.stats.Z;


}

void revertToCheckpoint(System &system) {
    // reverts variables of interest if needed
    
    // SUMS
        // ENERGY
    system.stats.current_rd_sum = system.last.current_rd_sum + system.last.rdU;
        system.stats.current_lj_lrc_sum = system.last.current_lj_lrc_sum + system.last.lj_lrc;
        system.stats.current_lj_self_lrc_sum = system.last.current_lj_self_lrc_sum + system.last.lj_self_lrc;
        system.stats.current_lj_sum = system.last.current_lj_sum + system.last.lj;
    system.stats.current_es_sum = system.last.current_es_sum + system.last.esU;
        system.stats.current_coulombic_self_sum = system.last.current_coulombic_self_sum + system.last.coulombic_self;
        system.stats.current_coulombic_real_sum = system.last.current_coulombic_real_sum + system.last.coulombic_real;
        system.stats.current_coulombic_reciprocal_sum = system.last.current_coulombic_reciprocal_sum + system.last.coulombic_reciprocal;
    system.stats.current_polar_sum = system.last.current_polar_sum + system.last.polarU;
    system.stats.current_energy_sum = system.last.current_energy_sum + system.last.totalU;
    /**    // QST 
    system.stats.current_Nsq_sum = system.last.current_Nsq_sum;
    system.stats.current_NU_sum = system.last.current_NU_sum;
    system.stats.current_qst_sum = system.last.current_qst_sum;
        // VOLUME
    system.stats.current_volume_sum = system.last.current_volume_sum;
        // DENSITY
    system.stats.current_density_sum = system.last.current_density_sum;
        // Z
    system.stats.current_z_sum = system.last.current_z_sum;
    */

    // VALUES
        // ENERGY
    system.stats.rdU = system.last.rdU;
        system.stats.lj_lrc = system.last.lj_lrc;
        system.stats.lj_self_lrc = system.last.lj_self_lrc;
        system.stats.lj = system.last.lj;
    system.stats.esU = system.last.esU;
        system.stats.coulombic_self = system.last.coulombic_self;
        system.stats.coulombic_real = system.last.coulombic_real;
        system.stats.coulombic_reciprocal = system.last.coulombic_reciprocal;
    system.stats.polarU = system.last.polarU;
    system.stats.totalU = system.last.totalU;    
    /*    // QST
    system.stats.qst = system.last.qst;
        // VOLUME
    system.stats.volume = system.last.volume;
        // DENSITY
    system.stats.density = system.last.density;
        // CHEMICAL POTENTIAL
    system.stats.chemical_potential = system.last.chemical_potential;
        // Z
    system.stats.z = system.last.z;
    */
}
