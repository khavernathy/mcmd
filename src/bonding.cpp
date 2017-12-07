/* Douglas M Franz
 * Space group, USF, 2017
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

// function to determine UFF atom-type based on
// name of element and number of bonds
string getUFFlabel(System &system, string name, int num_bonds) {
    // starting with just the organic-y atom-types.
    if (name == "H") {
        return "H_"; // assume it's never H_b (borane hydrogen)
    } else if (name == "B") {
        if (num_bonds == 3) return "B_2";
        else if (num_bonds == 4) return "B_3";
    } else if (name == "C") {
        if (num_bonds == 2) return "C_1";
        else if (num_bonds == 3) return "C_2";
        else if (num_bonds == 4) return "C_3";
        // need to dynamically account for resonant C_R too...
    } else if (name == "N") {
        if (num_bonds == 1 || num_bonds == 2) return "N_1";
        else if (num_bonds == 3) return "N_2";
        else if (num_bonds == 4) return "N_3";
        // account for N_R...
    } else if (name == "O") {
       if (num_bonds == 1 || num_bonds == 2) return "O_1";
       else if (num_bonds == 3) return "O_2";
       else if (num_bonds == 4) return "O_3";
       // account for O_R...
    } else if (name == "F") {
        return "F_";
    } else if (name == "P") {
        // weird geometries
    } else if (name == "S") {
        // weird geometries
    } else if (name == "Cl") {
        return "Cl";
    } else if (name == "Br") {
        return "Br";
    } else if (name == "I") {
        return "I_";
    }
        
    return "NOTFOUND";
}

// function to find all bonds (and angles) for all atoms.
void findBonds(System &system) {
    
    int i,j,l,m;
    double r, ra; // bond r, angle bond r
    int local_bonds=0;
    int duplicateFlag=0; //bond dupes
    int duplicateAngleFlag=0; // angle dupes
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            local_bonds = 0;
            // for each atom, we find its bonded neighbors by a distance search
            // (by only searching atoms on this same molecule)
            for (l=0; l<system.molecules[i].atoms.size(); l++) {
               if (j==l) continue; // don't do self-atom
               double* distances = getDistanceXYZ(system, i,j,i,l);
               r = distances[3];
               if (r < system.constants.bondlength) {
                    local_bonds++;

                    system.molecules[i].atoms[j].bonds.push_back(l);

                    // check for duplicate bond.
                    duplicateFlag=0;
                    for (int n=0;n<system.constants.uniqueBonds.size();n++) {
                        if (system.constants.uniqueBonds[n].mol == i &&
                            system.constants.uniqueBonds[n].atom1 == l &&
                            system.constants.uniqueBonds[n].atom2 == j) {
                            
                            duplicateFlag=1; break;
                        }
                    }
                    if (!duplicateFlag) {
                        // this bond is unique
                        Constants::UniqueBond tmp; tmp.mol=i; tmp.atom1=j; tmp.atom2=l;         
                        system.constants.uniqueBonds.push_back(tmp);
                    }

                    // now loop for angles
                    for (m=0; m<system.molecules[i].atoms.size(); m++) {
                        if (j==m) continue; // don't do ABA angles (0 degrees)
                        if (l==m) continue; // don't do ABB angles (0 degrees)
                        double* distancesa = getDistanceXYZ(system, i,l,i,m);
                        ra = distancesa[3];                

                        if (ra < system.constants.bondlength) {

                            // check for duplicate angles
                            duplicateAngleFlag = 0;
                            for (int n=0;n<system.constants.uniqueAngles.size();n++) {
                                if (system.constants.uniqueAngles[n].mol == i &&
                                    system.constants.uniqueAngles[n].atom1 == m &&
                                    system.constants.uniqueAngles[n].atom2 == l &&
                                    system.constants.uniqueAngles[n].atom3 == j) {
                                    
                                    duplicateAngleFlag=1; break;

                                }   
                            }
                            if (!duplicateAngleFlag) {
                                // this angle is unique
                                Constants::UniqueAngle tmp; tmp.mol=i; tmp.atom1=j; tmp.atom2=l; tmp.atom3=m;
                                system.constants.uniqueAngles.push_back(tmp);
                            }
                        }
                    }
                    
               } // end if r < bond-length      
            } // end pair (i,j) -- (i,l)  
         
            // based on the total number of bonds to this atom, 
            // determine the atom-type from UFF.
            system.molecules[i].atoms[j].UFFlabel = getUFFlabel(system, system.molecules[i].atoms[j].name, system.molecules[i].atoms[j].bonds.size()); 

        } // end j
    } // end i
}

/*
 * Essentially everything below comes from the parameters and
 * functional forms prescribed in UFF, via
 * Rappe et al.
 * J. Am. Chem. Soc. Vol. 114, No. 25, 1992
 * */

// provide rij parameter given needed information (atom IDs)
double get_rij(System &system, int i, int j, int k, int l) {
    const double ri = system.constants.UFF_bonds[system.molecules[i].atoms[j].UFFlabel.c_str()];
    const double rj = system.constants.UFF_bonds[system.molecules[k].atoms[l].UFFlabel.c_str()];
    const double Xi = system.constants.UFF_electroneg[system.molecules[i].atoms[j].name.c_str()];
    const double Xj = system.constants.UFF_electroneg[system.molecules[k].atoms[l].name.c_str()];
    const double BO = 1.0; // assume all single bonds for now.
    const double lambda = 0.1332; // Bond-order correction parameter

    const double rBO = -lambda*(ri + rj)*log(BO);
    const double Xij = sqrt(Xi) - sqrt(Xj);
    const double rEN = ri*rj*(Xij*Xij)/(Xi*ri + Xj*rj);
    return ri + rj + rBO + rEN; // in Angstroms
}

// get Force constant kij for a bond.
double get_kij(System &system, int i, int j, int k, int l, double rij) {
    const double Zi = system.constants.UFF_Z[system.molecules[i].atoms[j].UFFlabel.c_str()];
    const double Zj = system.constants.UFF_Z[system.molecules[k].atoms[l].UFFlabel.c_str()];
    return 664.12*Zi*Zj/(rij*rij*rij); // in kcal/molA^2, as-is
}

// get the total potential from bond stretches
// via the Morse potential
double stretch_energy(System &system) {
    double potential = 0;
    double alpha,Dij,kij,rij; // bond params
    double mainterm; // main chunk of potential, to be squared..
    int i,j,l; // atom indices
    double r; // actual, current distance for pair.
    /* ============================ */
    double BO = 1.0; // assume single bonds for now!!!
    /* ============================ */

    // loop through bonds of this atom.
    for (int it=0; it<system.constants.uniqueBonds.size(); it++) {
        i = system.constants.uniqueBonds[it].mol;
        j = system.constants.uniqueBonds[it].atom1;
        l = system.constants.uniqueBonds[it].atom2;

        rij = get_rij(system,i,j,i,l); // in Angstroms
        kij = get_kij(system,i,j,i,l, rij); // in kcal mol^-1 A^-2
//          printf("rij = %f; kij = %f\n", rij, kij);

        Dij = BO*70.0; // in kcal/mol
        alpha = sqrt(0.5*kij/Dij); // in 1/A

        double* distances = getDistanceXYZ(system, i,j,i,l);
        r = distances[3];
        mainterm = exp(-alpha*(r-rij)) - 1.0; // unitless
        potential += Dij*(mainterm*mainterm); // in kcal/mol
    }

    return potential; // in kcal/mol
}

// get angle-force parameter for IJK triplet 
double get_Kijk(System &system, double rij, double rjk, double rik, double Zi, double Zk, double angle) {
    const double beta = 664.12/(rij*rjk);
    //printf("angle = %f\n", angle*180.0/M_PI);
    //printf("K_ijk = %f\n",  beta*Zi*Zk/(rik*rik*rik*rik*rik) * rij*rjk*(3.0*rij*rjk*(1-cos(angle)*cos(angle)) - rik*rik*cos(angle)));
    return  beta*Zi*Zk/(rik*rik*rik*rik*rik) * rij*rjk*(3.0*rij*rjk*(1-cos(angle)*cos(angle)) - rik*rik*cos(angle));
}

// get the angle ABC where B is center atom, on molecule i
double get_angle(System &system, int i, int A, int B, int C) {
    // https://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space
    double AB[3] = {0,0,0};
    double BC[3] = {0,0,0};

    for (int n=0;n<3;n++) {
        AB[n] = system.molecules[i].atoms[A].pos[n] - system.molecules[i].atoms[B].pos[n];
        BC[n] = system.molecules[i].atoms[C].pos[n] - system.molecules[i].atoms[B].pos[n];
    }
    
    const double dotprod = dddotprod(AB,BC);
    const double ABm = sqrt(dddotprod(AB,AB));
    const double BCm = sqrt(dddotprod(BC,BC));
    
    return acos(dotprod/(ABm*BCm)); // returns in radians
}

// get r_ik, a parameter for angle bends, ** different from r_ij (and r_jk) **
double get_rik(System &system, double rij, double rjk, double angle) {
    // angle is in radians
    return sqrt(rij*rij + rjk*rjk - 2*rij*rjk*cos(angle));
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_energy(System &system) {
    double potential=0;
    const double deg2rad = M_PI/180.0;
    int i,j,l,m;
    double rij, rjk, rik, K_ijk, C0, C1, C2, theta_ijk; // angle-bend params
    double angle; // the actual angle IJK
    
    for (int it=0; it<system.constants.uniqueAngles.size(); it++) {
        i = system.constants.uniqueAngles[it].mol;
        j = system.constants.uniqueAngles[it].atom1;
        l = system.constants.uniqueAngles[it].atom2;
        m = system.constants.uniqueAngles[it].atom3;


        rij = get_rij(system,i,j,i,l); // in Angstroms
        theta_ijk = deg2rad*system.constants.UFF_angles[system.molecules[i].atoms[l].UFFlabel.c_str()]; // in rads
        C2 = 1.0/(4.0*sin(theta_ijk)*sin(theta_ijk));      // 1/rad^2
        C1 = -4.0*C2*cos(theta_ijk);                       // 1/rad
        C0 = C2*(2.0*cos(theta_ijk)*cos(theta_ijk) + 1.0); // 1
        //printf("theta_0 = %f\n", theta_ijk/deg2rad);
            
        angle = get_angle(system, i, j, l, m);  
        //printf("Angle %i %i %i = %f; real angle = %f\n", j,l,m,theta_ijk/deg2rad, angle/deg2rad);
        rjk = get_rij(system,i,l,i,m);
        rik = get_rik(system, rij, rjk, angle); // r_ik (A-C) is computed differently than r_ij (A-B) and r_jk (B-C)
        K_ijk = get_Kijk(system, rij, rjk, rik, system.constants.UFF_Z[system.molecules[i].atoms[j].UFFlabel.c_str()], system.constants.UFF_Z[system.molecules[i].atoms[m].UFFlabel.c_str()], theta_ijk);      
        //printf("K_ijk = %f \n", K_ijk);
        //printf("rij = %f; rjk = %f; rik = %f\n", rij, rjk, rik);

        potential += K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)); // in kcal/mol
    
    }
    return potential; // in kcal/mol
} // end angle bend energy

// get the total potential from torsions
// via simple Fourier small cosine expansion
double torsions_energy(System &system) {
    double potential=0;
    const double deg2rad = M_PI/180.0;
    int i,j,l,m,n; // molecule i, atoms (j,l,m and n)
    



    return potential; // in kcal/mol
}





// Morse potential gradient for all bonded atoms, to minimize energy
double morse_gradient(System &system) {
    // x,y,z is the point of interest for this gradient
    // ij tells us whether it's the 1st atom or 2nd within definition of delta x (+1 or -1)
    double alpha,Dij,kij,rij; // bond params
    int i,j,l; // atom indices
    double r; // actual, current distance for pair.
    /* ============================ */
    double BO = 1.0; // assume single bonds for now!!!
    /* ============================ */
    double prefactor,grad,delta;
        // typical gradient elements (e.g. dE/dx_i) are ~10^2 in these units.


            for (int it=0; it<system.constants.uniqueBonds.size(); it++) {
                i = system.constants.uniqueBonds[it].mol;
                j = system.constants.uniqueBonds[it].atom1;
                l = system.constants.uniqueBonds[it].atom2;

                rij = get_rij(system,i,j,i,l); // in Angstroms
                kij = get_kij(system,i,j,i,l, rij); // in kcal mol^-1 A^-2
      //          printf("rij = %f; kij = %f\n", rij, kij);

                Dij = BO*70.0; // in kcal/mol
                alpha = sqrt(0.5*kij/Dij); // in 1/A

                double* distances = getDistanceXYZ(system, i,j,i,l);
                r = distances[3];
    
                // gradient for a single bond is 6D (3D on each atom, 1 for each D.O.F.)
                prefactor = 2*alpha*Dij*exp(alpha*(rij-r))/r;
                for (int n=0;n<3;n++) {
                    delta = system.molecules[i].atoms[j].pos[n] - system.molecules[i].atoms[l].pos[n];
                    grad = prefactor * delta;
                    grad *= (1 - exp(alpha*(rij-r)));
                    system.molecules[i].atoms[j].energy_grad[n] += grad;
                    system.molecules[i].atoms[l].energy_grad[n] -= grad;
                }
                // xj, yj, zj
                // since gradient of the other atom is just minus the other, we apply a Newton-pair style thing above
                // instead of recomputing.
                /*
                    for (int n=0;n<3;n++) {
                    delta = system.molecules[i].atoms[l].pos[n] - system.molecules[i].atoms[j].pos[n];
                    grad = prefactor * delta;
                    grad *= (1 - exp(alpha*(rij-r)));
                    printf("%f\n", grad);
                    // move the atom position element in direction of energy minimum
                    system.molecules[i].atoms[l].pos[n] += grad*move_factor;
                }
                */

            }

    return 0; //.5*potential; // in kcal/mol
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_gradient(System &system) {
    const double deg2rad = M_PI/180.0;
    int i,j,l,m;
    double rij, rjk, rik, K_ijk, C0, C1, C2, theta_ijk; // angle-bend params
    double angle; // the actual angle IJK
    
    for (int it=0; it<system.constants.uniqueAngles.size(); it++) {

        i = system.constants.uniqueAngles[it].mol;
        j = system.constants.uniqueAngles[it].atom1;
        l = system.constants.uniqueAngles[it].atom2;
        m = system.constants.uniqueAngles[it].atom3;
        
            
        rij = get_rij(system,i,j,i,l); // in Angstroms
        theta_ijk = deg2rad*system.constants.UFF_angles[system.molecules[i].atoms[l].UFFlabel.c_str()]; // in rads
        C2 = 1.0/(4.0*sin(theta_ijk)*sin(theta_ijk));      // 1/rad^2
        C1 = -4.0*C2*cos(theta_ijk);                       // 1/rad
        C0 = C2*(2.0*cos(theta_ijk)*cos(theta_ijk) + 1.0); // 1
        //printf("theta_0 = %f\n", theta_ijk/deg2rad);
            
        angle = get_angle(system, i, j, l, m);  
        //printf("Angle %i %i %i = %f; real angle = %f\n", j,l,m,theta_ijk/deg2rad, angle/deg2rad);
        rjk = get_rij(system,i,l,i,m);
        rik = get_rik(system, rij, rjk, angle); // r_ik (A-C) is computed differently than r_ij (A-B) and r_jk (B-C)
        K_ijk = get_Kijk(system, rij, rjk, rik, system.constants.UFF_Z[system.molecules[i].atoms[j].UFFlabel.c_str()], system.constants.UFF_Z[system.molecules[i].atoms[m].UFFlabel.c_str()], theta_ijk);      
        //printf("K_ijk = %f \n", K_ijk);
        //printf("rij = %f; rjk = %f; rik = %f\n", rij, rjk, rik);



         // HALF BECAUSE THE ANGLES WILL BE DOUBLE COUNTED.
         double POT=0; POT += 0.5 * K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)); // in kcal/mol
            
    }

    return 0.; // in kcal/mol
} // end angle bend energy


