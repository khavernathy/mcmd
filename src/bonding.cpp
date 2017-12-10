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
        if (num_bonds == 1) return "N_1";
        else if (num_bonds == 2) return "N_3";
        else if (num_bonds == 3) return "N_2";
        else if (num_bonds == 4) return "N_3";
        // account for N_R...
    } else if (name == "O") {
       if (num_bonds == 1) return "O_1";
       else if (num_bonds == 2) return "O_3";
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
    
    int i,j,l,m,p; // i=mol, j,l,m,p are atoms (conventionally IJKL)
    double r, ra, rh; // bond r, angle bond ra, dihedral bond rh
    int local_bonds=0;
    int duplicateFlag=0; //bond dupes
    int duplicateAngleFlag=0; // angle dupes
    int duplicateDihFlag=0; // dihedral dupes
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            local_bonds = 0;
            // for each atom, we find its bonded neighbors by a distance search
            // (by only searching atoms on this same molecule)
            for (l=0; l<system.molecules[i].atoms.size(); l++) {
               if (j==l) continue; // don't do self-atom
               double* distances = getDistanceXYZ(system, i,j,i,l);
               r = distances[3];
               
               // avoid H-H bonds, except when H2 is the molecule.
               if (r < system.constants.bondlength
                    && (!((system.molecules[i].atoms[j].name == "H" &&
                        system.molecules[i].atoms[l].name == "H") &&
                        system.molecules[i].atoms.size() != 2))) {
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
                    
                        // avoid H-H bonds, except when H2 is the molecule.
                        if (ra < system.constants.bondlength
                            && (!((system.molecules[i].atoms[m].name == "H" &&
                            system.molecules[i].atoms[l].name == "H") &&
                            system.molecules[i].atoms.size() != 2))) {

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

                            // now loop for dihedrals
                            for (p=0; p<system.molecules[i].atoms.size(); p++) {
                                double * distancesh = getDistanceXYZ(system,i,m,i,p);
                                rh = distancesh[3];
                                if (j==l) continue; // don't do AABC
                                if (j==p) continue; // don't do ABCA
                                if (l==m) continue; // don't do ABBC
                                if (m==p) continue; // don't do ABCC
                                if (j==m) continue; // don't do ABAC
                                if (l==p) continue; // don't do ABCB

                                // avoid H-H bonds, except when H2 is the molecule.
                                if (rh < system.constants.bondlength
                                    && (!((system.molecules[i].atoms[p].name == "H" &&
                                    system.molecules[i].atoms[m].name == "H") &&
                                    system.molecules[i].atoms.size() != 2))) {
                                    // check duplicate dihedral
                                    duplicateDihFlag = 0;
                                    for (int n=0;n<system.constants.uniqueDihedrals.size();n++) {
                                        if ((system.constants.uniqueDihedrals[n].mol==i &&
                                            system.constants.uniqueDihedrals[n].atom1==p &&
                                            system.constants.uniqueDihedrals[n].atom2==m &&
                                            system.constants.uniqueDihedrals[n].atom3==l &&
                                            system.constants.uniqueDihedrals[n].atom4==j) ||
                                            (system.constants.uniqueDihedrals[n].mol==i &&
                                            system.constants.uniqueDihedrals[n].atom1==p &&
                                            system.constants.uniqueDihedrals[n].atom2==m &&
                                            system.constants.uniqueDihedrals[n].atom3==l &&
                                            system.constants.uniqueDihedrals[n].atom4==j )) {

                                                duplicateDihFlag = 1;break;
                                        }
                                    }
                                    // this dihedral is unique
                                    if (!duplicateDihFlag) {
                                        Constants::UniqueDihedral tmp; tmp.mol=i; tmp.atom1=j;
                                        tmp.atom2=l; tmp.atom3=m; tmp.atom4=p;
                                        system.constants.uniqueDihedrals.push_back(tmp);
                                    }
                                }

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

    // check for unphysical H-H bonds
        for (int n=0;n<system.constants.uniqueBonds.size();n++) {
            int mol=system.constants.uniqueBonds[n].mol;
            int a1=system.constants.uniqueBonds[n].atom1;
            int a2=system.constants.uniqueBonds[n].atom2;
            if (system.molecules[mol].atoms[a1].name == "H" &&
                system.molecules[mol].atoms[a2].name == "H") {
                
                printf("H-H bond: %i %i (bondid %i)\n", a1,a2,n);
                // see if this H-H is not H2 molecule.
                if (system.molecules[mol].atoms[a1].bonds.size() > 1 ||
                    system.molecules[mol].atoms[a2].bonds.size() > 1) {

                }
            }

        }


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


double get_Vjk(double vj, double vk) {
    // sp3 -- sp3
    return sqrt(vj*vk);
}

double get_dihedral_angle(System &system, int mol, int i, int j, int k, int l) {
    // current dihedral angle for i--j--k--l atoms
    /*
     *          L    <- 0 degree dihedral (in plane of screen)..
     *         /
     *    J---K
     *   /
     *  I
     *
     * */
    // 4 components of a plane are A,B,C,D in Ax + By + Cz + D = 0
    // The D parameter is irrelevant here so we just need a vector holding A,B,C
    double Plane1[3]; // build from atoms i,j,k
    double Plane2[3]; // build from atoms j,k,l
    double v1a[3], v1b[3]; // for plane 1
    double v2a[3], v2b[3]; // for plane 2

    for (int n=0; n<3; n++) {
        v1a[n] = system.molecules[mol].atoms[j].pos[n] - system.molecules[mol].atoms[i].pos[n];
        v1b[n] = system.molecules[mol].atoms[k].pos[n] - system.molecules[mol].atoms[i].pos[n];
        
        v2a[n] = system.molecules[mol].atoms[k].pos[n] - system.molecules[mol].atoms[j].pos[n];
        v2b[n] = system.molecules[mol].atoms[l].pos[n] - system.molecules[mol].atoms[j].pos[n];
    }

    double* norm1 = crossprod(v1a, v1b);
    double* norm2 = crossprod(v2a, v2b);
    
    for (int n=0;n<3;n++) {
        Plane1[n] = norm1[n];
        Plane2[n] = norm2[n];
    }
    //Plane1[3] = -dddotprod(norm1, system.molecules[mol].atoms[i].pos);
    //Plane2[3] = -dddotprod(norm2, system.molecules[mol].atoms[j].pos);

    // both planes done; now get angle
    const double dotplanes = dddotprod(Plane1,Plane2);
    const double mag1 = sqrt(dddotprod(Plane1,Plane1));
    const double mag2 = sqrt(dddotprod(Plane2,Plane2));

    return acos(dotplanes/(mag1*mag2)); 
}

// get the total potential from torsions
// via simple Fourier small cosine expansion
double torsions_energy(System &system) {
    double potential=0;
    double vjk, vj, vk, n, dihedral, phi_ijkl; // n is periodicity (integer quantity)
    const double deg2rad = M_PI/180.0;
    int i,j,l,m,p; // molecule i, atoms (j,l,m and p)
    for (int it=0; it<system.constants.uniqueDihedrals.size(); it++) {
        i = system.constants.uniqueDihedrals[it].mol;
        j = system.constants.uniqueDihedrals[it].atom1;
        l = system.constants.uniqueDihedrals[it].atom2;
        m = system.constants.uniqueDihedrals[it].atom3;
        p = system.constants.uniqueDihedrals[it].atom4;

        vj = system.constants.UFF_torsions[system.molecules[i].atoms[l].UFFlabel.c_str()];
        vk = system.constants.UFF_torsions[system.molecules[i].atoms[m].UFFlabel.c_str()];
        vjk = get_Vjk(vj, vk);
        dihedral = get_dihedral_angle(system, i, j,l,m,p);

        phi_ijkl = 60.0*deg2rad; // ...
        n = 3.0; // ..

        potential += 0.5*vjk*(1.0 - cos(n*phi_ijkl)*cos(n*dihedral));//0.5*vjk;
    }

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

double simple_r(double xi, double xj, double yi, double yj, double zi, double zj) {
    return sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj));
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_gradient(System &system) {
    const double deg2rad = M_PI/180.0;
    int i,j,l,m;
    double rij, rjk, rik, K_ijk, C0, C1, C2, theta_ijk; // angle-bend params
    double angle; // the actual angle IJK
    double grad;
    double cos1, cos2, B, C, xi, yi, zi, xj, yj, zj, xk, yk, zk, dij, djk; 
        // cos-derivative terms in the gradient, + other terms.


    for (int it=0; it<system.constants.uniqueAngles.size(); it++) {

        i = system.constants.uniqueAngles[it].mol;
        j = system.constants.uniqueAngles[it].atom1;
        l = system.constants.uniqueAngles[it].atom2;
        m = system.constants.uniqueAngles[it].atom3;
        
            
        rij = get_rij(system,i,j,i,l); // in Angstroms
        // force field theta is that of the middle atom "l" in j-l-m
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

        // compute the gradient for all angle components (xyz for 3 atoms = 9)
        // gradient is [dE/dx_i...] which ends up as a sum of two cosine derivs (d/dx C0 term -> 0)
        // this gets funky b/c gradient looks different for each direction x,y,z
        // just do dE/dx_i now.
        xi = system.molecules[i].atoms[j].pos[0]; 
        yi = system.molecules[i].atoms[j].pos[1];
        zi = system.molecules[i].atoms[j].pos[2];
        xj = system.molecules[i].atoms[l].pos[0];
        yj = system.molecules[i].atoms[l].pos[1];
        zj = system.molecules[i].atoms[l].pos[2];
        xk = system.molecules[i].atoms[m].pos[0];
        yk = system.molecules[i].atoms[m].pos[1];
        zk = system.molecules[i].atoms[m].pos[2];

        B = (yi-yj)*(yk-yj) + (zi-zj)*(zk-zj);
        C = (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);

        dij = simple_r(xi,xj,yi,yj,zi,zj);
        djk = simple_r(xj,xk,yj,yk,zj,zk);

        cos1 = ((B*(xj-xi) + C*(xk-xj))/(djk*dij*dij*dij));
        cos2 = 2*(((xk-xj)/(djk*dij)) - (( (xi-xj)*(B + (xk-xj)*(xi-xj)))/(djk*dij*dij*dij))) * sin(2*acos((B + (xk-xj)*(xi-xj))/(djk*dij))) / sqrt(1 - pow((B + (xk-xj)*(xi-xj)),2)/(djk*djk*dij*dij));

        grad = K_ijk*(C1*cos1 + C2*cos2);
        system.molecules[i].atoms[j].energy_grad[0] += grad; // grad;
        

        double POT=K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)); // in kcal/mol
            
    }

    return 0.; // in kcal/mol
} // end angle bend energy


double totalBondedEnergy(System &system) {
    return stretch_energy(system) + angle_bend_energy(system) + torsions_energy(system);
}


