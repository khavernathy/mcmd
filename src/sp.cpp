#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

/*
 * Functions for single-point energy of a given molecule
*/

double avgr(double r) {
    // returns sqrt(<r^2>) of a sphere with radius r.
    return sqrt(3./5.*r*r);
}

void singlePointEnergy(System &system) {

    int i,j; // atom id
    int n; // dimensions
    double q, x, y, z, comx, comy, comz, ax,ay,az,a2; // temporary charge holder, pos. vector
    double multipole0 = 0; // charge sum
    double multipole1[3] = {0,0,0}; // dipole moment
    double dipole_magnitude = 0;
    double multipole2[6] = {0,0,0,0,0,0}; // quadrupole moment
        // 0: xx; 1: yy: 2: zz; 3: xy; 4: xz; 5: yz
    double o[18] = {0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0, 0,0,0}; // octapole moment, writing as "o" b/c it's faster.
        // 0: xxx; 1: xyy; 2: xzz; 3: yxx; 4: yyy; 5: yzz; 6: zxx; 7: zyy; 8: zzz  (trace)
        // 9: xxy; 10: xxz; 11: xyz; 12: yxy; 13: yxz; 14: yyz; 15: zxy; 16: zxz; 17: zyz

    printf("Input atoms: %s\n=============================================\n\n", system.constants.atom_file.c_str());

    system.molecules[0].calc_center_of_mass();

    // multiple 0 (total charge)
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        q = system.molecules[0].atoms[i].C/system.constants.E2REDUCED;
        // charge sum
        multipole0 += q;
        x = system.molecules[0].atoms[i].pos[0];
        y = system.molecules[0].atoms[i].pos[1];
        z = system.molecules[0].atoms[i].pos[2];

        comx = system.molecules[0].com[0];
        comy = system.molecules[0].com[1];
        comz = system.molecules[0].com[2];

        ax = x-comx; ay = y-comy; az=z-comz;
        a2 = ax*ax + ay*ay + az*az;

        // dipole 
        
        multipole1[0] += q*ax;
        multipole1[1] += q*ay;
        multipole1[2] += q*az;
        /*    
        multipole1[0] += q*x;
        multipole1[1] += q*y;
        multipole1[2] += q*z;
        */

        // quadrupole
        multipole2[0] += q*(1.5*(ax*ax) - 0.5*(a2));
        multipole2[1] += q*(1.5*(ay*ay) - 0.5*(a2));
        multipole2[2] += q*(1.5*(az*az) - 0.5*(a2));
        multipole2[3] += q*1.5*ax*ay;
        multipole2[4] += q*1.5*ax*az;
        multipole2[5] += q*1.5*ay*az;
    
        // octapole                        d_yz   d_xz   d_xy   // delta f'ns
        o[0]  += q*(2.5*ax*ax*ax - 0.5*a2*(ax*1 + ax*1 + ax*1)); 
        o[1]  += q*(2.5*ax*ay*ay - 0.5*a2*(ax*1 + ay*0 + ay*0));
        o[2]  += q*(2.5*ax*az*az - 0.5*a2*(ax*1 + az*0 + az*0));
        o[3]  += q*(2.5*ay*ax*ax - 0.5*a2*(ay*1 + ax*0 + ax*0));
        o[4]  += q*(2.5*ay*ay*ay - 0.5*a2*(ay*1 + ay*1 + ay*1));
        o[5]  += q*(2.5*ay*az*az - 0.5*a2*(ay*1 + az*0 + az*0));
        o[6]  += q*(2.5*az*ax*ax - 0.5*a2*(az*1 + ax*0 + ax*0));
        o[7]  += q*(2.5*az*ay*ay - 0.5*a2*(az*1 + ay*0 + ay*0));
        o[8]  += q*(2.5*az*az*az - 0.5*a2*(az*1 + az*1 + az*1)); // trace is done
        o[9]  += q*(2.5*ax*ax*ay - 0.5*a2*(ax*0 + ax*0 + ay*1));
        o[10] += q*(2.5*ax*ax*az - 0.5*a2*(ax*0 + ax*0 + az*1));
        o[11] += q*(2.5*ax*ay*az - 0.5*a2*(ax*0 + ay*0 + az*0));
        o[12] += q*(2.5*ay*ax*ay - 0.5*a2*(ay*0 + ax*1 + ay*0));
        o[13] += q*(2.5*ay*ax*az - 0.5*a2*(ay*0 + ax*0 + az*0));
        o[14] += q*(2.5*ay*ay*az - 0.5*a2*(ay*0 + ay*0 + az*1));
        o[15] += q*(2.5*az*ax*ay - 0.5*a2*(az*0 + ax*0 + ay*0));
        o[16] += q*(2.5*az*ax*az - 0.5*a2*(az*0 + ax*1 + az*0));
        o[17] += q*(2.5*az*ay*az - 0.5*a2*(az*0 + ay*1 + az*0));

    }   

    for (n=0; n<3; n++)
        multipole1[n] *= system.constants.eA2D;

    for (n=0; n<3; n++)
        dipole_magnitude += multipole1[n]*multipole1[n];
    dipole_magnitude = sqrt(dipole_magnitude);

    for (n=0; n<6; n++)
        multipole2[n] *= system.constants.eA2D;

    for (n=0; n<18; n++)
        o[n] *= system.constants.eA2D;

    printf("Multipole 0 (total charge), e\n");
    printf("{ %9.5f   } \n\n", multipole0);

    printf("Multipole 1 (dipole moment), Debye\n");
    printf("{ %9.5f %9.5f %9.5f   } \n", multipole1[0], multipole1[1], multipole1[2]);
    printf("(dipole magnitude)  = { %9.5f   } \n\n", dipole_magnitude);
    
    printf("Multipole 2 (quadrupole moment), Debye A\n");
    printf("{ %9.5f %9.5f %9.5f   } \n", multipole2[0], multipole2[3], multipole2[4]);
    printf("{ %9.5f %9.5f %9.5f   } \n", multipole2[3], multipole2[1], multipole2[5]);
    printf("{ %9.5f %9.5f %9.5f   } \n\n", multipole2[4], multipole2[5], multipole2[2]);

    printf("Multipole 3 (octapole moment), Debye A^2\n");
    printf("{ %9.5f %9.5f %9.5f\n", o[0], o[9], o[10]);
    printf("{ %9.5f %9.5f %9.5f\n", o[9], o[1], o[11]);
    printf("{ %9.5f %9.5f %9.5f\n", o[10], o[11], o[2]);
    printf("                                  %9.5f %9.5f %9.5f\n", o[3], o[12], o[13]);
    printf("                                  %9.5f %9.5f %9.5f\n", o[12], o[4], o[14]);
    printf("                                  %9.5f %9.5f %9.5f\n", o[13], o[14], o[5]);
    printf("                                                                  %9.5f %9.5f %9.5f }\n", o[6], o[15], o[16]);
    printf("                                                                  %9.5f %9.5f %9.5f }\n", o[15], o[7], o[17]);
    printf("                                                                  %9.5f %9.5f %9.5f }\n\n", o[16], o[17], o[8]);


    double r, potential=0, sig, eps, sr6;
    double lj=0,es=0;
    
    // Repulsion-Dispersion from long-range (non-bonded)
    // LJ RD
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        for (j=i+1; j<system.molecules[0].atoms.size(); j++) {
            sig = lj_lb_sig(system.molecules[0].atoms[i].sig, system.molecules[0].atoms[j].sig);
            eps = lj_lb_eps(system.molecules[0].atoms[i].eps, system.molecules[0].atoms[j].eps);
            if (sig == 0 || eps == 0) continue;
            double* distances = getDistanceXYZ(system, 0,i,0,j);
            r = distances[3];
            if (r >= system.pbc.mincutoff) {
                sr6 = sig/r;
                sr6 *= sr6;
                sr6 *= sr6*sr6;
                lj += 4.0*eps*(sr6*sr6 - sr6);
            } // end if r >= rmin
        } // end j
    } // end i

    // Coulombic interactions
    double chargeprod;
    // ES
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        for (j=i+1; j<system.molecules[0].atoms.size(); j++) {
            chargeprod = system.molecules[0].atoms[i].C * system.molecules[0].atoms[j].C;

            if (!chargeprod) continue;

            double* distances = getDistanceXYZ(system,0,i,0,j);
            r = distances[3];
            if (1) { //system.pbc.mincutoff) {
                es += chargeprod/r;
            }
        }
    } 

    // Nuclear repulsions
    double nucrepul=0;
    double Zi, Zj; // atomic num's
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        for (j=i+1; j<system.molecules[0].atoms.size(); j++) {
            Zi = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[i].name];
            Zj = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[j].name];
            
            double* distances = getDistanceXYZ(system,0,i,0,j);
            r = distances[3];

            nucrepul += (Zi)*(Zj)/r;
        }
    }

    // Electron-nucleus attractions
    // THIS ASSUMES NEUTRALITY OF THE MOLECULE
    double elec=0, negi, posj;
    //double bohr = system.constants.bohr;
    //double reducefactor;
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        for (j=0; j<system.molecules[0].atoms.size(); j++) {
            Zi = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[i].name];
            Zj = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[j].name];
            negi = -Zi;
            posj = Zj;

            double* distances = getDistanceXYZ(system,0,i,0,j);
            r = distances[3];

            if (i==j) {
                // 1 - (x-1)/x is a decreasing factor from 1 to 0, for Hydrogen, this factor is 1,
                // for large elements, it is closer to zero.
                // we use this to reduce the effective average distance between electrons within atoms.
                double Ni = (double)system.constants.elements[system.molecules[0].atoms[i].name];
                //reducefactor = 1.0 - (Ni - 1.0)/Ni;
                double radius = system.constants.radii[system.molecules[0].atoms[i].name];
                r = avgr(radius);//*reducefactor;
                elec += (negi*posj/r); //  approx r as fraction of bohr radius
            } else elec += negi*posj/r;


        }
    }

    double elecrepul=0; double Ni; double Ninteracts;
    // each atom's electrons repel their local neighbors (on the same atom)
    // THIS ASSUMES NEUTRALITY OF THE MOLECULE
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        Ni = (double)system.constants.elements[system.molecules[0].atoms[i].name];
        Ninteracts = (Ni*Ni - Ni)/2.;
        elecrepul += system.constants.E2REDUCED*-1*system.constants.E2REDUCED*-1*Ninteracts;
    }
    // now contributions from elec's on different atoms
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        for (j=i+1; j<system.molecules[0].atoms.size(); j++) {
            Zi = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[i].name];
            Zj = system.constants.E2REDUCED*system.constants.elements[system.molecules[0].atoms[j].name];
            double* distances = getDistanceXYZ(system, 0,i,0,j);
            r = distances[3];
            elecrepul += Zi*Zj/r;
        
        }
    }    


    double toh = system.constants.K2Eh;
    potential = lj + es + nucrepul + elec + elecrepul;
   
    printf("========================================================================================\n"); 
    printf("The below is a very crude, pseudo-quantum mechanical approximation to the molecular energy:\n\n");
    
    printf("LJ RD (>= %f)  = %f Eh;\nES                   = %f Eh\nNuc. Repul.          = %f Eh\nElectronic Disp.     = %f Eh\nElectronic Repul.    = %f Eh\n\n", system.pbc.mincutoff, lj*toh, es*toh, nucrepul*toh, elec*toh, elecrepul*toh);
    
    printf("Total Energy         = %9.5f Eh\n                     = %9.5f K\n                     = %9.5f kJ/mol\n", potential*system.constants.K2Eh, potential, potential*system.constants.K2KJMOL);

    return;
}
