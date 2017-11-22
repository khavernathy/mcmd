#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

using namespace std;

/* POTENTIAL ENERGY STUFF */
// Credit to Adam Hogan for originally coding this
// algorithm in MPMC.

// Long range energy for a pair
double tt_lrc(System &system, double c6, double c8, double c10) {
    const double rc = system.pbc.cutoff;
    const double rc2 = rc*rc;
    const double rc5 = rc2*rc2*rc;
    return -4.0*M_PI*(c6/(3.0*rc*rc2) + c8/(5.0*rc5) + c10/(7.0*rc5*rc2))/system.pbc.volume;
}

// self energy for an atom
double tt_self(System &system, int i, int j) {
    const double rc = system.pbc.cutoff;
    const double rc2 = rc*rc;
    const double rc5 = rc2*rc2*rc;
    const double c6 = system.molecules[i].atoms[j].c6;
    const double c8 = system.molecules[i].atoms[j].c8;
    const double c10 = system.molecules[i].atoms[j].c10;
    return -4.0*M_PI*(c6/(3.0*rc2*rc) + c8/(5.0*rc5) + c10/(7.0*rc5*rc2))/system.pbc.volume;
}

// f_2n(bR) damping function in the paper.
double tt_damp(int n, double br) {
    double sum=0;
    int i;
    for (i=0;i<=n;i++)
        sum += pow(br,i)/factorial(i);

    const double result = 1.0 - exp(-br)*sum;

    if (result>0.000000001)
        return result;
    else
        return 0.0;
}

/* the Tang-Toennies potential for the entire system */
double tt(System &system) {
    double potential=0, repulsive=0, attractive = 0;
    double c6,c8,c10,sig,b; // for the pair itself.

    // energy from atom pairs
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            for (int k=i; k<system.molecules.size(); k++) {
                for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                    // skip frozens
                    if (system.molecules[i].frozen && system.molecules[k].frozen) continue;
                    
                    // do mixing
                    c6  = tt_c6(system.molecules[i].atoms[j].c6, system.molecules[k].atoms[l].c6);
                    c8  = tt_c8(system.molecules[i].atoms[j].c8, system.molecules[k].atoms[l].c8);
                    c10 = tt_c10(c6,c8);
                    sig = tt_sigma(system.molecules[i].atoms[j].sig, system.molecules[k].atoms[l].sig);
                    b   = tt_b(system.molecules[i].atoms[j].eps, system.molecules[k].atoms[l].eps);  

                    // LRC applies to all unique pairs (and should include intramolec.)
                    if (system.constants.rd_lrc && ((i==k && j<l) || (i<k))   ) {
                        potential += tt_lrc(system, c6, c8, c10); 
                    }

                    double* distances = getDistanceXYZ(system, i,j,k,l);
                    const double r = distances[3];
                    const double r2 = r*r;
                    const double r4 = r2*r2;
                    const double r6 = r4*r2;
                    const double r8 = r6*r2;
                    const double r10 = r8*r2;
                    //if (r <= system.pbc.cutoff && i<k) { // not intramolecular
                    if (i<k) { // not intramolecular 
                        repulsive = 315.7750382111558307123944638*exp(-b*(r-sig));

                        attractive = -tt_damp(6,b*r)*c6/r6 - tt_damp(8,b*r)*c8/r8 - tt_damp(10,b*r)*c10/r10;

                        // total potential from this pair with LRC.
                        potential += attractive + repulsive;

                    }// end if within r_c

                }// end l
            }// end k
        }// end j
    }// end i

    // and calculate self LRC contribution to energy if needed.
    if (system.constants.rd_lrc) {
        for (int i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue; // skip frozen
            for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                potential += tt_self(system, i, j);
            }
        }
    }
    
    return potential;
}



/* FORCE STUFF, -dU/dr */
// Douglas Franz, 2017

void tt_forces(System &system) {
    double c6,c8,c10,sig,b; // for the pair itself.
    double holder; // temp to store forces
    double sum,innersum,A,B; // temp to store damping sum for dispersion forces
    double ctmp; // C parameter tmp

    // energy from atom pairs
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            for (int k=i; k<system.molecules.size(); k++) {
                for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                    // skip frozens
                    if (system.molecules[i].frozen && system.molecules[k].frozen) continue;
                    
                    // do mixing
                    c6  = tt_c6(system.molecules[i].atoms[j].c6, system.molecules[k].atoms[l].c6);
                    c8  = tt_c8(system.molecules[i].atoms[j].c8, system.molecules[k].atoms[l].c8);
                    c10 = tt_c10(c6,c8);
                    sig = tt_sigma(system.molecules[i].atoms[j].sig, system.molecules[k].atoms[l].sig);
                    b   = tt_b(system.molecules[i].atoms[j].eps, system.molecules[k].atoms[l].eps);  

                    double* distances = getDistanceXYZ(system, i,j,k,l);
                    const double r = distances[3];
                    double ir_array[13]; // from r^0 to r^-12
                    ir_array[0] = 1.0;
                    ir_array[1] = 1./distances[3];
                    for (int z=2;z<13;z++) ir_array[z] = ir_array[z-1]*ir_array[1]; // assign r^-2, r^-3, ... r^-12
                    
                    //if (r <= system.pbc.cutoff && i<k) { // not intramolecular
                    if (i<k) { // not intramolecular 
                        for (int q=0; q<3; q++) { // 3 DIMS
                            // repulsion force
                            holder += 315.7750382111558307123944638*b*distances[q]*exp(b*sig)*exp(-b*r)/r;

                            // damped dispersion force (attractive)
                            sum = 0;
                            for (int d=6; d<=10; d+=2) { // so we cover d6, d8, d10 terms.
                                if (d==6) ctmp = c6;
                                else if (d==8) ctmp = c8;
                                else if (d==10) ctmp = c10;

                                    sum += -d*distances[q]*ctmp*ir_array[d+2]; // e.g. for d10, we divide by ir_array[12], which is r^-12.
                                    innersum = 0;
                                    for (int n=0; n<d+2; n++) {  // e.g., for 10-damping, this sum goes from n=0 -> n=11
                                        if (n==0) {
                                            A = (d-n)*pow(b,n)/factorial(n);
                                            B = 0;
                                        } else if (n==d+1) {
                                            A = 0;
                                            B = pow(b,n)/factorial(n-1);
                                        } else {
                                            A = (d-n)*pow(b,n)/factorial(n);
                                            B = pow(b,n)/factorial(n-1);
                                        }
                                        innersum += (A+B)*ir_array[d+1-n]; // so for n=0 and 10-damping, the last term in sum, we take ir_array[11] = r^-11
                                    }
                                    sum += innersum*ctmp*distances[q]*exp(-b*r)/r;
                                    
                            } // end damping functions
                            holder += sum;    
                        

                            // Apply Newton pair for this dimension
                            system.molecules[i].atoms[j].force[q] += holder;
                            system.molecules[k].atoms[l].force[q] -= holder; 
                        } // end 3D
                        
                    }// end if not intramolec.
                }// end l
            }// end k
        }// end j
    }// end i

    return;
}




