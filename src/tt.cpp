#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

using namespace std;

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
    double A; // leading coefficient for repulsive energy
    double lrc; // long-range contribution

    // energy from atom pairs
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            // skip 0-energy
            if (system.molecules[i].atoms[j].eps == 0 || system.molecules[i].atoms[j].sig == 0) continue;
            for (int k=i+1; k<system.molecules.size(); k++) {
                for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                    // skip 0-energy
                    if (system.molecules[k].atoms[l].eps == 0 || system.molecules[k].atoms[l].sig ==0) continue;
                    // skip frozens
                    if (system.molecules[i].frozen && system.molecules[k].frozen) continue;
                    
                    double* distances = getDistanceXYZ(system, i,j,k,l);
                    const double r = distances[3];
                    const double r2 = r*r;
                    const double r4 = r2*r2;
                    const double r6 = r4*r2;
                    const double r8 = r6*r2;
                    const double r10 = r8*r2;
                    if (r <= system.pbc.cutoff) {
                        printf("r = %f\n", r);
                        // do mixing
                        c6  = tt_c6(system.molecules[i].atoms[j].c6, system.molecules[k].atoms[l].c6);
                        c8  = tt_c8(system.molecules[i].atoms[j].c8, system.molecules[k].atoms[l].c8);
                        c10 = tt_c10(c6,c8);
                        sig = tt_sigma(system.molecules[i].atoms[j].sig, system.molecules[k].atoms[l].sig);
                        b   = tt_b(system.molecules[i].atoms[j].eps, system.molecules[k].atoms[l].eps);  
                
                        repulsive = 315.7750382111558307123944638*exp(-b*(r-sig));

                        attractive = -tt_damp(6,b*r)*c6/r6 - tt_damp(8,b*r)*c8/r8 - tt_damp(10,b*r)*c10/r10;

                        // total potential from this pair with LRC.
                        potential += repulsive + attractive; // + tt_lrc(system, c6, c8, c10);

                    }// end if within r_c

                }// end l
            }// end k
        }// end j
    }// end i
/*
    // and calculate self contribution to energy
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) continue; // skip frozen
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            potential += tt_self(system, i, j);
        }
    }
*/
    return potential;
}


