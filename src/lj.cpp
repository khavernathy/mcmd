#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

double self_lj_lrc(System &system) {
    double potential=0;
    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
           if (system.molecules[i].MF != "F") {
            double sig = system.molecules[i].atoms[j].sig;
            double eps = system.molecules[i].atoms[j].eps;
    
            if (!(sig == 0 || eps == 0)) {
            double sig3 = fabs(sig);
            sig3 *= sig3*sig3;
            double sigcut = fabs(sig)/cutoff;
            double sigcut3 = sigcut*sigcut*sigcut;
            double sigcut9 = sigcut3 * sigcut3 * sigcut3;


            double this_self_lrc = (16.0/3.0)*M_PI*eps*sig3*((1.0/3.0)*sigcut9 - sigcut3)/volume;
            potential += this_self_lrc;
            } // if nonzero sig/eps
        }
    } // end for j atom
    } // end for i molecule
    return potential;
}

double lj(System &system) {

    double total_pot=0, total_lj=0, total_rd_lrc=0, total_rd_self_lrc = 0;
    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i+1; k < system.molecules.size(); k++) {
    for (int l =0; l < system.molecules[k].atoms.size(); l++) {

        // do mixing rules
        double eps = system.molecules[i].atoms[j].eps,sig=system.molecules[i].atoms[j].sig;
        if (eps != system.molecules[k].atoms[l].eps)
            eps = sqrt(system.molecules[i].atoms[j].eps * system.molecules[k].atoms[l].eps);
       
        if (sig != system.molecules[k].atoms[l].sig)
         sig = 0.5 * (system.molecules[i].atoms[j].sig + system.molecules[k].atoms[l].sig);

        if (sig == 0 || eps == 0) continue; // skip 0 energy interactions

        // calculate distance between atoms
        double r,sr,sr2,sr6;
        double* distances = getDistanceXYZ(system, i, j, k, l);
        r = distances[3];

    
        sr = sig/r; //printf("r=%f\n",r);
        sr2 = sr*sr;
        sr6 = sr2*sr2*sr2; //;

        // ============================ LJ potential =============================

        // 1) Normal LJ: only apply if long range corrections are off, or if on and r<cutoff
        if ((system.constants.rd_lrc == "off" || r <= cutoff)) {
            double this_lj = 4.0*eps*(sr6*sr6 - sr6);
            total_lj += this_lj;    //;
            total_pot += this_lj;
        }
        // 2) Long range corr.: apply RD long range correction if needed
        // http://www.seas.upenn.edu/~amyers/MolPhys.pdf
        if (system.constants.rd_lrc == "on") {
         
            double sig3 = fabs(sig);
            sig3 *= sig3*sig3;
            double sigcut = fabs(sig)/cutoff;
            double sigcut3 = sigcut * sigcut * sigcut;
            double sigcut9 = sigcut3 * sigcut3 * sigcut3;

            double this_rd_lrc = (16.0/3.0)*M_PI*eps*sig3*((1.0/3.0)*sigcut9 - sigcut3)/volume;
            total_rd_lrc += this_rd_lrc;
            total_pot += this_rd_lrc;
        } // end RD LRC

    }  // loop l
    } // loop k 
    } //loop j
    } // loop i
    // DONE WITH PAIR INTERACTIONS

    // 3) LJ LRC self energy
    // only do for individual non-frozen atoms
    if (system.constants.rd_lrc == "on") {  
        total_rd_self_lrc = self_lj_lrc(system);
        total_pot += total_rd_self_lrc;
    } // end LRC self contribution.

    system.stats.lj_lrc.value = total_rd_lrc; 
    system.stats.lj_self_lrc.value = total_rd_self_lrc; 
    system.stats.lj.value = total_lj; 

//printf("potential: %f, self = %f, lrc = %f, lj = %f\n", total_pot, total_rd_self_lrc, total_rd_lrc, total_lj);
    return total_pot;

}


double lj_force(System &system) {

    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    double d[3], sr, eps, sig, sr2, sr6, r,rsq,r6,s2,s6, f[3];

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i+1; k < system.molecules.size(); k++) {
    for (int l =0; l < system.molecules[k].atoms.size(); l++) {

        // do mixing rules
        double eps,sig;
        eps = sqrt(system.molecules[i].atoms[j].eps * system.molecules[k].atoms[l].eps);
        sig = 0.5 * (system.molecules[i].atoms[j].sig + system.molecules[k].atoms[l].sig);

        if (!(sig == 0 || eps == 0)) {
        // calculate distance between atoms
        double* distances = getDistanceXYZ(system, i, j, k, l);
        r = distances[3];    
        rsq=r*r;
        for (int n=0; n<3; n++) d[n] = distances[n];

        r6 = rsq*rsq*rsq;
        s2 = sig*sig;
        s6 = s2*s2*s2;
    
                if (i != k) { // don't do self-interaction for potential.
                    sr = sig/r;
                    sr2 = sr*sr;    
                    sr6 = sr2*sr2*sr2;
                }

        if ((system.constants.rd_lrc == "off" || r <= cutoff)) {
            for (int n=0; n<3; n++) {
                f[n] = 24.0*d[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                system.molecules[i].atoms[j].force[n] += f[n];
                system.molecules[k].atoms[l].force[n] -= f[n];
            }
        }

        system.molecules[i].atoms[j].V += 4.0*eps*(sr6*sr6 - sr6);
        } // if nonzero sig/eps

    }  // loop l
    } // loop k 
    } //loop j
    } // loop i
    // DONE WITH PAIR INTERACTIONS
}
