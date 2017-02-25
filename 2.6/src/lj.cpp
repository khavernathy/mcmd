#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

double self_lj_lrc(System &system) {
    double potential=0;
    double cutoff = system.constants.cutoff;
    double volume = system.constants.volume;

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
           if (system.molecules[i].MF != "F") {
            double sig = system.molecules[i].atoms[j].sig;
            double eps = system.molecules[i].atoms[j].eps;
    
            if (!(sig == 0 || eps == 0)) {
            double sig3 = sig*sig*sig;
            double sigcut = fabs(sig)/cutoff;
            double sigcut3 = sigcut*sigcut*sigcut;
            double sigcut9 = sigcut3 * sigcut3 * sigcut3;

//            printf("sig= %f, eps= %f, sr = %f, s3 = %f, sr3 = %f, sr9 = %f, V = %f\n", sig, eps, sigcut, sig3, sigcut3, sigcut9, system.constants.volume);

            double this_self_lrc = (16.0/3.0)*M_PI*eps*sig3*(sigcut9 - sigcut3)/volume;
            potential += this_self_lrc;
            } // if nonzero sig/eps
        }
    } // end for j atom
    } // end for i molecule
    return potential;
}

double lj(System &system) {

    double total_pot=0, total_lj=0, total_rd_lrc=0, total_rd_self_lrc = 0;
    double cutoff = system.constants.cutoff;
    double volume = system.constants.volume;

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
        double r,sr,sr2,sr6;
        r = getDistance(system, i, j, k, l);
    
        sr = sig/r;
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
         
            double sig3 = sig*sig*sig;
            double sigcut = fabs(sig)/cutoff;
            double sigcut3 = sigcut * sigcut * sigcut;
            double sigcut9 = sigcut3 * sigcut3 * sigcut3;

            double this_rd_lrc = (16.0/3.0)*M_PI*eps*sig3*(sigcut9 - sigcut3)/volume;
            total_rd_lrc += this_rd_lrc;
            total_pot += this_rd_lrc;
        } // end RD LRC
        } // if nonzero sig/eps

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

    system.constants.lj_lrc = total_rd_lrc;
    system.constants.lj_self_lrc = total_rd_self_lrc;
    system.constants.lj = total_lj;
    system.constants.total_rd = total_pot;

//printf("potential: %f, self = %f, lrc = %f, lj = %f\n", total_pot, total_rd_self_lrc, total_rd_lrc, total_lj);
    return total_pot;

}
