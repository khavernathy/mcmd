#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

#define HBAR2 1.11211999e-68
#define HBAR4 1.23681087e-136
#define KB2 1.90619525e-46
#define KB 1.3806503e-23

// get the Feynman-Hibbs correction for a pair of atoms
double lj_fh_corr(System &system, int i,int k, double r, double term12, double term6, double sig, double eps) {
    double reduced_mass;
    double dE, d2E, d3E, d4E; //energy derivatives
    double corr;
    double ir = 1.0/r;
    double ir2 = ir*ir;
    double ir3 = ir2*ir;
    double ir4 = ir3*ir;
    int order = system.constants.fh_order;

    if (order != 2 && order != 4) return NAN;

    reduced_mass = system.molecules[i].mass*system.molecules[k].mass / (system.molecules[i].mass + system.molecules[k].mass);
    dE = -24.0*eps*(2.0*term12 - term6)*ir;
    d2E = 24.0*eps*(26.0*term12 - 7.0*term6)*ir2;

    // 2nd order corr
    corr = 1.0e20 *
        (HBAR2/(24.0*KB*system.constants.temp * reduced_mass)) *
        (d2E + 2.0*dE/r);

    if (order == 4) {
        d3E = -1344.0*eps*(6.0*term12 - term6) * ir3;
        d4E = 12096.0*eps*(10.0*term12 - term6) * ir4;

        // 4th order corr
        corr += 1.0e40 *
            (HBAR4/(1152.0*KB2*system.constants.temp*system.constants.temp*reduced_mass*reduced_mass)) *
            ( 15.0*dE*ir3 + 4.0*d3E*ir + d4E);
    }

    return corr;
}

double self_lj_lrc(System &system) {
    double potential=0;
    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    double sig, eps, sig3,sigcut,sigcut3,sigcut9;
    double this_self_lrc;

    if (system.stats.MCstep == 0 || system.constants.ensemble == ENSEMBLE_UVT || system.constants.ensemble == ENSEMBLE_NPT) { // only changes if N or V changes
    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
           if (!system.molecules[i].frozen) {
            sig = system.molecules[i].atoms[j].sig;
            eps = system.molecules[i].atoms[j].eps;
    
            if (!(sig == 0 || eps == 0)) {
            sig3 = fabs(sig);
            sig3 *= sig3*sig3;
            sigcut = fabs(sig)/cutoff;
            sigcut3 = sigcut*sigcut*sigcut;
            sigcut9 = sigcut3 * sigcut3 * sigcut3;


            this_self_lrc = (16.0/3.0)*M_PI*eps*sig3*((1.0/3.0)*sigcut9 - sigcut3)/volume;
            potential += this_self_lrc;
            } // if nonzero sig/eps
        }
    } // end for j atom
    } // end for i molecule
    } // end if self LRC changes
    else {
        potential = system.stats.lj_self_lrc.value;
    }
    return potential;
}



double lj(System &system) {
    double total_pot=0, total_lj=0, total_rd_lrc=0, total_rd_self_lrc = 0;
    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    int i,j,k,l,index;
    double this_lj;
    double r,sr6;
    
    for (i = 0; i < system.molecules.size(); i++) {
    for (j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (k = i+1; k < system.molecules.size(); k++) {
    for (l =0; l < system.molecules[k].atoms.size(); l++) {

        //if (system.pairs[i][j][k][l].recalculate) {

        // do mixing rules
        double eps = system.molecules[i].atoms[j].eps,sig=system.molecules[i].atoms[j].sig;
        if (eps != system.molecules[k].atoms[l].eps)
            eps = sqrt(eps * system.molecules[k].atoms[l].eps);
       
        if (sig != system.molecules[k].atoms[l].sig)
         sig = 0.5 * (sig + system.molecules[k].atoms[l].sig);

        if (sig == 0 || eps == 0) continue; // skip 0 energy interactions

        // calculate distance between atoms
        double* distances = getDistanceXYZ(system, i, j, k, l);
        r = distances[3];
        //r = system.pairs[i][j][k][l].r;    

        sr6 = sig/r; //printf("r=%f\n",r);
        sr6 *= sr6;
        sr6 *= sr6*sr6; //;

        // ============================ LJ potential =============================

        // 1) Normal LJ: only apply if long range corrections are off, or if on and r<cutoff
        if ((!system.constants.rd_lrc || r <= cutoff)) {
            this_lj = 4.0*eps*(sr6*sr6 - sr6);
            total_lj += this_lj;    //;
            total_pot += this_lj;

            if (system.constants.feynman_hibbs)
                total_pot += lj_fh_corr(system, i,k, r, sr6*sr6, sr6, sig, eps);
        }
        
        /*
        // end if recalculate
        } else {
            if ((system.constants.rd_lrc=="off" || r <= cutoff)) {
                this_lj = system.pairs[i][j][k][l].lj;
                total_lj += this_lj;
                total_pot += this_lj;
            }
        }
        */
    }  // loop l
    } // loop k 
    } //loop j
    } // loop i


    // 2) Long range corr.: apply RD long range correction if needed
        // http://www.seas.upenn.edu/~amyers/MolPhys.pdf
    if (system.constants.rd_lrc) {
        if (system.stats.MCstep == 0 || system.constants.ensemble == ENSEMBLE_NPT || system.constants.ensemble == ENSEMBLE_UVT) { // lrc only changes if volume or N changes.
        for (i=0; i < system.molecules.size(); i++) {
        for (j=0; j< system.molecules[i].atoms.size(); j++) {
        for (k=0; k <system.molecules.size(); k++) {
        for (l=0; l <system.molecules[k].atoms.size(); l++) {

        if (system.molecules[i].frozen && system.molecules[k].frozen) continue; // skip frozens
        if (i<k || (i==k && j<l)) {

        // do mixing rules
        double eps = system.molecules[i].atoms[j].eps,sig=system.molecules[i].atoms[j].sig;
        if (eps != system.molecules[k].atoms[l].eps)
            eps = sqrt(eps * system.molecules[k].atoms[l].eps);
        if (sig != system.molecules[k].atoms[l].sig)
         sig = 0.5 * (sig + system.molecules[k].atoms[l].sig);
        if (sig == 0 || eps == 0) continue; // skip 0 energy interactions         

            double sig3 = fabs(sig);
            sig3 *= sig3*sig3;
            double sigcut = fabs(sig)/cutoff;
            double sigcut3 = sigcut * sigcut * sigcut;
            double sigcut9 = sigcut3 * sigcut3 * sigcut3;

            double this_rd_lrc = (16.0/3.0)*M_PI*eps*sig3*((1.0/3.0)*sigcut9 - sigcut3)/volume;
            total_rd_lrc += this_rd_lrc;
            total_pot += this_rd_lrc;
        } // pair condition
        }
        }
        }
        } // end 4 atom loops.
        } // end if recalculate lrc
        else {
           total_rd_lrc = system.stats.lj_lrc.value;
           total_pot += total_rd_lrc;
        }
    } // end if RD LRC is on
    // DONE WITH PAIR INTERACTIONS

    // 3) LJ LRC self energy
    // only do for individual non-frozen atoms
    if (system.constants.rd_lrc) {  
        total_rd_self_lrc = self_lj_lrc(system);
        total_pot += total_rd_self_lrc;
    } // end LRC self contribution.

    system.stats.lj_lrc.value = total_rd_lrc; 
    system.stats.lj_self_lrc.value = total_rd_self_lrc; 
    system.stats.lj.value = total_lj; 

//printf("potential: %f, self = %f, lrc = %f, lj = %f\n", total_pot, total_rd_self_lrc, total_rd_lrc, total_lj);
    return total_pot;

}


void lj_force(System &system) {

    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    double d[3], sr, eps, sig, sr2, sr6, r,rsq,r6,s2,s6, f[3];
    int count=0; // for the pair values

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

        if ((!system.constants.rd_lrc || r <= cutoff)) {
            for (int n=0; n<3; n++) {
                f[n] = 24.0*d[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                system.molecules[i].atoms[j].force[n] += f[n];
                system.molecules[k].atoms[l].force[n] -= f[n];
            }

/*
            // this dot product is calc'd to get pressure.
            if (system.molecules[i].MF == "M" && system.molecules[k].MF == "M") { // make sure it's gas-to-gas
            system.pairs[count].fdotr = (f[0] * distances[0] +
            f[1] * distances[1] +
            f[2] * distances[2]); // F.r for the sum for pressure calculator.
            count++;       
            }
*/
        }

        //system.molecules[i].atoms[j].V += 4.0*eps*(sr6*sr6 - sr6);
        } // if nonzero sig/eps
    }  // loop l
    } // loop k 
    } //loop j
    } // loop i
    // DONE WITH PAIR INTERACTIONS
}

void lj_force_nopbc(System &system) {

    double cutoff = system.pbc.cutoff;
    double volume = system.pbc.volume;
    double d[3], sr, eps, sig, sr2, sr6, r,rsq,r6,s2,s6, f[3];
    int count=0; // for the pair values

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

            for (int n=0; n<3; n++) {
                f[n] = 24.0*d[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                system.molecules[i].atoms[j].force[n] += f[n];
                system.molecules[k].atoms[l].force[n] -= f[n];
            }

            system.molecules[i].atoms[j].V += 4.0*eps*(sr6*sr6 - sr6);
        } // if nonzero sig/eps
    }  // loop l
    } // loop k 
    } //loop j
    } // loop i
    // DONE WITH PAIR INTERACTIONS
}
