#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

/* THE THREE FUNCTIONS coulombic_self, coulombic_real,
coulombic_reciprocal are using EWALD method for computation
of electrostatic energy. 
*/

#define SQRTPI 1.77245385091
#define HBAR2 1.11211999e-68
#define HBAR4 1.23681087e-136
#define KB2 1.90619525e-46
#define KB 1.3806503e-23

double es_fh_corr(System &system, int i, int k, double r, double gaussian_term, double erfc_term) {
    double dE, d2E, d3E, d4E; 
    double corr;
    double rr = r*r;
    double ir = 1.0/r;
    double ir2 = ir*ir;
    double ir3 = ir*ir2;
    double ir4 = ir2*ir2;
    double order = system.constants.fh_order;
    double alpha = system.constants.ewald_alpha;
    double a2 = alpha*alpha;
    double a3 = a2*alpha;
    double a4 = a3*alpha;
    double reduced_mass = (system.molecules[i].mass * system.molecules[k].mass)/(system.molecules[i].mass + system.molecules[k].mass);

    if (order != 2 && order != 4) return NAN;

    dE = -2.0*alpha*gaussian_term/(r*SQRTPI) - erfc_term*ir2;
    d2E = (4.0/SQRTPI)*gaussian_term*(a3 + 1.0*ir2) + 2.0*erfc_term*ir3;

    corr = 1.0e20 * (HBAR2/(24.0*KB*system.constants.temp*reduced_mass)) * (d2E + 2.0*dE/r);

    if (order == 4) {
        d3E = (gaussian_term/SQRTPI) * (-8.0*(a3*a2)*r - 8.0*(a3)/r - 12.0*alpha*ir3) 
            - 6.0*erfc(alpha*r)*ir4;
        d4E = (gaussian_term/SQRTPI) * (-8.0*a3*a2 + 16.0*a3*a4*rr + 32.0*a3
            *ir2 + 48.0*ir4 ) + 24.0*erfc_term*(ir4*ir);
    
        corr += 1.0e40*(HBAR4/(1152.0*(KB2*system.constants.temp*system.constants.temp * reduced_mass*reduced_mass))) * (15.0*dE*ir3 + 4.0*d3E/r + d4E);
    }

    return corr;

}


/* entire system self potential sum */
// only changes when N changes.
double coulombic_self(System &system) {
    
    double potential=0.0, charge;
    double alpha=system.constants.ewald_alpha;
    double sqrtPI = sqrt(M_PI);
 
    // loop all atoms but skip frozen atoms   
    if (system.stats.MCstep == 0 || system.constants.ensemble == ENSEMBLE_UVT) { // only changes if N changes
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {

            if (!system.molecules[i].atoms[j].frozen &&
                system.molecules[i].atoms[j].C != 0) {
                charge = system.molecules[i].atoms[j].C; 
                potential -= alpha* charge * charge / sqrtPI; 
            }        
        } // end for atom i in molecule j
    } // end for molecule j

    } // end if re-calculate
    else {
        potential = system.stats.es_self.value;
    }
    return potential;
}

/* coloumbic_real Ewald result */
double coulombic_real(System &system) {
    
    double potential=0.0, pair_potential=0.0;
    double alpha=system.constants.ewald_alpha;
    double erfc_term; // = erfc(alpha*r);
    double r;  //  int count =0;
    double gaussian_term;

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i; k < system.molecules.size(); k++) {
    for (int l = 0; l < system.molecules[k].atoms.size(); l++) {    
        if (system.molecules[i].frozen && system.molecules[k].frozen) continue; // skip frozens
        if (system.molecules[i].atoms[j].C == 0 || system.molecules[k].atoms[l].C == 0) continue; // skip 0-energy
       
//        count++;
        pair_potential = 0; 
                
        // calculate distance between atoms
        double* distances = getDistanceXYZ(system,i,j,k,l);
        r = distances[3];
        //r = system.pairs[i][j][k][l].r;

        if (r < system.pbc.cutoff && (i < k)) { // only pairs and not beyond cutoff
            erfc_term = erfc(alpha*r);
            pair_potential += system.molecules[i].atoms[j].C * system.molecules[k].atoms[l].C * erfc_term / r;  // positive (inter)
        
            if (system.constants.feynman_hibbs) {
                gaussian_term = exp(-alpha*alpha*r*r);
                pair_potential += es_fh_corr(system, i, k, r, gaussian_term, erfc_term);
            }
        
        } else if (i == k && j < l) { // self molecule interaction
            pair_potential -= (system.molecules[i].atoms[j].C * system.molecules[k].atoms[l].C * erf(alpha*r) / r); // negative (intra)
        }
        if (std::isnan(potential) == 0) { // CHECK FOR NaN
            potential += pair_potential;
        }

    } // end l
    } // end k
    } // end j
    } // end i 
//    printf("alpha = %f; es_real = %f; count = %i\n", alpha, potential, count);
    return potential; 
}

void coulombic_real_force(System &system) {
    
    double alpha=system.constants.ewald_alpha;
    double erfc_term; // = erfc(alpha*r);
    double charge1, charge2, r,rsq;
    double u[3]; double holder;

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = 0; k < system.molecules.size(); k++) {
    for (int l = 0; l < system.molecules[k].atoms.size(); l++) {
    if (!(system.molecules[i].frozen && system.molecules[k].frozen) &&
        !(system.molecules[i].atoms[j].C == 0 || system.molecules[i].atoms[j].C == 0) ) { // don't do frozen-frozen or zero charge

        charge1 = system.molecules[i].atoms[j].C;
        charge2 = system.molecules[k].atoms[l].C;

        // calculate distance between atoms
        double* distances = getDistanceXYZ(system,i,j,k,l);
        r = distances[3];

        rsq = r*r;
        for (int n=0; n<3; n++) u[n] = distances[n]/r;

        if (r < system.pbc.cutoff && (i < k)) { // only pairs and not beyond cutoff
            erfc_term = erfc(alpha*r);
            for (int n=0; n<3; n++) {
                holder = charge1*charge2*erfc_term/rsq * u[n];
                system.molecules[i].atoms[j].force[n] += holder;
                system.molecules[k].atoms[l].force[n] -= holder;

            }
            //system.molecules[i].atoms[j].V += charge1 * charge2 * erfc_term / r;
        } else if (i == k && j != l) { // self molecule interaction
            for (int n=0; n<3; n++) {
                holder = charge1*charge2*erf(alpha*r) / rsq * u[n];
                system.molecules[i].atoms[j].force[n] += holder;
                system.molecules[k].atoms[l].force[n] -= holder;
                
            }

            //system.molecules[i].atoms[j].V -= (charge1 * charge2 * erf(alpha*r) / r); // negative (intra)
        }
    } // end if not frozen
    } // end l
    } // end k
    } // end j
    } // end i 
}

// Coulombic reciprocal electrostatic energy from Ewald //
double coulombic_reciprocal(System &system) {   
    int p, q, kmax, l[3], i, j;
    double alpha, k[3], k_sq, position_product;
    double SF_re=0, SF_im=0;
    double potential = 0.0;

    alpha = system.constants.ewald_alpha;
    kmax = system.constants.ewald_kmax;   

   // get recip (re-calc only needed for NPT)
   if (system.constants.ensemble == ENSEMBLE_NPT) {
        system.pbc.calcVolume();
        system.pbc.calcRecip(); 
   }

    // fourier sum over a hemisphere (skipping certain points to avoid overcounting the face //
    for (l[0] = 0; l[0] <= kmax; l[0]++) {
        for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
            for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {
        
                // skip if norm is out of sphere
                if (l[0]*l[0] + l[1]*l[1] + l[2]*l[2] > kmax*kmax) continue;   

                // get reciprocal lattice vectors
                for (p=0; p<3; p++) {
                    for (q=0, k[p] = 0; q < 3; q++) {
                        k[p] += 2.0*M_PI*system.pbc.reciprocal_basis[p][q] * l[q];
                    }
                }
                k_sq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                // Structure factor. Loop all atoms.
                SF_re =0.0; SF_im = 0;
                for (i=0; i<system.molecules.size(); i++) {
                    for (j=0; j<system.molecules[i].atoms.size(); j++) {
                        if (system.molecules[i].atoms[j].frozen) continue;
                        if (system.molecules[i].atoms[j].C == 0) continue;

                        // inner product of position vector and k vector
                        position_product = (k[0]*system.molecules[i].atoms[j].pos[0] + k[1]*system.molecules[i].atoms[j].pos[1] + k[2]*system.molecules[i].atoms[j].pos[2]);

                        SF_re += system.molecules[i].atoms[j].C * cos(position_product);
                        SF_im += system.molecules[i].atoms[j].C * sin(position_product);

                    } // end for atom j in molecule i
                } // end for molecule i           

                potential += exp(-k_sq/(4.0*alpha*alpha)) / k_sq * (SF_re*SF_re + SF_im*SF_im);

            } // end for l[2], n
        } // end for l[1], m
    } // end for l[0], l

    potential *= 4.0 * M_PI / system.pbc.volume;

    //printf("coulombic_reciprocal: %f K\n",potential);
   return potential;
}



double coulombic_ewald(System &system) {
    // Ewald method for coulombic
    double potential;
    
    
    double self, real, recip;
    self = coulombic_self(system);
    real = coulombic_real(system);
    recip = coulombic_reciprocal(system);

    system.stats.es_self.value = self; 
    system.stats.es_real.value = real; 
    system.stats.es_recip.value = recip; 

    potential = self + real + recip;
    
    return potential;
}

double coulombic(System &system) { // old super basic coulombic
   // plain old coloumb
   double potential = 0;
   double r, q1, q2;
   
    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i+1; k < system.molecules.size(); k++) {
    for (int l =0; l < system.molecules[k].atoms.size(); l++) { 
        if (!(system.molecules[i].atoms[j].C == 0 || system.molecules[k].atoms[l].C ==0)) {
        
        double* distances = getDistanceXYZ(system,i,j,k,l);
        r = distances[3];
        //r = system.pairs[i][j][k][l].r;

        q1 = system.molecules[i].atoms[j].C;
        q2 = system.molecules[k].atoms[l].C;
        
        potential += q1*q2/r;
        } // end if nonzero charges
    } // end l
    } // end k
    } // end j
    } // end i
    return potential;
}



