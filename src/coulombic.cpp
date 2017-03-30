#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

/* THE THREE FUNCTIONS coulombic_self, coulombic_real,
coulombic_reciprocal are using EWALD method for computation
of electrostatic energy. 

I borrowed these functions, more or less, from Belof & McLaughlin et. al.

*/


/* entire system self potential sum */
// only changes when N changes.
double coulombic_self(System &system) {
    
    double potential=0.0, charge;
    double alpha=system.constants.ewald_alpha;
    double sqrtPI = sqrt(M_PI);
 
    // loop all atoms but skip frozen atoms   
    if (system.stats.MCstep == 0 || system.constants.ensemble == "uvt") { // only changes if N changes
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
    return potential; // i dont think it should be multiplied by ke
}

/* coloumbic_real Ewald result */
double coulombic_real(System &system) {
    
    double potential=0.0, pair_potential=0.0;
    double alpha=system.constants.ewald_alpha;
    double erfc_term; // = erfc(alpha*r);
    double r;    

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i; k < system.molecules.size(); k++) {
    for (int l = 0; l < system.molecules[k].atoms.size(); l++) {    
        if (system.molecules[i].frozen && system.molecules[k].frozen) continue; // skip frozens
        if (system.molecules[i].atoms[j].C == 0 || system.molecules[k].atoms[l].C == 0) continue; // skip 0-energy
       
        pair_potential = 0; 
                
        // calculate distance between atoms
        double* distances = getDistanceXYZ(system,i,j,k,l);
        r = distances[3];
        //r = system.pairs[i][j][k][l].r;

        if (r < system.pbc.cutoff && (i < k)) { // only pairs and not beyond cutoff
            erfc_term = erfc(alpha*r);
            pair_potential += system.molecules[i].atoms[j].C * system.molecules[k].atoms[l].C * erfc_term / r;  // positive (inter)
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
    return potential; // I think it's double counted. Compared to MPMC this is correct. 
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
            system.molecules[i].atoms[j].V += charge1 * charge2 * erfc_term / r;
        } else if (i == k && j != l) { // self molecule interaction
            for (int n=0; n<3; n++) {
                holder = charge1*charge2*erf(alpha*r) / rsq * u[n];
                system.molecules[i].atoms[j].force[n] += holder;
                system.molecules[k].atoms[l].force[n] -= holder;
                
            }

            system.molecules[i].atoms[j].V -= (charge1 * charge2 * erf(alpha*r) / r); // negative (intra)
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
   if (system.constants.ensemble == "npt") {
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



