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
/* this quantity is the same every step, later I can save it and not re-calculate */
double coulombic_self(System &system) {
    //initialize to zero
    system.constants.coulombic_self =0.0;
    
    double potential=0.0, charge=0;
    double alpha=system.constants.ewald_alpha;
 
    //printf("alpha = %f\n",alpha);

    // loop all atoms but skip frozen atoms   
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {

            if (system.molecules[i].atoms[j].MF == "M" &&
                system.molecules[i].atoms[j].C != 0) {
                charge = system.molecules[i].atoms[j].C * system.constants.E2REDUCED; 
                   // system.molecules[i].atoms[j].printAll(); 
                 potential -= alpha* charge * charge / sqrt(M_PI); 
                //printf("self pot on mol %i atom %i = %f\n", i, j, alpha*charge*charge/sqrt(M_PI));   
            }        
        } // end for atom i in molecule j
    } // end for molecule j

   //     printf("coulombic_self: %f K\n", system.constants.ke*potential);

   // system.constants.coulombic_self = potential;
    return potential; // i dont think it should be multiplied by ke
}

/* coloumbic_real Ewald result */
double coulombic_real(System &system) {
    
    double potential=0.0, pair_potential=0.0;
    double alpha=system.constants.ewald_alpha;
    double erfc_term; // = erfc(alpha*r);
    double charge1, charge2, r;    

    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = 0; k < system.molecules.size(); k++) {
    for (int l = 0; l < system.molecules[k].atoms.size(); l++) {    
    if (!(system.molecules[i].MF =="F" && system.molecules[k].MF =="F") &&
        !(system.molecules[i].atoms[j].C == 0 || system.molecules[i].atoms[j].C == 0) ) { // don't do frozen-frozen or zero charge
       
        pair_potential = 0; 
        charge1 = system.molecules[i].atoms[j].C*system.constants.E2REDUCED;
        charge2 = system.molecules[k].atoms[l].C*system.constants.E2REDUCED;
                
        // calculate distance between atoms
        r = getDistance(system,i,j,k,l);

//printf("r = %f, r_c = %f, i = %i, j = %i, k = %i, l = %i, q1 = %f, q2 = %f\n", r, system.constants.cutoff, i, j , k, l, charge1, charge2);

        if (r < system.constants.cutoff && (i != k)) { // only pairs and not beyond cutoff
            erfc_term = erfc(alpha*r);
            pair_potential += charge1 * charge2 * erfc_term / r;  // positive (inter)
  //          printf("nonfrozen and r < cutoff\n");
        } else if (i == k && j != l) { // self molecule interaction
    //        printf("frozen real pair\n");
            pair_potential -= (charge1 * charge2 * erf(alpha*r) / r); // negative (intra)
        }
            //printf("isnan? = %i = %d\n",std::isnan(potential),std::isnan(potential));
        if (std::isnan(potential) == 0) { // CHECK FOR NaN
            potential += pair_potential;
        }
        
      //  printf("pair potential: %f\n=============\n", pair_potential);

    } // end if not frozen
    } // end l
    } // end k
    } // end j
    } // end i 
    return 0.5*potential; // I think it's double counted. Compared to MPMC this is correct. 
}

// Coulombic reciprocal electrostatic energy from Ewald //
double coulombic_reciprocal(System &system) {   
    system.constants.coulombic_reciprocal=0;
    int p, q, kmax, l[3];
    double alpha, k[3], k_sq, position_product;
    double SF_re=0, SF_im=0;
    double potential = 0.0;

    alpha = system.constants.ewald_alpha;
    kmax = system.constants.ewald_kmax;   

    // define cubic reciprocals right quick.
    double inverse_volume = 1.0/system.constants.volume;
    double reciprocal_basis[3][3];   
    reciprocal_basis[0][0] = inverse_volume * system.constants.y_length * system.constants.z_length;
    reciprocal_basis[1][1] = inverse_volume * system.constants.x_length * system.constants.z_length;
    reciprocal_basis[2][2] = inverse_volume * system.constants.x_length * system.constants.y_length;

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            if (i != j) reciprocal_basis[i][j] = 0.0;
        }
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
                        k[p] += 2.0*M_PI*reciprocal_basis[p][q] * l[q];
                    }
                }
                k_sq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                // Structure factor. Loop all atoms.
                SF_re =0.0; SF_im = 0;
                for (int i=0; i<system.molecules.size(); i++) {
                    for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                        if (system.molecules[i].atoms[j].MF == "F") continue;
                        if (system.molecules[i].atoms[j].C == 0) continue;

                        // inner product of position vector and k vector
                        position_product = (k[0]*system.molecules[i].atoms[j].pos[0] + k[1]*system.molecules[i].atoms[j].pos[1] + k[2]*system.molecules[i].atoms[j].pos[2]);

                        SF_re += system.molecules[i].atoms[j].C*system.constants.E2REDUCED * cos(position_product);
                        SF_im += system.molecules[i].atoms[j].C*system.constants.E2REDUCED * sin(position_product);

                    } // end for atom j in molecule i
                } // end for molecule i           

                potential += exp(-k_sq/(4.0*alpha*alpha)) / k_sq * (SF_re*SF_re + SF_im*SF_im);

            } // end for l[2], n
        } // end for l[1], m
    } // end for l[0], l

    potential *= 4.0 * M_PI / system.constants.volume;

    //printf("coulombic_reciprocal: %f K\n",potential);
 //   system.constants.coulombic_reciprocal = potential; 
   return potential;
}


double coulombic_ewald(System &system) {
    // Ewald method for coulombic
    double potential;
    
    system.constants.coulombic_self =0;
    system.constants.coulombic_real = 0;
    system.constants.coulombic_reciprocal = 0;
    
    double self, real, recip;
    self = coulombic_self(system);
    real = coulombic_real(system);
    recip = coulombic_reciprocal(system);

    system.constants.coulombic_self = self;
    system.constants.coulombic_real = real;
    system.constants.coulombic_reciprocal = recip;

    potential = self + real + recip;
    
    return potential;
}

double coulombic(System &system) {
   // plain old coloumb
   double potential = 0;
   double r, q1, q2;
   
    for (int i = 0; i < system.molecules.size(); i++) {
    for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
    for (int k = i+1; k < system.molecules.size(); k++) {
    for (int l =0; l < system.molecules[k].atoms.size(); l++) { 
        if (!(system.molecules[i].atoms[j].C == 0 || system.molecules[k].atoms[l].C ==0)) {
        r = getDistance(system, i, j, k, l);
        q1 = system.molecules[i].atoms[j].C * system.constants.E2REDUCED;
        q2 = system.molecules[k].atoms[l].C * system.constants.E2REDUCED;
        
        potential += q1*q2/r;
        } // end if nonzero charges
    } // end l
    } // end k
    } // end j
    } // end i
    return potential;
}



