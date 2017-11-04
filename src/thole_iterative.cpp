#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#define MAX_ITERATION_COUNT 128

//set dipoles to alpha*E_static
void init_dipoles (System &system) {
	unsigned int i, j, p;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (p=0; p<3; p++) {
                system.molecules[i].atoms[j].dip[p] =
                system.molecules[i].atoms[j].polar *
                (system.molecules[i].atoms[j].efield[p] + system.molecules[i].atoms[j].efield_self[p]);
                // improve convergence
                system.molecules[i].atoms[j].dip[p] *= system.constants.polar_gamma;
            }
        }
    }
    return;
}

void contract_dipoles (System &system, int * ranked_array ) {
    unsigned int i, j, ii, jj, p, index, n, ti, tj, tk, tl;
      
    for(i = 0; i < system.constants.total_atoms; i++) {
        index = ranked_array[i]; //do them in the order of the ranked index
        ii = index*3;

        ti = system.atommap[index][0]; tj = system.atommap[index][1];

        if ( system.molecules[ti].atoms[tj].polar == 0 ) { //if not polar
            for (n=0; n<3; n++) {
                system.molecules[ti].atoms[tj].newdip[n] = 0;
                system.molecules[ti].atoms[tj].dip[n] = 0;
            }
            continue;
        }

        for(j = 0; j < system.constants.total_atoms; j++) {
            jj = j*3;
            if(index != j) { 
                tk = system.atommap[j][0]; tl = system.atommap[j][1];
                for(p = 0; p < 3; p++) {
                    for (int q=0; q<3; q++) {
                        // account for the 1/2 matrix
                        if (!system.constants.full_A_matrix_option) {
                            if (j>index) {system.molecules[ti].atoms[tj].efield_induced[p] -=
                                system.constants.A_matrix[jj+p][ii+q] * system.molecules[tk].atoms[tl].dip[q];
                            } else { system.molecules[ti].atoms[tj].efield_induced[p] -=
                                system.constants.A_matrix[ii+p][jj+q] * system.molecules[tk].atoms[tl].dip[q];
                            }
                        // full matrix computation
                        } else {
                            system.molecules[ti].atoms[tj].efield_induced[p] -=
                                system.constants.A_matrix_old[ii+p][jj+q] * system.molecules[tk].atoms[tl].dip[q];
                        }
                    } // end q
                    // correct old matrix, written out dot prod
                    //(system.constants.A_matrix[ii+p]+jj)[0] * system.molecules[tk].atoms[tl].dip[0] + 
                    //(system.constants.A_matrix[ii+p]+jj)[1] * system.molecules[tk].atoms[tl].dip[1] +
                    //(system.constants.A_matrix[ii+p]+jj)[2] * system.molecules[tk].atoms[tl].dip[2];
                } // end p
            }
        } /* end j */

        /* dipole is the sum of the static and induced parts */
        for(p = 0; p < 3; p++) {
            system.molecules[ti].atoms[tj].newdip[p] = system.molecules[ti].atoms[tj].polar *
            (system.molecules[ti].atoms[tj].efield[p] + 
                system.molecules[ti].atoms[tj].efield_self[p] + 
                system.molecules[ti].atoms[tj].efield_induced[p]);
        
            if (system.constants.polar_gs || system.constants.polar_gs_ranked) {
                system.molecules[ti].atoms[tj].dip[p] = system.molecules[ti].atoms[tj].newdip[p];
            }
        }
    } /* end matrix multiply */

    return;
}

void calc_dipole_rrms (System &system) {
    unsigned int i, j, p;
    double carry;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            // mean square distance
            system.molecules[i].atoms[j].dipole_rrms=0;
            for (p=0; p<3; p++) {
                carry = system.molecules[i].atoms[j].newdip[p] - system.molecules[i].atoms[j].olddip[p];
                system.molecules[i].atoms[j].dipole_rrms += carry*carry;
            }
                // normalize
                system.molecules[i].atoms[j].dipole_rrms /= dddotprod(system.molecules[i].atoms[j].newdip, system.molecules[i].atoms[j].newdip);
                system.molecules[i].atoms[j].dipole_rrms = sqrt(system.molecules[i].atoms[j].dipole_rrms);

                if (system.molecules[i].atoms[j].dipole_rrms != system.molecules[i].atoms[j].dipole_rrms)
                    system.molecules[i].atoms[j].dipole_rrms=0;
        }
    }
    return;
}

int are_we_done_yet (System &system, int iteration_counter ) {
	unsigned int i, p, ti, tj;
	double allowed_sqerr, error;
    int N = system.constants.total_atoms;	

	if (system.constants.polar_precision == 0.0) {	/* DEFAULT ... by fixed iteration ... */
		if(iteration_counter != system.constants.polar_max_iter)
			return 1;
	} 
	else { /* ... or by dipole precision */
		allowed_sqerr = system.constants.polar_precision*system.constants.polar_precision*
                    system.constants.DEBYE2SKA*system.constants.DEBYE2SKA;
		
        for(i = 0; i < N; i++) { //check the change in each dipole component
            ti= system.atommap[i][0]; tj = system.atommap[i][1];
            for(p = 0; p < 3; p++) {
                error = system.molecules[ti].atoms[tj].newdip[p] - 
                    system.molecules[ti].atoms[tj].olddip[p];
                if(error*error > allowed_sqerr)
                    return 1; //we broke tolerance
            }
        }


    }
	return 0;
}

void palmo_contraction (System &system, int * ranked_array ) {
    unsigned int i, j, ii, jj, index, p, q,ti,tj, tk, tl;
    int N = system.constants.total_atoms;

    /* calculate change in induced field due to this iteration */
    for(i = 0; i < N; i++) {
        index = ranked_array[i];
        ii = index*3;

        ti = system.atommap[index][0]; tj = system.atommap[index][1];

        for (p=0; p<3; p++ )
            system.molecules[ti].atoms[tj].efield_induced_change[p] = -system.molecules[ti].atoms[tj].efield_induced[p];

        for(j = 0; j < N; j++) {
            jj = j*3;
            if(index != j) {
                tk = system.atommap[j][0]; tl = system.atommap[j][1];
                for(p = 0; p < 3; p++) {
                    for (q=0; q<3; q++) {
                        // account for the 1/2 matrix
                        if (!system.constants.full_A_matrix_option) {
                            if (j>index) {system.molecules[ti].atoms[tj].efield_induced_change[p] -=
                                system.constants.A_matrix[jj+p][ii+q] * system.molecules[tk].atoms[tl].dip[q];
                            } else { system.molecules[ti].atoms[tj].efield_induced_change[p] -=
                                system.constants.A_matrix[ii+p][jj+q] * system.molecules[tk].atoms[tl].dip[q];
                            }
                        // full matrix
                        } else {
                            system.molecules[ti].atoms[tj].efield_induced_change[p] -=
                                system.constants.A_matrix_old[ii+p][jj+q] * system.molecules[tk].atoms[tl].dip[q];
                        }

                    // old correct full A matrix
                    //system.molecules[ti].atoms[tj].efield_induced_change[p] -=
                    //dddotprod((system.constants.A_matrix[ii+p]+jj), system.molecules[tk].atoms[tl].dip); 
                } // end q 
                }// end p
            } 
        }
    }
    return;
}


void update_ranking (System &system, int * ranked_array ) {
	unsigned int i, j, k,l,sorted, tmp, trk, trl, trk1, trl1, rankedj, rankedj1;
	const double rmin = system.constants.polar_rmin*1.5; 
    double r;

    int N = system.constants.total_atoms;

    // initialize rank metrics
    for (i=0; i<system.molecules.size(); i++) 
        for (j=0; j<system.molecules[i].atoms.size(); j++)
            system.molecules[i].atoms[j].rank_metric = 0;
    

    // set rank metrics
    int hits=0;
    int paircount=0;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            if (system.molecules[i].atoms[j].polar == 0) continue;
            for (k=i;k<system.molecules.size(); k++) {
                for (l=0; l<system.molecules[k].atoms.size(); l++) {
                    if (system.molecules[k].atoms[l].polar ==0) continue;
                    if ((i==k && j>=l)) continue; // don't do self, or double count.
                    double* distances = getDistanceXYZ(system, i,j,k,l);
                    r = distances[3]; // in MPMC, this is NOT rimg, as it is here. There may be a reason for that.
                    if (r <= rmin) {
                        system.molecules[i].atoms[j].rank_metric += 1.0;
                        system.molecules[k].atoms[l].rank_metric += 1.0;
                        hits++;
                    }
                    paircount++;   
                } // l
            } // k
        } // j
    } // i end rank metric determination
    
    /* rank the dipoles by bubble sort */
        for(i = 0; i < N; i++) {
            for(j = 0, sorted = 1; j < (N-1); j++) {

                rankedj = ranked_array[j];
                trk = system.atommap[rankedj][0]; trl = system.atommap[rankedj][1];
    
                rankedj1 = ranked_array[j+1];
                trk1 = system.atommap[rankedj1][0]; trl1 = system.atommap[rankedj1][1];

                if(system.molecules[trk].atoms[trl].rank_metric < system.molecules[trk1].atoms[trl1].rank_metric) {
                    sorted = 0;
                    tmp = ranked_array[j];
                    ranked_array[j] = ranked_array[j+1];
                    ranked_array[j+1] = tmp;
                }
            }
            if(sorted) break;
        }
	return;
}


/* iterative solver of the dipole field tensor */
/* returns the number of iterations required */
int thole_iterative(System &system) {
    unsigned int i, j, N, p, ti, tj;
    unsigned int iteration_counter, keep_iterating;
    int *ranked_array;

    N = system.constants.total_atoms;

    /* array for ranking */
    ranked_array = (int *) calloc(N, sizeof(int));
    for(i = 0; i < N; i++) ranked_array[i] = i;

    //set all dipoles to alpha*E_static * polar_gamma
    init_dipoles(system);


    /* iterative solver of the dipole field equations */
    keep_iterating = 1;
    iteration_counter = 0;
    while (keep_iterating) {
        iteration_counter++;

        /* divergence detection */
        /* if we fail to converge, then return dipoles as alpha*E */
        if(iteration_counter >= MAX_ITERATION_COUNT && system.constants.polar_precision) // && system.constants.polar_precision != 0) {
        {
            for(i = 0; i < N; i++) {
                ti = system.atommap[i][0]; tj = system.atommap[i][1];
                for(p = 0; p < 3; p++) {
                    system.molecules[ti].atoms[tj].dip[p] = 
                    system.molecules[ti].atoms[tj].polar * 
                    (system.molecules[ti].atoms[tj].efield[p] + 
                        system.molecules[ti].atoms[tj].efield_self[p]);
                    system.molecules[ti].atoms[tj].efield_induced_change[p] = 0.0; //so we don't break palmo
                }
            }
            //set convergence failure flag
            system.constants.iter_success = 1;
            printf("POLAR CONVERGENCE FAILURE\n");           
 
            free(ranked_array);
            return(iteration_counter);
        }

        //zero out induced e-field
        for (i=0; i<system.molecules.size(); i++) {
            for (j=0; j<system.molecules[i].atoms.size(); j++) {
                for (p=0; p<3; p++)
                    system.molecules[i].atoms[j].efield_induced[p] = 0;
            }
        }

        //save the current dipole information if we want to calculate precision (or if needed for relaxation)
        if ( system.constants.polar_rrms || system.constants.polar_precision > 0)  { 
            for(i = 0; i < N; i++) {
                ti = system.atommap[i][0]; tj = system.atommap[i][1];
                for(p = 0; p < 3; p++) 
                    system.molecules[ti].atoms[tj].olddip[p] = system.molecules[ti].atoms[tj].dip[p];
            }
        }

        // contract the dipoles with the field tensor (gauss-seidel/gs-ranked optional)
        contract_dipoles(system, ranked_array);

        if ( system.constants.polar_rrms || system.constants.polar_precision > 0 )
            calc_dipole_rrms(system);

        /* determine if we are done... */
        keep_iterating = are_we_done_yet(system, iteration_counter);

       // if we would be finished, contract once more to get the next induced field for palmo
        if (system.constants.polar_palmo && !keep_iterating) {
            palmo_contraction(system, ranked_array);
        } 

        //new gs_ranking if needed
        if ( system.constants.polar_gs_ranked && keep_iterating )
            update_ranking(system, ranked_array);

        /* save the dipoles for the next pass */
        for(i = 0; i < N; i++) {
            ti = system.atommap[i][0]; tj= system.atommap[i][1];
            for(p = 0; p < 3; p++) {
                    system.molecules[ti].atoms[tj].dip[p] = system.molecules[ti].atoms[tj].newdip[p];
            }
        }

    } //end iterate
    free(ranked_array);

    /* return the iteration count */
    return(iteration_counter);
}


