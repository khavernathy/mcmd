#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include <thole_iterative.cpp>

using namespace std;

void makeAtomMap(System &system) {
    int count =0;
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            vector<int> local = {i,j};
            
            system.atommap.push_back(local);

        //printf("system.atommap[%i] = {%i,%i}\n", count, system.atommap[count][0], system.atommap[count][1]);
            count++;
        }
    }
    return;  
}

void print_matrix(System &system, int N, double **matrix) {
    system.checkpoint("PRINTING MATRIX");

    int i,j;
    printf("\n");
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            printf("%.3f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void zero_out_amatrix(System &system, int N) {
    int i,j;
    for (i=0; i<3*N; i++)
        for (j=0; j<3*N; j++)
            system.constants.A_matrix[i][j] = 0;
    return;
}

double get_dipole_rrms (System &system) {
    double N, dipole_rrms;
    dipole_rrms = N = 0;

    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            if (isfinite(system.molecules[i].atoms[j].dipole_rrms))
                dipole_rrms += system.molecules[i].atoms[j].dipole_rrms;
            N++;
        }
    }

    return dipole_rrms / N;
}

/* for uvt runs, resize the A matrix */
void thole_resize_matrices(System &system) {

    int i, N, dN, oldN;

    /* determine how the number of atoms has changed and realloc matrices */
    oldN = 3*system.last.thole_total_atoms; //will be set to zero if first time called
    N = 3*system.constants.total_atoms;
    dN = N-oldN;
    system.last.thole_total_atoms = system.constants.total_atoms;

    //printf("oldN: %i     N: %i     dN: %i\n",oldN,N,dN);

    if(!dN) { return; }

    // grow A matricies by free/malloc (to prevent fragmentation)
    //free the A matrix
    for (i=0; i < oldN; i++) free(system.constants.A_matrix[i]);
    free(system.constants.A_matrix);


    //(RE)allocate the A matrix
    system.constants.A_matrix= (double **) calloc(N,sizeof(double*));
//    memnullcheck(system.constants.A_matrix,N*sizeof(double*), __LINE__-1, __FILE__);

    for (i=0; i< N; i++ ) {
        system.constants.A_matrix[i]= (double *) malloc(N*sizeof(double));
  //      memnullcheck(system.constants.A_matrix[i],N*sizeof(double), __LINE__-1, __FILE__);
    }

     return;
}

/* calculate the dipole field tensor */
void thole_amatrix(System &system) {

    makeAtomMap(system); // re-calculates unique indices for atoms

    int i, j, ii, jj, N, p, q;
    int w, x, y, z;
    double damp1=0, damp2=0, wdamp1=0, wdamp2=0, v, s;
    double r, r2, ir3, ir5, ir=0;
    double rcut, rcut2, rcut3;
    rcut = system.pbc.cutoff;
    rcut2 = rcut*rcut; rcut3 = rcut2*rcut;
    double l, l2, l3;
    l = system.constants.polar_damp;
    l2 = l*l; l3 = l2*l;
    double explr; //exp(-l*r)
    double explrcut = exp(-l*rcut);
    double MAXVALUE = 1.0e40;
    N = (int)system.constants.total_atoms;

    system.checkpoint("in thole_amatrix() --> zeroing out");
    zero_out_amatrix(system,N);
    system.checkpoint("done with zero_out_amatrix()");

    //system.checkpoint("setting diagonals in A");
    /* set the diagonal blocks */
    for(i = 0; i < N; i++) {
        ii = i*3;
        w = system.atommap[i][0];
        x = system.atommap[i][1];

        for(p = 0; p < 3; p++) {
            if(system.molecules[w].atoms[x].polar != 0.0)
                system.constants.A_matrix[ii+p][ii+p] = 1.0/system.molecules[w].atoms[x].polar;
            else
                system.constants.A_matrix[ii+p][ii+p] = MAXVALUE;
        }
    }   
    //system.checkpoint("done setting diagonals in A");
   
    //system.checkpoint("starting Tij loop"); 
    /* calculate each Tij tensor component for each dipole pair */
    for(i = 0; i < (N - 1); i++) {
        ii = i*3;
        for(j = (i + 1);  j < N; j++) {
            jj = j*3;

            w = system.atommap[i][0]; x = system.atommap[i][1];
            y = system.atommap[j][0]; z = system.atommap[j][1];

            //printf("i %i j %i ======= w %i x %i y %i z %i \n",i,j,w,x,y,z);

            double* distances = getDistanceXYZ(system, w,x,y,z);
            r = distances[3];
            r2 = r*r;

            //system.checkpoint("got r.");
            //printf("distances: x %f y %f z %f r %f\n", distances[0], distances[1], distances[2], r);

            /* inverse displacements */
            if(r == 0.)
                ir3 = ir5 = MAXVALUE;
            else {
                ir = 1.0/r;
                ir3 = ir*ir*ir;
                ir5 = ir3*ir*ir;
            }

            //evaluate damping factors
                    explr = exp(-l*r);
                    damp1 = 1.0 - explr*(0.5*l2*r2 + l*r + 1.0);
                    damp2 = damp1 - explr*(l3*r2*r/6.0);

            //system.checkpoint("got damping factors.");

           // system.checkpoint("buildling tensor.");
            /* build the tensor */
            for(p = 0; p < 3; p++) {
                for(q = 0; q < 3; q++) {

                    system.constants.A_matrix[ii+p][jj+q] = -3.0*distances[p]*distances[q]*damp2*ir5;

                    /* additional diagonal term */
                    if(p == q) {
                        system.constants.A_matrix[ii+p][jj+q] += damp1*ir3;
                    }   
                    //printf("A[%i][%i] = %f\n", ii+p, jj+q, system.constants.A_matrix[ii+p][jj+q]);
                }
            }
            //system.checkpoint("done building tensor. setting lower half.");
    
            //printf("\n");
            /* set the lower half of the tensor component */
            for(p = 0; p < 3; p++) {
                for(q = 0; q < 3; q++) {
                   //printf("setting A_matrix[%i][%i] = A_matrix[%i][%i] = %f at 3*N = %i, i=%i, j=%i\n", jj+p, ii+q, ii+p, jj+q, system.constants.A_matrix[jj+p][ii+q], 3*N,i,j);
                    system.constants.A_matrix[jj+p][ii+q] = system.constants.A_matrix[ii+p][jj+q];
                } // end q
            } // end p
        } /* end j */
    } /* end i */
    //print_matrix(system, N*3, system.constants.A_matrix);
    //printf("\n===================================================\n\n");
    return;
}

void thole_field(System &system) {
    // wolf thole field

    // first zero-out field vectors
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            for (int p=0; p<3; p++) {
                system.molecules[i].atoms[j].efield[p] = 0;
                system.molecules[i].atoms[j].efield_self[p] = 0;
            }
        }
    }

    double OneOverSqrtPi = 1.0/sqrt(M_PI);
    int p; //dimensionality
    double r, rr; //r and 1/r (reciprocal of r)
    double R = system.pbc.cutoff;
    double rR = 1./R;
    //used for polar_wolf_alpha (aka polar_wolf_damp)
    double a = system.constants.polar_wolf_alpha;
    double erR; //complementary error functions
        erR=erfc(a*R);
    double cutoffterm = (erR*rR*rR + 2.0*a*OneOverSqrtPi*exp(-a*a*R*R)*rR);
    double bigmess=0;

    //printf("R = %f; a = %f; erR = %f; cutoffterm = %f\n", R, a, erR, cutoffterm);

    for(int i=0; i<system.molecules.size(); i++) {
        for(int j=0; j<system.molecules[i].atoms.size(); j++) {
            for(int k=i+1; k<system.molecules.size(); k++) { // molecules not allowed to self-polarize
            for (int l=0; l<system.molecules[k].atoms.size(); l++) {

                if ( system.molecules[i].MF == "F" && system.molecules[k].MF == "F" ) continue; //don't let the MOF polarize itself

                double* distances = getDistanceXYZ(system, i,j,k,l);
                r = distances[3];

                if((r  < system.pbc.cutoff) && (r != 0.)) {
                    rr = 1./r;

                    if ( a != 0 )   
                        bigmess=(erfc(a*r)*rr*rr+2.0*a*OneOverSqrtPi*exp(-a*a*r*r)*rr);
                
                    for ( p=0; p<3; p++ ) {
                        //see JCP 124 (234104)
                        if ( a == 0 ) {
                            system.molecules[i].atoms[j].efield[p] += 
                            (system.molecules[k].atoms[l].C)*
                            (rr*rr-rR*rR)*distances[p]*rr;
                        
                            system.molecules[k].atoms[l].efield[p] -= 
                            (system.molecules[i].atoms[j].C )*
                            (rr*rr-rR*rR)*distances[p]*rr;

                        } else {
                            system.molecules[i].atoms[j].efield[p] +=
                            (system.molecules[k].atoms[l].C )*
                            (bigmess-cutoffterm)*distances[p]*rr;

                            system.molecules[k].atoms[l].efield[p] -= 
                            (system.molecules[i].atoms[j].C )*
                            (bigmess-cutoffterm)*distances[p]*rr;
                         }
                      //      printf("efield[%i]: %f\n", p,system.molecules[i].atoms[j].efield[p]);

                    } // end p

                } //cutoff 
            } // end l
            }  // end k
        } // end j
    } // end i

    return;

/*
printf("THOLE FIELD RESULT: \n");
for (int i=0; i<system.molecules.size(); i++) {
    for (int j=0; j<system.molecules[i].atoms.size(); j++) {
        printf("
    }
}
*/

} // end thole_field()


// =========================== POLAR POTENTIAL ========================
double polarization(System &system) {

    // POLAR ITERATIVE METHOD IS WHAT I USE.
    // THERE ARE OTHERS, E.G. MATRIX INVERSION OR FULL EWALD
    // MPMC CAN DO THOSE TOO, BUT WE ALMOST ALWAYS USE ITERATIVE.
    double potential; int num_iterations;

    // initialize potentials to zero
	for (int j=0; j<system.molecules.size(); j++) {
		for (int i = 0; i < system.molecules[j].atoms.size(); i++) {
			system.molecules[j].atoms[i].V = 0.0;
		} // end atom loop i
	} // end molecule loop j
    system.checkpoint("done with zero-initiallization of e-fields");

    // 00) RESIZE THOLE A MATRIX IF NEEDED
    if (system.constants.ensemble == "uvt") {
        thole_resize_matrices(system);
    }
    system.checkpoint("done resizing thole matrix.");

    system.checkpoint("running thole_amatrix().");
    // 0) MAKE THOLE A MATRIX
    thole_amatrix(system); // this function also makes the i,j -> single-index atommap.
    system.checkpoint("done running thole_amatrix(). Running thole_field()");

    // 1) CALCULATE ELECTRIC FIELD AT EACH SITE
    thole_field(system);
    system.checkpoint("done with thole_field(). Doing dipole iterations");


    // 2) DO DIPOLE ITERATIONS
    num_iterations = thole_iterative(system);    
    //printf("num_iterations = %f\n", (double)num_iterations);
    system.stats.polar_iterations = (double)num_iterations;
    system.constants.dipole_rrms = get_dipole_rrms(system);
    system.checkpoint("done with dipole iters. Calculating polarization energy");


    // 3) CALCULATE POLARIZATION ENERGY 1/2 mu*E
    potential=0;
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            potential += (
            (system.molecules[i].atoms[j].dip[0] * system.molecules[i].atoms[j].efield[0]) +
            (system.molecules[i].atoms[j].dip[1] * system.molecules[i].atoms[j].efield[1]) +
            (system.molecules[i].atoms[j].dip[2] * system.molecules[i].atoms[j].efield[2])
            );
            
            if (system.constants.polar_palmo) { 
                potential += (
                    system.molecules[i].atoms[j].dip[0] * system.molecules[i].atoms[j].efield_induced_change[0] +
                    system.molecules[i].atoms[j].dip[1] * system.molecules[i].atoms[j].efield_induced_change[1] +
                    system.molecules[i].atoms[j].dip[2] * system.molecules[i].atoms[j].efield_induced_change[2]
                );
            }
        }
    }
    system.checkpoint("POLARIZATION POTENTIAL CALCULATED.");
    
    potential *= -0.5;
    //printf("POLARIZATION POTENTIAL: %f\n",potential);
    return potential;

} // end polarization() function
