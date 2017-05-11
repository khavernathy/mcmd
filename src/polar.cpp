#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "thole_iterative.cpp"

#define OneOverSqrtPi 0.56418958354

using namespace std;

void makeAtomMap(System &system) {
    //int count =0;
    int i,j;
    // delete all elements from the atommap!! (i wasn't doing this before. This was the bug)
    system.atommap.clear();

    vector<int> local = vector<int>(2);
    //int v[2];
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            local = {i,j};
            system.atommap.push_back(local);

        //printf("system.atommap[%i] = {%i,%i}\n", count, system.atommap[count][0], system.atommap[count][1]);
          //  count++;
        }
    }
    //printf("SIZE OF ATOM MAP: %i\n", (int)system.atommap.size());
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
    int i,j;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            //if (isfinite(system.molecules[i].atoms[j].dipole_rrms))
            if (system.molecules[i].atoms[j].dipole_rrms != system.molecules[i].atoms[j].dipole_rrms) 
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
    system.last.thole_total_atoms = system.constants.total_atoms;
    N = 3*system.last.thole_total_atoms;
    dN = N-oldN;

    //printf("oldN: %i     N: %i     dN: %i\n",oldN,N,dN);

    if (system.proto.size() > 1) makeAtomMap(system); // don't know why, but for LJESPOLAR this is needed for multi-sorbate simulations.

    if(!dN) { return; }

    // re-make the AtomMap if needed
    makeAtomMap(system); // re-calculates unique indices for atoms

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

    int i, j, ii, jj, N, p, q;
    int w, x, y, z;
    double damp1=0, damp2=0; //, wdamp1=0, wdamp2=0; // v, s; //, distancesp[3], rp;
    double r, r2, ir3, ir5, ir=0;
    const double rcut = system.pbc.cutoff;
    //const double rcut2=rcut*rcut;
    //const double rcut3=rcut2*rcut;
    const double l=system.constants.polar_damp;
    const double l2=l*l;
    const double l3=l2*l;
    double explr; //exp(-l*r)
    const double explrcut = exp(-l*rcut);
    const double MAXVALUE = 1.0e40;
    N = (int)system.constants.total_atoms;

    system.checkpoint("in thole_amatrix() --> zeroing out");
    zero_out_amatrix(system,N);
    //print_matrix(system, 3*N, system.constants.A_matrix);
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
        w = system.atommap[i][0]; x = system.atommap[i][1];
        for(j = (i + 1);  j < N; j++) {
            jj = j*3;
            y = system.atommap[j][0]; z = system.atommap[j][1];

            //printf("i %i j %i ======= w %i x %i y %i z %i \n",i,j,w,x,y,z);

            double* distances = getDistanceXYZ(system, w,x,y,z);
            r = distances[3];
            // this on-the-spot distance calculator works, but the new method
            // below does not work, even though everywhere else, it does...
            //rp = system.pairs[w][x][y][z].r;
            //for (int n=0;n<3;n++) distancesp[n] = system.pairs[w][x][y][z].d[n];
            r2 = r*r;

            //printf("r: %f; rp: %f\n", r, rp);

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
    int i,j,k,l,p;
    const double SMALL_dR = 1e-12;
    double r, rr; //r and 1/r (reciprocal of r)
    const double R = system.pbc.cutoff;
    const double rR = 1./R;
    //used for polar_wolf_alpha (aka polar_wolf_damp)
    const double a = system.constants.polar_wolf_alpha; //, distances[3];
    const double erR=erfc(a*R); //complementary error functions
    const double cutoffterm = (erR*rR*rR + 2.0*a*OneOverSqrtPi*exp(-a*a*R*R)*rR);
    double bigmess=0;


    // first zero-out field vectors
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (p=0; p<3; p++) {
                system.molecules[i].atoms[j].efield[p] = 0;
                system.molecules[i].atoms[j].efield_self[p] = 0;
            }
        }
    }
    

    for(i=0; i<system.molecules.size(); i++) {
        for(j=0; j<system.molecules[i].atoms.size(); j++) {
            for(k=i+1; k<system.molecules.size(); k++) { // molecules not allowed to self-polarize
            for (l=0; l<system.molecules[k].atoms.size(); l++) {

                if ( system.molecules[i].frozen && system.molecules[k].frozen ) continue; //don't let the MOF polarize itself

                double* distances = getDistanceXYZ(system, i,j,k,l);
                r = distances[3];
                //r = system.pairs[i][j][k][l].r;
                //for (int n=0;n<3;n++) distances[n] = system.pairs[i][j][k][l].d[n];

                if((r - SMALL_dR  < system.pbc.cutoff) && (r != 0.)) {
                    rr = 1./r;

                    if ( a != 0 )   
                        bigmess=(erfc(a*r)*rr*rr+2.0*a*OneOverSqrtPi*exp(-a*a*r*r)*rr);
                
                    for ( p=0; p<3; p++ ) {
                        //see JCP 124 (234104)
                        if ( a == 0 ) {

                            // the commented-out charge=0 check here doesn't save time really.
                            //if (system.molecules[k].atoms[l].C != 0) {
                                system.molecules[i].atoms[j].efield[p] += 
                                (system.molecules[k].atoms[l].C)*
                                (rr*rr-rR*rR)*distances[p]*rr;
                            //}
                            //if (system.molecules[i].atoms[j].C != 0) {
                                system.molecules[k].atoms[l].efield[p] -= 
                                (system.molecules[i].atoms[j].C )*
                                (rr*rr-rR*rR)*distances[p]*rr;
                            //}

                        } else {
                            //if (system.molecules[k].atoms[l].C != 0) {
                                system.molecules[i].atoms[j].efield[p] +=
                                (system.molecules[k].atoms[l].C )*
                                (bigmess-cutoffterm)*distances[p]*rr;
                            //}
                            //if (system.molecules[i].atoms[j].C != 0) {
                                system.molecules[k].atoms[l].efield[p] -= 
                                (system.molecules[i].atoms[j].C )*
                                (bigmess-cutoffterm)*distances[p]*rr;
                            //}
                         }
                      //      printf("efield[%i]: %f\n", p,system.molecules[i].atoms[j].efield[p]);

                    } // end p

                } //cutoff 
            } // end l
            }  // end k
        } // end j
    } // end i

    /*
    printf("THOLE ELECTRIC FIELD: \n");
    for (int i=0; i<system.molecules.size(); i++)
        for (int j=0; j<system.molecules[i].atoms.size(); j++) 
            printf("ij efield: %f %f %f\n", system.molecules[i].atoms[j].efield[0], system.molecules[i].atoms[j].efield[1], system.molecules[i].atoms[j].efield[2]);
    printf("=======================\n");
    */
    return;


} // end thole_field()

void thole_field_nopbc(System &system) {
    int p, i, j, k,l;
    double r;
    const double SMALL_dR = 1e-12; //, distances[3];

    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (k=i+1; k<system.molecules.size(); k++) {
                for (l=0; l<system.molecules[k].atoms.size(); l++) {
                    if (system.molecules[i].frozen && system.molecules[k].frozen) continue;
                    double* distances = getDistanceXYZ(system,i,j,k,l);
                    r = distances[3];
                    //r = system.pairs[i][j][k][l].r;
                    //for (int n=0;n<3;n++) distances[n] = system.pairs[i][j][k][l].d[n];

                    if ( (r-SMALL_dR < system.pbc.cutoff) && (r != 0.)) {
                        for (p=0; p<3; p++) {
                            system.molecules[i].atoms[j].efield[p] += system.molecules[k].atoms[l].C * distances[p]/(r*r*r);
                            system.molecules[k].atoms[l].efield[p] -= system.molecules[i].atoms[j].C * distances[p]/(r*r*r);
                        }
                    }
                }
            }
        }
    }
} // end thole_field_nopbc

// =========================== POLAR POTENTIAL ========================
double polarization(System &system) {

    // POLAR ITERATIVE METHOD IS WHAT I USE.
    // THERE ARE OTHERS, E.G. MATRIX INVERSION OR FULL EWALD
    // MPMC CAN DO THOSE TOO, BUT WE ALMOST ALWAYS USE ITERATIVE.
    double potential; int i,j,num_iterations;

    // zero out everything on the atoms that pertains to polarness.
    // except static vars (like charge, polarizability, name, blah)
    // hopefully this fixes the issue but it will be redundant.
    // it did not fix the issue but i'm leaving it 
/*
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            system.molecules[i].atoms[j].rank_metric=0;
            for (int p=0; p<3; p++) {
                system.molecules[i].atoms[j].dip[p] =0;
                system.molecules[i].atoms[j].newdip[p]=0;
                system.molecules[i].atoms[j].olddip[p]=0;
                system.molecules[i].atoms[j].efield[p]=0;
                system.molecules[i].atoms[j].efield_self[p]=0;
                system.molecules[i].atoms[j].efield_induced[p]=0;
                system.molecules[i].atoms[j].efield_induced_change[p]=0;
            }
            system.molecules[i].atoms[j].dipole_rrms=0;
        }
    }
*/

    // 00) RESIZE THOLE A MATRIX IF NEEDED
    if (system.constants.ensemble == ENSEMBLE_UVT) {
        thole_resize_matrices(system);
    }
    system.checkpoint("done resizing thole matrix.");

    system.checkpoint("running thole_amatrix().");
    // 0) MAKE THOLE A MATRIX
    thole_amatrix(system); // ***this function also makes the i,j -> single-index atommap.
    system.checkpoint("done running thole_amatrix(). Running thole_field()");

    // 1) CALCULATE ELECTRIC FIELD AT EACH SITE
    if (system.constants.mc_pbc)
        thole_field(system);
    else
        thole_field_nopbc(system); // maybe in wrong place? doubt it. 4-13-17
    system.checkpoint("done with thole_field(). Doing dipole iterations");


    // 2) DO DIPOLE ITERATIONS
    num_iterations = thole_iterative(system);    
//    printf("num_iterations = %f\n", (double)num_iterations);
    system.stats.polar_iterations = (double)num_iterations;
    system.constants.dipole_rrms = get_dipole_rrms(system);
    system.checkpoint("done with dipole iters. Calculating polarization energy");


    // 3) CALCULATE POLARIZATION ENERGY 1/2 mu*E
    potential=0;
    for (i=0; i<system.molecules.size(); i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            potential += 
                dddotprod(system.molecules[i].atoms[j].dip, system.molecules[i].atoms[j].efield);

            if (system.constants.polar_palmo) {
                potential += dddotprod(system.molecules[i].atoms[j].dip, system.molecules[i].atoms[j].efield_induced_change);
            }
        }
//        printf("molecule %i dip %f efield %f indchang %f\n", i, system.molecules[i].atoms[0].dip[0], system.molecules[i].atoms[0].efield[0], system.molecules[i].atoms[0].efield_induced_change[0]);
    }
    system.checkpoint("POLARIZATION POTENTIAL CALCULATED.");
    
    potential *= -0.5;
    //printf("POLARIZATION POTENTIAL: %f\n",potential);
    return potential;

} // end polarization() function



void polarization_force(System &system) {
    // gets force on atoms due to dipoles calculated before (via iterative method)
    // TODO
    // adam sent a code snippet that should help on slack
}
