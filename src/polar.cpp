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
    //system.checkpoint("PRINTING MATRIX");
    int i,j;
    printf("\n");
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            printf("%.3f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /* This is for the 1/2 A-matrix, deactivated for now. 
    // NEW ***
    int blocksize=3, inc=0;
    printf("\n");
    for (i=0; i<N; i++) {
        for (j=0; j<blocksize; j++) {
            printf("%.3f ", matrix[i][j]);
        }
        printf("\n");
        inc++;
        if (inc%3==0) blocksize+=3;
    }
    */
}

void zero_out_amatrix(System &system, int N) {
    int i,j;
    // half matrix
    if (!system.constants.full_A_matrix_option) {
        int blocksize=3, inc=0;
        for (i=0; i<3*N; i++) {
            for (j=0; j<blocksize; j++) {
                system.constants.A_matrix[i][j] = 0;    
            }
            inc++;
            if (inc%3==0) blocksize+=3;
        }
    // full matrix
    } else {
        for (i=0;i<3*N;i++) {
            for (j=0;j<3*N;j++) {
                system.constants.A_matrix_old[i][j]=0;
            }
        }
    }
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

    //just make the atom map of indices no matter what.
    // for whatever reason things get buggy when I try to
    // minimize calls to this function. It's not expensive anyway.
    makeAtomMap(system);

    /* determine how the number of atoms has changed and realloc matrices */
    oldN = 3*system.last.thole_total_atoms; //will be set to zero if first time called
    system.last.thole_total_atoms = system.constants.total_atoms;
    N = 3*system.last.thole_total_atoms;
    dN = N-oldN;

    //printf("oldN: %i     N: %i     dN: %i\n",oldN,N,dN);

    if(!dN) { return; }

    // 1/2 matrix
    if (!system.constants.full_A_matrix_option) {
        for (i=0; i< oldN; i++) free(system.constants.A_matrix[i]);
        free(system.constants.A_matrix);
        system.constants.A_matrix = (double **) calloc(N, sizeof(double*));
        int blocksize=3, inc=0;
        for (i=0; i<N; i++) {
            system.constants.A_matrix[i] = (double *) malloc(blocksize*sizeof(double));
            inc++;
            if (inc%3==0) blocksize+=3;
        }
    // full matrix
    } else {
        for (i=0; i<oldN; i++) free(system.constants.A_matrix_old[i]);
        free(system.constants.A_matrix_old);
        system.constants.A_matrix_old = (double **) calloc(N, sizeof(double*));
        for (i=0; i<N; i++) {
            system.constants.A_matrix_old[i] = (double *) malloc(N*sizeof(double));
        }
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
    double rmin = 1.0e40;

    system.checkpoint("in thole_amatrix() --> zeroing out");
    zero_out_amatrix(system,N);
    system.checkpoint("done with zero_out_amatrix()");

    system.checkpoint("setting diagonals in A");
    /* set the diagonal blocks */
    for(i = 0; i < N; i++) {
        ii = i*3;
        w = system.atommap[i][0];
        x = system.atommap[i][1];

        // 1/2 matrix
        if (!system.constants.full_A_matrix_option) {
            for (p=0; p<3; p++) {
                if (system.molecules[w].atoms[x].polar != 0.0)
                    system.constants.A_matrix[ii+p][ii+p] = 1.0/system.molecules[w].atoms[x].polar;
                else
                    system.constants.A_matrix[ii+p][ii+p] = MAXVALUE;
            }
        // full matrix
        } else {
            for (p=0;p<3;p++) {
                if (system.molecules[w].atoms[x].polar != 0.0)
                    system.constants.A_matrix_old[ii+p][ii+p] = 1.0/system.molecules[w].atoms[x].polar;
                else
                    system.constants.A_matrix_old[ii+p][ii+p] = MAXVALUE;
            }
        }

    }
    system.checkpoint("done setting diagonals in A");

    system.checkpoint("starting Tij loop");
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

            if (r<rmin && (system.molecules[w].atoms[x].polar!=0 && system.molecules[y].atoms[z].polar!=0))
                rmin=r; // for ranking, later


            // this on-the-spot distance calculator works, but the new method
            // below does not work, even though everywhere else, it does...
            //rp = system.pairs[w][x][y][z].r;
            //for (int n=0;n<3;n++) distancesp[n] = system.pairs[w][x][y][z].d[n];
            r2 = r*r;

            //printf("r: %f; rp: %f\n", r, rp);

            ////system.checkpoint("got r.");
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

            system.checkpoint("got damping factors.");

            system.checkpoint("buildling tensor.");
            /* build the tensor */
            // 1/2 MATRIX
            if (!system.constants.full_A_matrix_option) {
                for (p=0; p<3; p++) {
                    for (q=0; q<3; q++) {
                        system.constants.A_matrix[jj+p][ii+q] = -3.0*distances[p]*distances[q]*damp2*ir5;
                        // additional diagonal term
                        if (p==q)
                            system.constants.A_matrix[jj+p][ii+q] += damp1*ir3;
                    }
                }
            // full matrix
            } else {
                for (p=0; p<3; p++) {
                    for (q=0; q<3; q++) {
                        system.constants.A_matrix_old[ii+p][jj+q] = -3.0 * distances[p]*distances[q] * damp2 * ir5;
                        // additional diagonal term
                        if (p==q)
                            system.constants.A_matrix_old[ii+p][jj+q] += damp1*ir3;
                    }
                }
            } // end full matrix


            // fill in other half of full matrix if app.
            if (system.constants.full_A_matrix_option) {
                for (p=0; p<3; p++) {
                    for (q=0; q<3; q++) {
                        system.constants.A_matrix_old[jj+p][ii+q] = system.constants.A_matrix_old[ii+p][jj+q];
                    }
                }
            } // end full matrix
        } /* end j */
    } /* end i */

    


    system.constants.polar_rmin = rmin;
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

    // POLAR ITERATIVE METHOD FROM THOLE/APPLEQUIST IS WHAT I USE.
    // THERE ARE OTHERS, E.G. MATRIX INVERSION OR FULL EWALD
    // MPMC CAN DO THOSE TOO, BUT WE ALMOST ALWAYS USE ITERATIVE.
    double potential; int i,j,num_iterations;

    // 00) RESIZE THOLE A MATRIX IF NEEDED
    if (system.constants.ensemble == ENSEMBLE_UVT) { // uVT is the only ensemble that changes N
        thole_resize_matrices(system);
    }

    // 0) MAKE THOLE A MATRIX
    thole_amatrix(system); // ***this function also makes the i,j -> single-index atommap.

    // 1) CALCULATE ELECTRIC FIELD AT EACH SITE
    if (system.constants.mc_pbc)
        thole_field(system);
    else
        thole_field_nopbc(system); // maybe in wrong place? doubt it. 4-13-17


    // 2) DO DIPOLE ITERATIONS
    num_iterations = thole_iterative(system);
    system.stats.polar_iterations = (double)num_iterations;
    system.constants.dipole_rrms = get_dipole_rrms(system);


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

    }

    potential *= -0.5;
    return potential;

} // end polarization() function


void polarization_force(System &system) {
    // gets force on atoms due to dipoles calculated via iterative method
    int i,j,k,l,n;
    double common_factor, r, rinv, r2, r2inv, r3, r3inv, r5inv, r7inv;
    double x2,y2,z2,xy,yz,xz,x,y,z; // distance relations
    double udotu, ujdotr, uidotr; // dot products
    const double damp = system.constants.polar_damp;
    const double cc2inv = (1.0/system.pbc.cutoff)*(1.0/system.pbc.cutoff); // coulombic cutoff is same as LJ; this is -1*f_shift
    double f_local[3]; // temp forces
    double u_i[3], u_j[3]; // temp. dipoles
    double q_i,q_j; // temp. charges
    double t1, t2, t3, p1, p2, p3, p4, p5; // terms and prefactors

    if (system.constants.ensemble == ENSEMBLE_UVT) { // uVT is the only ensemble that changes N
        thole_resize_matrices(system);
    } 

    thole_amatrix(system); // fill in the A-matrix
    thole_field(system); // calculate electric field at each atom (assume PBC)
    int num_iterations = thole_iterative(system); // iteratively solve the dipoles
       system.stats.polar_iterations = (double)num_iterations;
       system.constants.dipole_rrms = get_dipole_rrms(system);
    
    // ready for forces; loop all atoms
    for (i=0; i<system.molecules.size(); i++) {
    for (j=0; j<system.molecules[i].atoms.size(); j++) {
        // initializers for atom
        q_i = system.molecules[i].atoms[j].C;
        for (n=0;n<3;n++) u_i[n] = system.molecules[i].atoms[j].dip[n];

        // loop pairs
        for (k=i+1; k<system.molecules.size(); k++) { // no self-molecule interactions
        for (l=0; l<system.molecules[k].atoms.size(); k++) {
            // there are 3 pairwise contributions to polar force:
            // (1) u_i -- q_j  ||  (2) u_j -- q_i  ||  (3) u_i -- u_j
            for (n=0;n<3;n++) f_local[n]=0;

            double* distances = getDistanceXYZ(system,i,j,k,l);
            r = distances[3];
            if (r > system.pbc.cutoff) continue; // only within r_cc
            x = distances[0]; y = distances[1]; z = distances[2];
            xy = x*y;
            yz = y*z;
            xz = x*z;
            x2 = x*x;
            y2 = y*y;
            z2 = z*z;
            r2 = r*r;
            r3 = r2*r;
            rinv = 1./r;
            r2inv = rinv*rinv;
            r3inv = r2inv*rinv;

            // (1) u_i -- q_j
            if (system.molecules[k].atoms[l].C != 0 && system.molecules[i].atoms[j].polar != 0) {
                q_j = system.molecules[k].atoms[l].C;
                common_factor = q_j*r3inv;
    
                f_local[0] += common_factor*((u_i[0]*(r2inv*(-2*x2 + y2 + z2) - cc2inv*(y2 + z2))) + (u_i[1]*(r2inv*(-3*x*y) + cc2inv*x*y)) + (u_i[2]*(r2inv*(-3*x*z) + cc2inv*x*z)));
     
                f_local[1] += common_factor*(u_i[0]*(r2inv*(-3*x*y) + cc2inv*x*y) + u_i[1]*(r2inv*(-2*y2 + x2 + z2) - cc2inv*(x2 + z2)) + u_i[2]*(r2inv*(-3*y*z) + cc2inv*y*z));
                
                f_local[2] += common_factor*(u_i[0]*(r2inv*(-3*x*z) + cc2inv*x*z) + u_i[1]*(r2inv*(-3*y*z) + cc2inv*y*z) + u_i[2]*(r2inv*(-2*z2 + x2 + y2) - cc2inv*(x2 + y2)));
            }

            // (2) u_j -- q_i
            if (q_i != 0 && system.molecules[k].atoms[l].polar != 0) {
                for (n=0;n<3;n++) u_j[n] = system.molecules[k].atoms[l].dip[n];
                common_factor = q_i*r3inv;

                f_local[0] -= common_factor*((u_j[0]*(r2inv*(-2*x2 + y2 + z2) - cc2inv*(y2 + z2))) + (u_j[1]*(r2inv*(-3*x*y) + cc2inv*x*y)) + (u_j[2]*(r2inv*(-3*x*z) + cc2inv*x*z)));

                f_local[1] -= common_factor*(u_j[0]*(r2inv*(-3*x*y) + cc2inv*x*y) + u_j[1]*(r2inv*(-2*y2 + x2 + z2) - cc2inv*(x2 + z2)) + u_j[2]*(r2inv*(-3*y*z) + cc2inv*y*z));

                f_local[2] -= common_factor*(u_j[0]*(r2inv*(-3*x*z) + cc2inv*x*z) + u_j[1]*(r2inv*(-3*y*z) + cc2inv*y*z) + u_j[2]*(r2inv*(-2*z2 + x2 + y2) - cc2inv*(x2 + y2)));
            }
            
            // (3) u_i -- u_j  -- assume exponential damping.
            if (system.molecules[i].atoms[j].polar != 0 && system.molecules[k].atoms[l].polar !=0) {
                r5inv = r2inv*r3inv;
                r7inv = r5inv*r2inv;
                udotu = dddotprod(u_i,u_j);
                uidotr = dddotprod(u_i,distances);
                ujdotr = dddotprod(u_j,distances);

                t1 = exp(-damp*r);
                t2 = 1. + damp*r + 0.5*damp*damp*r2;
                t3 = t2 + damp*damp*damp*r3/6.;
                p1 = 3*r5inv*udotu*(1. - t1*t2) - r7inv*15.*uidotr*ujdotr*(1. - t1*t3);
                p2 = 3*r5inv*ujdotr*(1. - t1*t3);
                p3 = 3*r5inv*uidotr*(1. - t1*t3);
                p4 = -udotu*r3inv*(-t1*(damp*rinv + damp*damp) + rinv*t1*damp*t2);
                p5 = 3*r5inv*uidotr*ujdotr*(-t1*(rinv*damp + damp*damp + 0.5*r*damp*damp*damp) + rinv*t1*damp*t3);

                f_local[0] += p1*x + p2*u_i[0] + p3*u_j[0] + p4*x + p5*x;
                f_local[1] += p1*y + p2*u_i[1] + p3*u_j[1] + p4*y + p5*y;
                f_local[2] += p1*z + p2*u_i[2] + p3*u_j[2] + p4*z + p5*z;
            }
            // done with this (i,j)--(k,l) pair. Apply Newton's law
            for (n=0;n<3;n++) {
                system.molecules[i].atoms[j].force[n] += f_local[n];
                system.molecules[k].atoms[l].force[n] -= f_local[n];
            }
        } // end loop l atoms in k
        } // end loop k molecules
    } // end loop j atoms in i
    } // end loop i molecules
} // end polarization forces function
