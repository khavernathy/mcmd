#include <stdio.h>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#if __CUDA_ARCH__ < 200
/* this is an explicit definition for atomicAdd, to be safe */
__device__ double atomicAdd(double* address, double val)
{
 unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do { assumed = old;
  old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
  }
  while (assumed != old);
  return __longlong_as_double(old);
}
#endif

// minimal data to send to GPU. this is all that's needed to calc forces.
typedef struct atom_t {
    double pos[3]={0,0,0};
    double eps=0; // lj
    double sig=0; // lj
    double charge=0;
    double f[3]={0,0,0}; // force
    int molid=0;
    int frozen=0;
    double u[3] = {0,0,0}; // dipole
    double polar=0; // polarizability
} cuda_atom;


__global__
void calculateForceKernel(cuda_atom * atom_list, int N, double cutoffD, double * basis, double * reciprocal_basis, int pformD, double ewald_alpha, int kmax, int kspace, double polar_damp) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    // only run for real atoms (no ghost threads)
    if(i<N){   
        const register cuda_atom anchoratom = atom_list[i];
        const int pform = pformD;
        const double damp = polar_damp;
        const double alpha = ewald_alpha;
        const double cutoff=cutoffD;
        register double rimg, rsq;
        const double sqrtPI=sqrt(M_PI);
        double d[3], di[3], img[3], dimg[3],r,r2,ri,ri2;
        int q,j,n;
        double sig,eps,r6,s6,u[3]={0,0,0};
        register double af[3] = {0,0,0}; // accumulated forces for anchoratom
        double holder,chargeprod; // for ES force    
        __syncthreads();
        // order N (generally, IF N_atoms < N_cores) instead of N^2 bc 
        // this runs on all GPU cores at once 

        // if LJ 
        if (pform == 0 || pform == 1 || pform == 2) {
        for (j=i+1;j<N;j++) {

           if (anchoratom.molid == atom_list[j].molid) continue; // skip same molecule 
            if (anchoratom.frozen && atom_list[j].frozen) continue; // skip frozens            
            
            // LB mixing
            sig = anchoratom.sig;
            if (sig != atom_list[j].sig) sig = 0.5*(sig+atom_list[j].sig);
            eps = anchoratom.eps;
            if (eps != atom_list[j].eps) eps = sqrt(eps * atom_list[j].eps);
            if (sig == 0 || eps == 0) continue;

            // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            for (n=0;n<3;n++) {
                img[n]=0;
                for (q=0;q<3;q++) {
                    img[n] += reciprocal_basis[n*3+q]*d[q];
                }
                img[n] = rint(img[n]);
            }
            for (n=0;n<3;n++) {
                di[n] = 0;
                for (q=0;q<3;q++) {
                    di[n] += basis[n*3+q]*img[q];
                }
            }
            for (n=0;n<3;n++) di[n] = d[n] - di[n];
            r2=0;ri2=0;
            for (n=0;n<3;n++) {
                r2 += d[n]*d[n];
                ri2 += di[n]*di[n];
            }
            r = sqrt(r2);
            ri = sqrt(ri2);
            if (ri != ri) {
                rimg=r;
                rsq=r2;
                for (n=0;n<3;n++) dimg[n] = d[n];
            } else {
                rimg=ri;
                rsq=ri2;
                for (n=0;n<3;n++) dimg[n] = di[n];
            }
            // distance is now rimg
    
            if (rimg <= cutoff) {
                r6 = rsq*rsq*rsq;
                s6 = sig*sig;
                s6 *= s6 * s6;
        
                for (n=0;n<3;n++) {
                    holder = 24.0*dimg[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                    atomicAdd(&(atom_list[j].f[n]), -holder); 
                    af[n] += holder;      
                }
            }
        } // end pair j
        
        // finally add the accumulated forces (stored on register) to the anchor atom
        for (n=0;n<3;n++)
            atomicAdd(&(atom_list[i].f[n]), af[n]);
        
        } // end if LJ
        // ==============================================================================
        // Now handle electrostatics
        if (pform == 1 || pform == 2) {
            for (n=0;n<3;n++) af[n]=0; // reset register-stored force for anchoratom.
            double invV;
            int l[3], p, q; double k[3], k_sq, fourPI = 4.0*M_PI;
            invV =  basis[0] * (basis[4]*basis[8] - basis[7]*basis[5] );
            invV += basis[3] * (basis[7]*basis[2] - basis[1]*basis[8] );
            invV += basis[6] * (basis[1]*basis[5] - basis[5]*basis[2] );
            invV = 1.0/invV;

           for (j=0;j<N;j++) {
                if (anchoratom.frozen && atom_list[j].frozen) continue; // don't do frozen pairs
                if (anchoratom.charge == 0 || atom_list[j].charge == 0) continue; // skip 0-force
                if (i==j) continue; // don't do atom with itself

               // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            for (n=0;n<3;n++) {
                img[n]=0;
                for (q=0;q<3;q++) {
                    img[n] += reciprocal_basis[n*3+q]*d[q];
                }
                img[n] = rint(img[n]);
            }
            for (n=0;n<3;n++) {
                di[n] = 0;
                for (q=0;q<3;q++) {
                    di[n] += basis[n*3+q]*img[q];
                }
            }
            for (n=0;n<3;n++) di[n] = d[n] - di[n];
            r2=0;ri2=0;
            for (n=0;n<3;n++) {
                r2 += d[n]*d[n];
                ri2 += di[n]*di[n];
            }
            r = sqrt(r2);
            ri = sqrt(ri2);
            if (ri != ri) {
                rimg=r;
                rsq=r2;
                for (n=0;n<3;n++) dimg[n] = d[n];
            } else {
                rimg=ri;
                rsq=ri2;
                for (n=0;n<3;n++) dimg[n] = di[n];
            }

            // real-space
            if (rimg <= cutoff && (anchoratom.molid < atom_list[j].molid)) { // non-duplicated pairs, not intramolecular, not beyond cutoff
                chargeprod = anchoratom.charge * atom_list[j].charge;
                for (n=0;n<3;n++) u[n] = dimg[n]/rimg;
                for (n=0;n<3;n++) {
                    holder = -((-2.0*chargeprod*alpha*exp(-alpha*alpha*rsq))/(sqrtPI*rimg) - (chargeprod*erfc(alpha*rimg)/rsq))*u[n];
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);                
                }
            }
            // k-space
            if (kspace && (anchoratom.molid < atom_list[j].molid)) {
               chargeprod = anchoratom.charge * atom_list[j].charge;
                
               for (n=0; n<3; n++) {
                    for (l[0] = 0; l[0] <= kmax; l[0]++) {
                    for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
                    for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {
                        // skip if norm is out of sphere
                        if (l[0]*l[0] + l[1]*l[1] + l[2]*l[2] > kmax*kmax) continue;
                        /* get reciprocal lattice vectors */
                        for (p=0; p<3; p++) {
                            for (q=0, k[p] = 0; q < 3; q++) {
                            k[p] += 2.0*M_PI*reciprocal_basis[3*q+p] * l[q];
                            }
                        }
                        k_sq = k[0]*k[0] + k[1]*k[1] + k[2]*k[2];

                        holder = chargeprod * invV * fourPI * k[n] *
                        exp(-k_sq/(4*alpha*alpha))*
                        sin(k[0]*dimg[0] + k[1]*dimg[1] + k[2]*dimg[2])/k_sq * 2; // times 2 b/c half-Ewald sphere

                        af[n] += holder;
                        atomicAdd(&(atom_list[j].f[n]), -holder);

                    } // end for l[2], n
                    } // end for l[1], m
                    } // end for l[0], l
                } // end 3d 
            }

            } // end pair loop j 

            // finally add ES contribution to anchor-atom
            for (n=0;n<3;n++) atomicAdd(&(atom_list[i].f[n]), af[n]);
        } // end ES component
        // ============================================================
        // Polarization
        if (pform == 2) {
            double common_factor, r, rinv, r2, r2inv, r3, r3inv, r5inv, r7inv;
            double x2,y2,z2,x,y,z;
            double udotu, ujdotr, uidotr;
            const double cc2inv = 1.0/(cutoff*cutoff);
            double t1,t2,t3,p1,p2,p3,p4,p5;
            const double u_i[3] = {anchoratom.u[0], anchoratom.u[1], anchoratom.u[2]};
            double u_j[3];
            // loop all pair atoms
            for (int j=i+1; j<N; j++) {
                for (n=0;n<3;n++) af[n] = 0; // reset local force for this pair.
                if (anchoratom.molid == atom_list[j].molid) continue; // no same-molecule
                // get R (nearest image)
		        for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
		        for (n=0;n<3;n++) {
		            img[n]=0;
		            for (q=0;q<3;q++) {
		                img[n] += reciprocal_basis[n*3+q]*d[q];
		            }
		            img[n] = rint(img[n]);
		        }
		        for (n=0;n<3;n++) {
		            di[n] = 0;
		            for (q=0;q<3;q++) {
		                di[n] += basis[n*3+q]*img[q];
		            }
		        }
		        for (n=0;n<3;n++) di[n] = d[n] - di[n];
		        r2=0;ri2=0;
		        for (n=0;n<3;n++) {
		            r2 += d[n]*d[n];
		            ri2 += di[n]*di[n];
		        }
		        r = sqrt(r2);
		        ri = sqrt(ri2);
		        if (ri != ri) {
		            rimg=r;
		            rsq=r2;
		            for (n=0;n<3;n++) dimg[n] = d[n];
		        } else {
		            rimg=ri;
		            rsq=ri2;
		            for (n=0;n<3;n++) dimg[n] = di[n];
		        }
                // got pair displacements

                if (rimg > cutoff) continue; // skip outside cutoff
                r = rimg;
                x = dimg[0]; y = dimg[1]; z = dimg[2];
                x2 = x*x;
                y2 = y*y;
                z2 = z*z;
                r2 = r*r;
                r3 = r2*r;
                rinv = 1./r;
                r2inv = rinv*rinv;
                r3inv = r2inv*rinv;
                for (n=0;n<3;n++) u_j[n] = atom_list[j].u[n];

                // (1) u_i -- q_j
                if (atom_list[j].charge != 0 && anchoratom.polar != 0) {
                    common_factor = atom_list[j].charge * r3inv;

		            af[0] += common_factor*((u_i[0]*(r2inv*(-2*x2 + y2 + z2) - cc2inv*(y2 + z2))) + (u_i[1]*(r2inv*(-3*x*y) + cc2inv*x*y)) + (u_i[2]*(r2inv*(-3*x*z) + cc2inv*x*z)));

		            af[1] += common_factor*(u_i[0]*(r2inv*(-3*x*y) + cc2inv*x*y) + u_i[1]*(r2inv*(-2*y2 + x2 + z2) - cc2inv*(x2 + z2)) + u_i[2]*(r2inv*(-3*y*z) + cc2inv*y*z));

		            af[2] += common_factor*(u_i[0]*(r2inv*(-3*x*z) + cc2inv*x*z) + u_i[1]*(r2inv*(-3*y*z) + cc2inv*y*z) + u_i[2]*(r2inv*(-2*z2 + x2 + y2) - cc2inv*(x2 + y2)));
 
                }

                // (2) u_j -- q_i
                if (anchoratom.charge != 0 && atom_list[j].polar != 0) {
                    common_factor = anchoratom.charge * r3inv;

                    af[0] -= common_factor*((u_j[0]*(r2inv*(-2*x2 + y2 + z2) - cc2inv*(y2 + z2))) + (u_j[1]*(r2inv*(-3*x*y) + cc2inv*x*y)) + (u_j[2]*(r2inv*(-3*x*z) + cc2inv*x*z)));

                    af[1] -= common_factor*(u_j[0]*(r2inv*(-3*x*y) + cc2inv*x*y) + u_j[1]*(r2inv*(-2*y2 + x2 + z2) - cc2inv*(x2 + z2)) + u_j[2]*(r2inv*(-3*y*z) + cc2inv*y*z));

                    af[2] -= common_factor*(u_j[0]*(r2inv*(-3*x*z) + cc2inv*x*z) + u_j[1]*(r2inv*(-3*y*z) + cc2inv*y*z) + u_j[2]*(r2inv*(-2*z2 + x2 + y2) - cc2inv*(x2 + y2)));
                }

                // (3) u_i -- u_j
                if (anchoratom.polar != 0 && atom_list[j].polar != 0) {
                    r5inv = r2inv*r3inv;
                    r7inv = r5inv*r2inv;
                    udotu = u_i[0]*u_j[0] + u_i[1]*u_j[1] + u_i[2]*u_j[2];
                    uidotr = u_i[0]*dimg[0] + u_i[1]*dimg[1] + u_i[2]*dimg[2];
                    ujdotr = u_j[0]*dimg[0] + u_j[1]*dimg[1] + u_j[2]*dimg[2];                    

                    t1 = exp(-damp*r);
                    t2 = 1. + damp*r + 0.5*damp*damp*r2;
                    t3 = t2 + damp*damp*damp*r3/6.;
                    p1 = 3*r5inv*udotu*(1. - t1*t2) - r7inv*15.*uidotr*ujdotr*(1. - t1*t3);
                    p2 = 3*r5inv*ujdotr*(1. - t1*t3);
                    p3 = 3*r5inv*uidotr*(1. - t1*t3);
                    p4 = -udotu*r3inv*(-t1*(damp*rinv + damp*damp) + rinv*t1*damp*t2);
                    p5 = 3*r5inv*uidotr*ujdotr*(-t1*(rinv*damp + damp*damp + 0.5*r*damp*damp*damp) + rinv*t1*damp*t3);

                    af[0] += p1*x + p2*u_i[0] + p3*u_j[0] + p4*x + p5*x;
                    af[1] += p1*y + p2*u_i[1] + p3*u_j[1] + p4*y + p5*y;
                    af[2] += p1*z + p2*u_i[2] + p3*u_j[2] + p4*z + p5*z;
                }

                // apply Newton for pair.
                for (n=0;n<3;n++) {
                    atomicAdd(&(atom_list[i].f[n]), af[n]);
                    atomicAdd(&(atom_list[j].f[n]), -af[n]);    
                }

		    } // end pair loop with atoms j
        } // end polarization forces

    } // end if i<n (all threads)
}


__global__
void calculateForceNopbcKernel(cuda_atom * atom_list, int N, int pformD) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    // only run for real atoms (no ghost threads)
    if(i<N){   
        const register cuda_atom anchoratom = atom_list[i];
       const int pform = pformD;
        const double cutoff=10.; // default 10 A for no-pbc cutoff.
        double d[3], r,r2;
        int j,n;
        double sig,eps,r6,s6,u[3]={0,0,0};
        register double af[3] = {0,0,0}; // accumulated forces for anchoratom
        double holder,chargeprod; // for ES force    
        __syncthreads();

        // if LJ 
        if (pform == 0 || pform == 1) {
        for (j=i+1;j<N;j++) {

           if (anchoratom.molid == atom_list[j].molid) continue; // skip same molecule 
            if (anchoratom.frozen && atom_list[j].frozen) continue; // skip frozens            

           
            // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            r2=0;
            for (n=0;n<3;n++) {
                r2 += d[n]*d[n];
            }
            r = sqrt(r2);
               
                if (r <= cutoff) {
           
                 sig = anchoratom.sig;
                if (sig != atom_list[j].sig) sig = 0.5*(sig+atom_list[j].sig);
                eps = anchoratom.eps;
                if (eps != atom_list[j].eps) eps = sqrt(eps * atom_list[j].eps);

                if (sig == 0 || eps == 0) continue;
                
                r6 = r2*r2*r2;
                s6 = sig*sig;
                s6 *= s6 * s6;
        
                    for (n=0;n<3;n++) {
                        holder = 24.0*d[n]*eps*(2*(s6*s6)/(r6*r6*r2) - s6/(r6*r2));
                        atomicAdd(&(atom_list[j].f[n]), -holder); 
                        af[n] += holder;      
                    }
                }

        } // end pair j
        
        // finally add the accumulated forces (stored on register) to the anchor atom
        for (n=0;n<3;n++)
            atomicAdd(&(atom_list[i].f[n]), af[n]);
        
        } // end if LJ
        // ==============================================================================
        // Now handle electrostatics
        if (pform == 1) {
            for (n=0;n<3;n++) af[n]=0; // reset register-stored force for anchoratom.
           for (j=i+1;j<N;j++) {
                if (anchoratom.frozen && atom_list[j].frozen) continue; // don't do frozen pairs
                if (anchoratom.charge == 0 || atom_list[j].charge == 0) continue; // skip 0-force
                if (anchoratom.molid == atom_list[j].molid) continue; // don't do molecule with itself

               // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            r2=0;
            for (n=0;n<3;n++) {
                r2 += d[n]*d[n];
            }
            r = sqrt(r2);

            if (r <= cutoff)  { //&& (anchoratom.molid < atom_list[j].molid)) { // non-duplicated pairs, not intramolecular, not beyond cutoff
                chargeprod = anchoratom.charge * atom_list[j].charge;
                for (n=0;n<3;n++) u[n] = d[n]/r;
                for (n=0;n<3;n++) {
                    holder = chargeprod/r2 * u[n];
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);                
                }
            }

            } // end pair loop j 

            // finally add ES contribution to anchor-atom
            for (n=0;n<3;n++) atomicAdd(&(atom_list[i].f[n]), af[n]);
        } // end ES component

    } // end if i<n (all threads)
} // end no-pbc force


void CUDA_force(System &system) {

    const int N = (int)system.constants.total_atoms;
    const int block_size = system.constants.cuda_block_size; 
    const int atoms_array_size=sizeof(cuda_atom)*N;
    int index=0;

    // if polarization force needed, get dipoles on CPU first
    if (system.constants.potential_form == POTENTIAL_LJESPOLAR) {
        if (system.constants.ensemble == ENSEMBLE_UVT) {
            thole_resize_matrices(system); // only if N can change
        }
        thole_amatrix(system); // populate A matrix
        thole_field(system); // calculate electric field
        int num_iterations = thole_iterative(system); // calculate dipoles
            system.stats.polar_iterations.value = (double)num_iterations;
            system.stats.polar_iterations.calcNewStats();
            system.constants.dipole_rrms = get_dipole_rrms(system);
    }

    cuda_atom H[N]; // host atoms
    cuda_atom *D; // device atoms (gpu)
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            H[index].molid = i;
            H[index].sig = system.molecules[i].atoms[j].sig;
            H[index].eps = system.molecules[i].atoms[j].eps;
            H[index].charge = system.molecules[i].atoms[j].C;
            if (system.constants.potential_form == POTENTIAL_LJESPOLAR)
                H[index].polar = system.molecules[i].atoms[j].polar;
            for (int n=0;n<3;n++) {
                H[index].pos[n] = system.molecules[i].atoms[j].pos[n];       
                H[index].f[n] = 0; // initialize to zero
                if (system.constants.potential_form == POTENTIAL_LJESPOLAR) {
                    H[index].u[n] = system.molecules[i].atoms[j].dip[n];
                }
            }
            H[index].frozen = system.molecules[i].atoms[j].frozen;     
            index++;       
        }
    }

    int bs = sizeof(double)*9;
    double *basis;
    double *reciprocal_basis;
    basis = (double*)malloc(bs);
    reciprocal_basis = (double*)malloc(bs);
    double *dbasis;
    double *dreciprocal_basis;

    for (int p=0;p<3;p++) {
        for (int q=0;q<3;q++) {
            basis[3*q+p] = system.pbc.basis[p][q]; 
            reciprocal_basis[3*q+p] = system.pbc.reciprocal_basis[p][q]; 
        }
    }

    // allocate memory on GPU
    cudaMalloc((void**) &dbasis, bs);
    cudaMemcpy(dbasis, basis, bs, cudaMemcpyHostToDevice);
    cudaMalloc((void**) &dreciprocal_basis, bs);
    cudaMemcpy(dreciprocal_basis, reciprocal_basis, bs, cudaMemcpyHostToDevice); 
    cudaMalloc((void**) &D, atoms_array_size);
    cudaMemcpy(D, H, atoms_array_size, cudaMemcpyHostToDevice);

    // grid elements
    int dimGrid = ceil((double)N/block_size);
    int dimBlock = block_size;   

    // assign potential form for force calculator
    int pform,theval=system.constants.potential_form;
    if (theval == POTENTIAL_LJ || theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=0;
    if (theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=1;
    if (theval == POTENTIAL_LJESPOLAR)
        pform=2;

    // calculate forces on GPU !!
    calculateForceKernel<<< dimGrid, dimBlock >>>(D,N,system.pbc.cutoff, dbasis, dreciprocal_basis, pform, system.constants.ewald_alpha, system.constants.ewald_kmax, system.constants.kspace_option, system.constants.polar_damp);
    // make sure the threads are synced so we don't overflow
    cudaThreadSynchronize();
    // copy device data back to host
    cudaMemcpy(H, D, atoms_array_size, cudaMemcpyDeviceToHost);

    index=0;
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            for (int n=0;n<3;n++) {
                system.molecules[i].atoms[j].force[n] = H[index].f[n];
            }     
            index++;       
        }
    }

    // clean up -- so we don't have a memory leak
     cudaFree(D);
     cudaFree(dbasis);
     cudaFree(dreciprocal_basis);
     free(basis);
     free(reciprocal_basis);

}

void CUDA_force_nopbc(System &system) {

/*
    needs works to be fully functional.
*/


    const int N = (int)system.constants.total_atoms;
    const int block_size = system.constants.cuda_block_size; 
    const int atoms_array_size=sizeof(cuda_atom)*N;
    int index=0;

    cuda_atom H[N]; // host atoms
    cuda_atom *D; // device atoms (gpu)
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            H[index].molid = i;
            H[index].sig = system.molecules[i].atoms[j].sig;
            H[index].eps = system.molecules[i].atoms[j].eps;
            H[index].charge = system.molecules[i].atoms[j].C;
            for (int n=0;n<3;n++) {
                H[index].pos[n] = system.molecules[i].atoms[j].pos[n];       
                H[index].f[n] = 0; // initialize to zero
            }
            H[index].frozen = system.molecules[i].atoms[j].frozen;     
            index++;       
        }
    }

    // allocate memory on GPU
    cudaMalloc((void**) &D, atoms_array_size);
    cudaMemcpy(D, H, atoms_array_size, cudaMemcpyHostToDevice);

    // grid elements
    int dimGrid = ceil((double)N/block_size);
    int dimBlock = block_size;   

    // assign potential form for force calculator
    int pform,theval=system.constants.potential_form;
    if (theval == POTENTIAL_LJ || theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=0;
    if (theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=1;

    calculateForceNopbcKernel<<< dimGrid, dimBlock >>>(D,N, pform);
    // make sure the threads are synced so we don't overflow
    cudaThreadSynchronize();
    // copy device data back to host
    cudaMemcpy(H, D, atoms_array_size, cudaMemcpyDeviceToHost);

    //for (int i=0;i<N;i++) printf("H[%i] force0 = %f\n", i, H[i].f[0]);
    index=0;
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            for (int n=0;n<3;n++) {
                system.molecules[i].atoms[j].force[n] = H[index].f[n];
            }     
            index++;       
        }
    }

     cudaFree(D);
}
