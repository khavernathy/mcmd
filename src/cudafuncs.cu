#include <stdio.h>
#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>


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


// minimal data to send to GPU. this is all that's needed to calc forces.
typedef struct atom_t {
    double pos[3]={0,0,0};
    double eps=0; // lj
    double sig=0; // lj
    double charge=0;
    double f[3]={0,0,0}; // force
    int molid=0;
    int frozen=0;
} cuda_atom;

/*
// same but for molecule
typedef struct molecule_t {
    double old_ang_acc[3]={0,0,0};
    double ang_acc[3]={0,0,0};
    double ang_vel[3]={0,0,0};
    double ang_pos[3]={0,0,0};
    double torque[3]={0,0,0};
    double inertia=0;
    double mass=0;
    double old_acc[3]={0,0,0};
    double acc[3]={0,0,0};
    double vel[3]={0,0,0};
    double com[3]={0,0,0};
    double force[3]={0,0,0};
} cuda_molecule;
*/

__global__
void calculateForceKernel(cuda_atom * atom_list, int N, double cutoffD, double * basis, double * reciprocal_basis, int pformD, double ewald_alpha, int ewald_num_k, double kmax) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    // only run for real atoms (no ghost threads)
    if(i<N){   
        const register cuda_atom anchoratom = atom_list[i];
        //printf("I AM THREAD %i\n", i);
        //atom_list[i].pos[0] += cutoff;
       const int pform = pformD;
         const double alpha = ewald_alpha;
        const double cutoff=cutoffD;
        register double rimg, rsq;
        const double sqrtPI=sqrt(M_PI);
        double d[3], di[3], img[3], dimg[3],r,r2,ri,ri2;
        int q,j,n;
        double sig,eps,r6,s6,u[3]={0,0,0};
        //int count=0;
        register double af[3] = {0,0,0}; // accumulated forces for anchoratom
        double holder,chargeprod; // for ES force    
        //printf("basis[3] = %f\n", basis[3]);
        __syncthreads();
        // order N instead of N^2 bc this runs on all GPU cores at once (basically)

        // if LJ 
        if (pform == 0 || pform == 1) {
        for (j=i+1;j<N;j++) {

           if (anchoratom.molid == atom_list[j].molid) continue; // skip same molecule 
            if (anchoratom.frozen && atom_list[j].frozen) continue; // skip frozens            

           
            // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            for (n=0;n<3;n++) {
                img[n]=0;
                for (q=0;q<3;q++) {
                    img[n] += reciprocal_basis[n*3+q]*d[q];
                    //if (i==0 && j==1188) printf("img[%i] = reciprocal_basis[%i]*d[%i] = %f\n",p,p*3+q,q,reciprocal_basis[p*3+q]*d[q]);
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
               
            //if (i==0) {
              //  printf("r[%i].%i = %f\n", i,j,rimg);
                //printf("CUTOFF: %f\n", cutoff);
                //for (int h=0;h<9;h++) {
                  //  printf("basis[%i] = %f\n", h, basis[h]);
                //}
            //}

                if (rimg <= cutoff) {
           
                 sig = anchoratom.sig;
                if (sig != atom_list[j].sig) sig = 0.5*(sig+atom_list[j].sig);
                eps = anchoratom.eps;
                if (eps != atom_list[j].eps) eps = sqrt(eps * atom_list[j].eps);

                if (sig == 0 || eps == 0) continue;

     
                
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
        if (pform == 1) {
            for (n=0;n<3;n++) af[n]=0; // reset register-stored force for anchoratom.
          
                // get inverse volume
                double invV = basis[0]*(basis[4]*basis[8] - basis[7]*basis[5]);
                invV +=       basis[3]*(basis[7]*basis[2] - basis[1]*basis[8]);
                invV +=       basis[6]*(basis[1]*basis[5] - basis[5]*basis[2]);
                invV = 1.0/invV;
 
                double k_sq; double fourPI = M_PI*4;

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


            // real-space forces
            if (rimg <= cutoff && (anchoratom.molid < atom_list[j].molid)) { // non-duplicated pairs, not intramolecular, not beyond cutoff
                chargeprod = anchoratom.charge * atom_list[j].charge;
                for (n=0;n<3;n++) u[n] = dimg[n]/rimg;
                for (n=0;n<3;n++) {
                    holder = -((-2.0*chargeprod*alpha*exp(-alpha*alpha*rsq))/(sqrtPI*rimg) - (chargeprod*erfc(alpha*rimg)/rsq))*u[n];
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);                
                }
            }
            // k-space forces
			if (anchoratom.molid < atom_list[j].molid) { // no cutoff for this.
            for (int n=0;n<3;n++) { //x,y,z
              // loop k vectors
				// EWALD k-vectors     
					int l[3],p,q; double k[3] = {0,0,0};
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
                    holder = chargeprod*invV*fourPI*k[n]*
                        exp(-k_sq/(4*alpha*alpha))*
                        sin(k[0]*dimg[0]+
                            k[1]*dimg[1]+
                            k[2]*dimg[2])/k_sq;
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);
								
								
							} // end for l[2], n
						} // end for l[1], m
					} // end for l[0], l
	            } // end 3D
            } // end k-space if not duplicate

            } // end pair loop j for anchoratom

            // finally add ES contribution to anchor-atom
            for (n=0;n<3;n++) atomicAdd(&(atom_list[i].f[n]), af[n]);
        } // end ES component

        //if (i==0) printf("COUNT: %i\n",count);
    } // end if i<n (all threads)
}


__global__
void calculateForceNopbcKernel(cuda_atom * atom_list, int N, int pformD) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    // only run for real atoms (no ghost threads)
    if(i<N){   
        const register cuda_atom anchoratom = atom_list[i];
        //printf("I AM THREAD %i\n", i);
        //atom_list[i].pos[0] += cutoff;
       const int pform = pformD;
        const double cutoff=10.; // default 10 A for no-pbc cutoff.
        double d[3], r,r2;
        int j,n;
        double sig,eps,r6,s6,u[3]={0,0,0};
        //int count=0;
        register double af[3] = {0,0,0}; // accumulated forces for anchoratom
        double holder,chargeprod; // for ES force    
        //printf("basis[3] = %f\n", basis[3]);
        __syncthreads();
        // order N instead of N^2 bc this runs on all GPU cores at once (basically)

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


/*
__global__
void velocityVerletKernel(cuda_molecule * molecule_list, int N, int md_mode) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    register cuda_molecule anchormolecule = molecule_list[i];

    // I ain't afraid o' no ghost 
    if(i<N){
            



    }
}



void CUDA_verlet(System &system) {
    const int N = (int)system.stats.count_movables;
    const int block_size = system.constants.cuda_block_size;
    const int molecules_array_size=sizeof(cuda_molecule)*N;
    cuda_atom H[N]; // host atoms
    cuda_atom *D; // device atoms (gpu)
    for (int i=0;i<system.molecules.size();i++) {
        for (int n=0;n<3;n++) {
            H[i].old_acc[n] = system.molecules[i].old_acc[n];
            H[i].acc[n] = system.molecules[i].acc[n];
            H[i].old_ang_acc[n] = system.molecules[i].old_ang_acc[n];
            H[i].ang_acc[n] = system.molecules[i].ang_acc[n];
            H[i].vel[n] = system.molecules[i].vel[n];
            H[i].ang_vel[n] = system.molecules[i].ang_vel[n];
            
        }
    }

    // allocate memory on GPU
    cudaMalloc((void**) &D, molecules_array_size);
    cudaMemcpy(D, H, molecules_array_size, cudaMemcpyHostToDevice);

    // grid elements
    int dimGrid = ceil((double)N/block_size);
    int dimBlock = block_size;

    // determine molecular or atomic motion
    int md_mode = system.constants.md_mode;
    if (md_mode == MD_ATOMIC) md_mode = 0;
    else if (md_mode == MD_MOLECULAR) md_mode = 1;

    velocityVerletKernel<<< dimGrid, dimBlock >>>(D,N,md_mode);
    // make sure the threads are synced so we don't overflow
    cudaThreadSynchronize();
    // copy device data back to host
    cudaMemcpy(H, D, molecules_array_size, cudaMemcpyDeviceToHost);

    index=0;
    for (int i=0;i<system.molecules.size();i++) {
        for (int j=0;j<system.molecules[i].atoms.size();j++) {
            for (int n=0;n<3;n++) {
                system.molecules[i].atoms[j].force[n] = H[index].f[n];
            }
            index++;
        }
    }

    //printf("H[0] force = %f %f %f\n",system.molecules[0].atoms[0].force[0], system.molecules[0].atoms[0].force[1], system.molecules[0].atoms[0].force[2]);


     cudaFree(D);



}
*/

void CUDA_force(System &system) {

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

    int bs = sizeof(double)*9;
    double *basis;
    double *reciprocal_basis;
    basis = (double*)malloc(bs);
    reciprocal_basis = (double*)malloc(bs);
    double *dbasis;
    double *dreciprocal_basis;

    for (int p=0;p<3;p++) {
        for (int q=0;q<3;q++) {
            basis[3*q+p] = system.pbc.basis[p][q]; // quite sure correct.
            reciprocal_basis[3*q+p] = system.pbc.reciprocal_basis[p][q]; // quite sure correct.
        }
    }
    //system.pbc.printBasis();

    

    //for (int l=0;l<9;l++) printf("basis[%i] = %f\n", l,basis[l]);

    // allocate memory on GPU
    cudaMalloc((void**) &dbasis, bs);
    cudaMemcpy(dbasis, basis, bs, cudaMemcpyHostToDevice);
    cudaMalloc((void**) &dreciprocal_basis, bs);
    cudaMemcpy(dreciprocal_basis, reciprocal_basis, bs, cudaMemcpyHostToDevice); 
    cudaMalloc((void**) &D, atoms_array_size);
    cudaMemcpy(D, H, atoms_array_size, cudaMemcpyHostToDevice);
	//cudaMalloc((void**) &Dkvecs, kvecs_size);
	//cudaMemcpy(Dkvecs, Hkvecs, kvecs_size, cudaMemcpyHostToDevice);

    // grid elements
    int dimGrid = ceil((double)N/block_size);
    int dimBlock = block_size;   

    // assign potential form for force calculator
    int pform,theval=system.constants.potential_form;
    if (theval == POTENTIAL_LJ || theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=0;
    if (theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=1;

    calculateForceKernel<<< dimGrid, dimBlock >>>(D,N,system.pbc.cutoff, dbasis, dreciprocal_basis, pform, system.constants.ewald_alpha, system.constants.ewald_num_k, system.constants.ewald_kmax);
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

    //printf("H[0] force = %f %f %f\n",system.molecules[0].atoms[0].force[0], system.molecules[0].atoms[0].force[1], system.molecules[0].atoms[0].force[2]);
 

    // clean up -- so we don't have a memory leak
     cudaFree(D);
     cudaFree(dbasis);
     cudaFree(dreciprocal_basis);
     free(basis);
     free(reciprocal_basis);


    // we're done. forces have been calc'd on GPU and written to local mem.
}

void CUDA_force_nopbc(System &system) {

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

    // we're done. forces have been calc'd on GPU and written to local mem.
}
