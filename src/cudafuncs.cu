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
void calculateForceKernel(cuda_atom * atom_list, int N, double cutoff, double * basis, double * reciprocal_basis, int pform, double ewald_alpha) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    register cuda_atom anchoratom = atom_list[i];

    // only run for real atoms (no ghost threads)
    if(i<N){   
        //printf("I AM THREAD %i\n", i);
        //atom_list[i].pos[0] += cutoff;
        const double alpha = ewald_alpha;
        register double rimg, rsq;
        const double sqrtPI=sqrt(M_PI);
        double d[3], di[3], img[3], dimg[3],r,r2,ri,ri2;
        int p,q,j,n;
        double sig,eps,r6,s6,f[3]={0,0,0},u[3]={0,0,0};
        //int count=0;
        register double af[3] = {0,0,0}; // accumulated forces for anchoratom
        double holder,erfc_term,chargeprod; // for ES force    
        //printf("basis[3] = %f\n", basis[3]);
        __syncthreads();
        // order N instead of N^2 bc this runs on all GPU cores at once (basically)

        // if LJ 
        if (pform == 0 || pform == 1) {
        for (j=i+1;j<N;j++) {

           if (anchoratom.molid == atom_list[j].molid) continue; // skip same molecule 
            // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            for (p=0;p<3;p++) {
                img[p]=0;
                for (q=0;q<3;q++) {
                    img[p] += reciprocal_basis[p*3+q]*d[q];
                    //if (i==0 && j==1188) printf("img[%i] = reciprocal_basis[%i]*d[%i] = %f\n",p,p*3+q,q,reciprocal_basis[p*3+q]*d[q]);
                }
                img[p] = rint(img[p]);
            }
            for (p=0;p<3;p++) {
                di[p] = 0;
                for (q=0;q<3;q++) {
                    di[p] += basis[p*3+q]*img[q];
                }
            }
            for (p=0;p<3;p++) di[p] = d[p] - di[p];
            r2=0;ri2=0;
            for (p=0;p<3;p++) {
                r2 += d[p]*d[p];
                ri2 += di[p]*di[p];
            }
            r = sqrt(r2);
            ri = sqrt(ri2);
            if (ri != ri) {
                rimg=r;
                for (p=0;p<3;p++) dimg[p] = d[p];
            } else {
                rimg=ri;
                for (p=0;p<3;p++) dimg[p] = di[p];
            }
            // distance is now rimg
               
            //if (i==0) {
              //  printf("r[%i].%i = %f\n", i,j,rimg);
                //printf("CUTOFF: %f\n", cutoff);
                //for (int h=0;h<9;h++) {
                  //  printf("basis[%i] = %f\n", h, basis[h]);
                //}
            //}
            
            rsq=rimg*rimg;

                //if (i==0) printf("hi\n");

                sig = anchoratom.sig;
                if (sig != atom_list[j].sig) sig = 0.5*(sig+atom_list[j].sig);
                eps = anchoratom.eps;
                if (eps != atom_list[j].eps) eps = sqrt(anchoratom.eps * atom_list[j].eps);

                if (sig == 0 || eps == 0) continue;

                r6 = rsq*rsq*rsq;
                s6 = sig*sig;
                s6 *= s6 * s6;

                if (rimg <= cutoff) {
                    for (n=0;n<3;n++) {
                        f[n] = 24.0*dimg[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                        atomicAdd(&(atom_list[j].f[n]), -f[n]); 
                        af[n] += f[n];      
                        
                        //af[n] += f[n];
                    }
                    //if (i==0) count++;
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
           for (j=0;j<N;j++) {
                if (i==j) continue; // don't do atom with itself
                if (anchoratom.frozen && atom_list[j].frozen) continue; // don't do frozen pairs
                if (anchoratom.charge == 0 && atom_list[j].charge == 0) continue; // skip 0-force

                chargeprod = anchoratom.charge * atom_list[j].charge;
            

               // get R (nearest image)
            for (n=0;n<3;n++) d[n] = anchoratom.pos[n] - atom_list[j].pos[n];
            for (p=0;p<3;p++) {
                img[p]=0;
                for (q=0;q<3;q++) {
                    img[p] += reciprocal_basis[p*3+q]*d[q];
                }
                img[p] = rint(img[p]);
            }
            for (p=0;p<3;p++) {
                di[p] = 0;
                for (q=0;q<3;q++) {
                    di[p] += basis[p*3+q]*img[q];
                }
            }
            for (p=0;p<3;p++) di[p] = d[p] - di[p];
            r2=0;ri2=0;
            for (p=0;p<3;p++) {
                r2 += d[p]*d[p];
                ri2 += di[p]*di[p];
            }
            r = sqrt(r2);
            ri = sqrt(ri2);
            if (ri != ri) {
                rimg=r;
                for (p=0;p<3;p++) dimg[p] = d[p];
            } else {
                rimg=ri;
                for (p=0;p<3;p++) dimg[p] = di[p];
            }

            rsq=rimg*rimg;
            for (n=0;n<3;n++) u[n] = dimg[n]/rimg;

            if (r <= cutoff && (anchoratom.molid < atom_list[j].molid)) { // non-duplicated pairs, not intramolecular, not beyond cutoff
                erfc_term = erfc(alpha*r);
                for (n=0;n<3;n++) {
                    holder = -((-2.0*chargeprod*alpha*exp(-alpha*alpha*r*r))/(sqrtPI*r) - (chargeprod*erfc_term/rsq))*u[n];
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);                
                }
            } else if (anchoratom.molid == atom_list[j].molid && i != j) { // intramolecular interaction
                for (n=0;n<3;n++) {
                    holder = -((chargeprod*erf(alpha*r))/rsq - (2*chargeprod*alpha*exp(-alpha*alpha*r*r)/(sqrtPI*r)))*u[n];
                    af[n] += holder;
                    atomicAdd(&(atom_list[j].f[n]), -holder);
                }
            }

            } // end pair loop j 

            // finally add ES contribution to anchor-atom
            for (n=0;n<3;n++) atomicAdd(&(atom_list[i].f[n]), af[n]);
        } // end ES component

        //if (i==0) printf("COUNT: %i\n",count);
    } // end if i<n (all threads)
}


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

    // grid elements
    int dimGrid = ceil((double)N/block_size);
    int dimBlock = block_size;   

    // assign potential form for force calculator
    int pform,theval=system.constants.potential_form;
    if (theval == POTENTIAL_LJ || theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=0;
    if (theval == POTENTIAL_LJES || theval == POTENTIAL_LJESPOLAR)
        pform=1;

    calculateForceKernel<<< dimGrid, dimBlock >>>(D,N,system.pbc.cutoff, dbasis, dreciprocal_basis, pform, system.constants.ewald_alpha);
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
 

     cudaFree(D);

    // we're done. forces have been calc'd on GPU and written to local mem.
}
