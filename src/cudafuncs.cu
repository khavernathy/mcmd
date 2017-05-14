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
} cuda_atom;

__global__
void calculateForceKernel(cuda_atom * atom_list, int N, double cutoff, double * basis, double * reciprocal_basis, int pform) {
    // define thread id
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    register cuda_atom anchoratom = atom_list[i];

    // PROBLEM WITH GLOBAL MEMORY ACCESS CONFLICT?


    // only run for real atoms (no ghost threads)
    if(i<N){   
        //printf("I AM THREAD %i\n", i);
        //atom_list[i].pos[0] += cutoff;
        register double rimg, rsq;
        double d[3], di[3], img[3], dimg[3],r,r2,ri,ri2;
        int p,q,j,n;
        double sig,eps,r6,s6,f[3]={0,0,0};//u[3];
        //int count=0;
        //register double af[3] = {0,0,0}; // accumulated forces
            //printf("basis[3] = %f\n", basis[3]);
        __syncthreads();
        // order N instead of N^2 bc this runs on all GPU cores at once (basically)
        for (j=i+1;j<N;j++) {

            if (anchoratom.molid == atom_list[j].molid) continue; // skip self-molecule interactions

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

            // 0 is LJ, 1 is LJ+ES
            if (pform == 0 || pform == 1) {
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
                        atomicAdd(&(atom_list[i].f[n]), f[n]);
                        //af[n] += f[n];
                    }
                    //if (i==0) count++;

                }

            }
            if (pform == 1) {
                //for (n=0;n<3;n++) u[n] = dimg[n]/r;
            }     


            //__syncthreads();
        } // end pair j

        //if (i==0) printf("COUNT: %i\n",count);
    } // end if i<n (all threads)
}

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

    calculateForceKernel<<< dimGrid, dimBlock >>>(D,N,system.pbc.cutoff, dbasis, dreciprocal_basis, pform);
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
