#include <stdio.h>
#include <cuda.h>

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

__global__
void saxpy(int n, float a, float *x, float *y)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if (i < n) {

        y[i] = a*x[i] + y[i];
        printf("y[%i] = %f\n", i, y[i]);
    }
}

void CUDA_force(System &system) {

    // first send information to the device
    int system_size = sizeof(system);
    printf("SYSTEM SIZE: %i\n", system_size); 

    int atoms_size = system.constants.total_atoms * sizeof(system.molecules[0].atoms[0]);
    printf("molecules vector size: %i\n", atoms_size);

    printf("sizeof(Atom): %i\n", (int)sizeof(Atom));
}
