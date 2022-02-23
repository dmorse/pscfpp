#ifndef PSPG_KERNEL_WRAPPERS_CU
#define PSPG_KERNEL_WRAPPERS_CU

#include "KernelWrappers.h"

namespace Pscf {
namespace Pspg {

// __host__ cudaReal gpuSum(cudaReal* d_in, int size) 
// {

// }

__host__ cudaReal gpuInnerProduct(cudaReal* d_a, cudaReal* d_b, int size)
{

   // Check to see if data size is a power of two
   int nPow2 = size;
   int nExcess = 0;
   if ( (size & (size - 1)) != 0 ) {
      // if not, handle only up to the nearest lower power of two with parallel reduction
      // handle the rest with the CPU
      nPow2 = pow(2,floor(log2(size)));
      nExcess = size-nPow2;
   }
   
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   setGpuBlocksThreads(nBlocks,nThreads,nPow2/2);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks + nExcess];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   std::cout << nPow2 << "   " << nExcess << "   " << nBlocks << "   " << nThreads << std::endl;
   // Perform parallel reduction inner product
   reductionInnerProduct<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_a, d_b, nPow2);

   // Checks first for kernel launch errors and then for kernel runtime errors
   gpuErrchk( cudaPeekAtLastError() );
   gpuErrchk( cudaDeviceSynchronize() );

   // Load the partially reduced inner product result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

   // If there are excess, load those from device and perform elementwise multiplication on CPU
   if (nExcess != 0) {
      cudaReal* a_excess = new cudaReal[nExcess];
      cudaReal* b_excess = new cudaReal[nExcess];
      cudaMemcpy(a_excess, d_a + nPow2, nExcess*sizeof(cudaReal), cudaMemcpyDeviceToHost);
      cudaMemcpy(b_excess, d_b + nPow2, nExcess*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      for (int i = 0; i < nExcess; i++) {
         temp_[nPow2+i] = a_excess[i] * b_excess[i];
      }
   }

   // Perform final sum on CPU using kahan summation to reduce accumulation of error
   cudaReal sum = 0, tempVal, tempSum;
   cudaReal err = 0;
   for (int i = 0; i < nBlocks; ++i) {
      tempVal = temp_[i] - err;
      tempSum = sum + tempVal;
      err = tempSum - sum - tempVal;
      sum = tempSum;  
   }

   return sum;
}

// __host__ cudaReal gpuMax(cudaReal* d_in, int size)
// {

// }

// __host__ cudaReal gpuMaxAbs(cudaReal* d_in, int size)
// {

// }

// __host__ cudaReal gpuMin(cudaReal* d_in, int size)
// {

// }

// __host__ cudaReal gpuMinAbs(cudaReal* d_in, int size)
// {

// }

}
}
#endif
