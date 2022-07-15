#ifndef PSPG_KERNEL_WRAPPERS_CU
#define PSPG_KERNEL_WRAPPERS_CU

#include "KernelWrappers.h"
#include "ThreadGrid.h"

namespace Pscf {
namespace Pspg {

__host__ cudaReal gpuSum(cudaReal const * d_in, int size) 
{
   
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction sum operation
   reductionSum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_in, size);

   // Load the partially reduced sum result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

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

__host__ cudaReal gpuInnerProduct(cudaReal const * d_a, cudaReal const * d_b, int size)
{
   
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction inner product operation
   reductionInnerProduct<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_a, d_b, size);

   // Load the partially reduced inner product result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

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

__host__ cudaReal gpuMax(cudaReal const * d_in, int size)
{
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction maximum operation
   reductionMax<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_in, size);

   // Load the partially reduced maximum result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

   // Perform final comparison on CPU
   cudaReal max = 0;
   for (int i = 0; i < nBlocks; i++) {
      if (temp_[i] > max) max = temp_[i];
   }

   return max;
}

__host__ cudaReal gpuMaxAbs(cudaReal const * d_in, int size)
{
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction maximum of the absolute value
   reductionMaxAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_in, size);

   // Load the partially reduced maximum of the absolute value result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

   // Perform final comparison on CPU. Absolute values already taken in kernel, so 
   // no need to take them here.
   cudaReal max = 0;
   for (int i = 0; i < nBlocks; i++) {
      if (temp_[i] > max) max = temp_[i];
   }

   return max;
}

__host__ cudaReal gpuMin(cudaReal const * d_in, int size)
{
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction minimum operation
   reductionMin<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_in, size);

   // Load the partially reduced minimum result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

   // Perform final comparison on CPU
   cudaReal min = temp_[0];
   for (int i = 1; i < nBlocks; i++) {
      if (temp_[i] < min) min = temp_[i];
   }

   return min;
}

__host__ cudaReal gpuMinAbs(cudaReal const * d_in, int size)
{
   // Establish GPU resources for this parallel reduction. Divided by two because
   // of the global memory load in the kernel performing the first level of reduction!
   int nBlocks, nThreads;
   int halvedsize = ceil((float)size/2);

   ThreadGrid::setThreadsLogical(halvedsize,nBlocks,nThreads);

   // Set up temporary, small device and host arrays for storing reduced data
   // temp_ will handle excess and reduced and sum them up
   cudaReal* temp_ = new cudaReal[nBlocks];
   cudaReal* d_temp_;
   gpuErrchk(cudaMalloc((void**) &d_temp_, nBlocks*sizeof(cudaReal)));

   // Perform parallel reduction minimum of the absolute value
   reductionMinAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp_, d_in, size);

   // Load the partially reduced minimum of the absolute value result to temp_
   gpuErrchk(cudaMemcpy(temp_, d_temp_, nBlocks*sizeof(cudaReal), cudaMemcpyDeviceToHost));

   // Perform final comparison on CPU
   cudaReal min = temp_[0];
   for (int i = 1; i < nBlocks; i++) {
      if (temp_[i] < min) min = temp_[i];
   }

   return min;
}

}
}
#endif
