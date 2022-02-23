#ifndef PSPG_GPU_TYPES_CU
#define PSPG_GPU_TYPES_CU

#include "GpuTypes.h"

namespace Pscf {
namespace Pspg {

int THREADS_PER_BLOCK;
int NUMBER_OF_BLOCKS;
int MAX_THREADS_PER_BLOCK;

void setGpuBlocksThreads(int & nBlocks, int & nThreads, int const & datasize)
{
   // Check if datasize is a power of two
   bool isPowerOfTwo = (datasize & (datasize - 1))==0;
   int roundedSize;

   if (!isPowerOfTwo) {
      // if not a power of two, compute the closest greater power of two
      roundedSize = pow(2, ceil(log2(datasize)));
   }
   else {
      roundedSize = datasize;
   }

   if (roundedSize % MAX_THREADS_PER_BLOCK == 0) {
      // If size divided by max threads gives an integer,
      // then the number of blocks can be that integer
      nBlocks = roundedSize/MAX_THREADS_PER_BLOCK;
      nThreads = MAX_THREADS_PER_BLOCK;
   } else {
      // If size divided by max threads does not give an integer,
      // then because both are a power of two in this situation, this
      // means that size is smaller than max threads. 
      nBlocks = 1;
      nThreads = roundedSize;
   }
}

}
}
#endif