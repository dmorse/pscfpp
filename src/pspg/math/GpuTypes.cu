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
   nThreads = MAX_THREADS_PER_BLOCK;
   nBlocks = ceil((float)datasize/MAX_THREADS_PER_BLOCK);
}

}
}
#endif