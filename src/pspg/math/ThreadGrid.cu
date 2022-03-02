#ifndef PSPG_THREADGRID_CU
#define PSPG_THREADGRID_CU

#include "ThreadGrid.h"

namespace {
   // Maximum threads per block, either set by querying hardware or by user.
   int MAX_THREADS_PER_BLOCK = -1;
   // Number of threads per block for execution. Determined by setThreadsLogical.
   int THREADS_PER_BLOCK = -1;
   // Number of blocks for execution. Determined by setThreadsLogical.
   int BLOCKS = -1;
   // Total number of threads requested for execution. Set by setThreadsLogical.
   int THREADS_LOGICAL = -1;
   // Will threads go unused?
   bool UNUSED_THREADS;
}

namespace Pscf {
namespace Pspg {
namespace ThreadGrid {

   using namespace Util;

   void init()
   {
      // Check that a CUDA device is available.
      int count = 0;
      cudaGetDeviceCount(&count);

      if (count == 0) 
         UTIL_THROW("No CUDA devices found.");

      // Set a default maximum threads per block by querying hardware.
      setThreadsPerBlock();
   }

   void setThreadsPerBlock()
   {
      cudaDeviceProp dprop;
      // get properties, assuming one GPU.
      cudaGetDeviceProperties(&dprop, 0);
      int maxThPerSM = dprop.maxThreadsPerMultiProcessor;

      // Find the highest power of two that evenly divides into the
      // maximum number of threads per streaming multiprocessor
      // This will lead to the highest occupancy!

      int threadsPerBlock = (maxThPerSM & (~(maxThPerSM - 1)));
      
      // Check for validity:
      while (true) {
         if (threadsPerBlock > dprop.maxThreadsPerBlock)
            threadsPerBlock /= 2;
         else
            break;
      }

      setThreadsPerBlock(threadsPerBlock);
   }

   void setThreadsPerBlock(int const nThreadsPerBlock)
   {
      MAX_THREADS_PER_BLOCK = nThreadsPerBlock;
      BLOCKS = 0;
      THREADS_LOGICAL = 0;
      checkExecutionConfig();
   }

   void setThreadsLogical(int const nThreadsLogical)
   {
      // Verify that requested threads is valid (greater than 0).
      UTIL_ASSERT(nThreadsLogical > 0);
      
      if (MAX_THREADS_PER_BLOCK == -1) // if max_threads_per_block hasn't been set at all, initialize.
         init();

      // Check if the number of requested threads matches the previous number of request threads
      if (THREADS_LOGICAL == nThreadsLogical) {
         // Do nothing. Previous execution configuration will be used.
         return;
      }

      // Set the number of total requested threads.
      THREADS_LOGICAL = nThreadsLogical;

      // Compute the execution configuration. Number of blocks rounded up to the nearest integer.
      THREADS_PER_BLOCK = MAX_THREADS_PER_BLOCK;
      BLOCKS = ceil(double(nThreadsLogical)/double(THREADS_PER_BLOCK));

      // Determine if there will be unused threads
      UNUSED_THREADS = (BLOCKS*THREADS_PER_BLOCK > THREADS_LOGICAL);

   }

   void setThreadsLogical(int const nThreadsLogical, int & nBlocks)
   {
      setThreadsLogical(nThreadsLogical);

      nBlocks = BLOCKS;
   }

   void setThreadsLogical(int const nThreadsLogical, int & nBlocks, int & nThreads)
   {
      setThreadsLogical(nThreadsLogical);

      nBlocks = BLOCKS;
      nThreads = THREADS_PER_BLOCK;
   }

   void checkExecutionConfig()
   {
      // get relevant device hardware properties, assuming one device.
      cudaDeviceProp dprop;
      cudaGetDeviceProperties(&dprop, 0);
      int warpSize = dprop.warpSize;
      int maxThreadsPerMultiProcessor = dprop.maxThreadsPerMultiProcessor;

      // Check that threads per block is a power of two. This is required for parallel reductions.
      if ((MAX_THREADS_PER_BLOCK & (MAX_THREADS_PER_BLOCK - 1)) != 0) {
         UTIL_THROW("Set number of threads per block must be a power of two.");
      }

      // Check that threads per block is multiple of warpSize. This is required
      // because a warp is generally 32.
      if (MAX_THREADS_PER_BLOCK%warpSize != 0)
      {
         char buffer[100];
         sprintf(buffer, "Number of threads per block must be a multiple of warp size %d.\n",warpSize);
         UTIL_THROW(buffer);
      }

      // Check that the maximum number of threads per multiprocessor is an integer multiple
      // of the threads per block. This is not required for validity, but performance will be suboptimal
      // if not the case as it will limit the total number of threads that can be scheduled at any
      // given time.
      if (maxThreadsPerMultiProcessor % MAX_THREADS_PER_BLOCK%warpSize != 0)
      {
         std::cerr << "WARNING: The specified number of threads per block (" << MAX_THREADS_PER_BLOCK 
                  << ") is not an even divisor of the maximum number of threads per streaming multiprocessor ("
                  << maxThreadsPerMultiProcessor << "). Performance will be suboptimal." << std::endl;
      }

   }

   // Accessors

   int nBlocks()
   { return BLOCKS; }

   int nThreads()
   { return THREADS_PER_BLOCK; }

   int nThreadsLogical()
   { return THREADS_LOGICAL; }

   bool hasUnusedThreads()
   { return UNUSED_THREADS; }

}
}
}
#endif