/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ThreadMesh.h"
#include <cuda_runtime.h>

namespace Pscf {
namespace ThreadMesh {

   using namespace Util;

   namespace {

      /*
      * Anonymous namespace containing "static" variables and functions only 
      * used by global functions defined in namespace ThreadMesh. These are 
      * thus persistent pseudo-private variables and functions, much like 
      * private static class variables.
      */

      // Max threads per streaming multiprocessor, set by querying hardware.
      int MAX_THREADS_PER_SM = -1;

      // Max threads per block, set by querying hardware and optimizing.
      int MAX_THREADS_PER_BLOCK = -1;

      // Threads per block in current configuration.
      // This may be either user-defined or equal to MAX_THREADS_PER_BLOCK.
      int THREADS_PER_BLOCK = -1;

      // Number of threads per warp.
      int WARP_SIZE = -1;

      // Multidimensional array of threads requested for execution.
      dim3 MESH_DIMS(0, 0, 0);

      // Dimensions of multidimensional thread block for execution. 
      dim3 BLOCK_DIMS(0, 0, 0);

      // Dimensions of multidimensional grid of blocks for execution. 
      dim3 GRID_DIMS(0, 0, 0);

      // Will threads go unused?
      bool UNUSED_THREADS;

      /*
      * Initialize the ThreadMesh namespace
      */
      void init()
      {
         // Check that a CUDA device is available.
         int count = 0;
         cudaGetDeviceCount(&count);

         if (count == 0) {
            UTIL_THROW("No CUDA devices found.");
         } else if (count > 1) {
            Log::file() << "\nWarning: multiple GPUs detected.\n"
               << "This program is not compatible with multiple devices.\n"
               << "Only the first device will be used." << std::endl;
         }

         // Get properties, assuming one GPU.
         cudaDeviceProp dprop;
         cudaGetDeviceProperties(&dprop, 0);
         WARP_SIZE = dprop.warpSize;
         MAX_THREADS_PER_SM = dprop.maxThreadsPerMultiProcessor;

         // Find the highest power of two that evenly divides into the
         // maximum number of threads per streaming multiprocessor
         // This will lead to the highest occupancy!
         MAX_THREADS_PER_BLOCK = (MAX_THREADS_PER_SM & 
                                 (~(MAX_THREADS_PER_SM - 1)));
         
         // If this value is > maxThreadsPerBlock, reduce
         while (MAX_THREADS_PER_BLOCK > dprop.maxThreadsPerBlock) {
            MAX_THREADS_PER_BLOCK /= 2;
         }
      }

   }

   template <int D>
   void setConfig(IntVec<D> const & meshDims, bool invert, int blockSize)
   {
      // Verify that requested thread grid is valid
      for (int i = 0; i < D; i++) {
         UTIL_CHECK(meshDims[i] > 0);
      }
      
      // If MAX_THREADS_PER_BLOCK hasn't been set at all, initialize
      if (MAX_THREADS_PER_BLOCK == -1) {
         init();
      }

      // Check if requested number of threads matches the previous request
      bool match = true;
      if ((blockSize > 0) && (blockSize != THREADS_PER_BLOCK)) {
         match = false;
      } else if (invert) {
         if (meshDims[D-1] != MESH_DIMS.x) match = false;
         if (D > 1) { if (meshDims[D-2] != MESH_DIMS.y) match = false; }
         if (D > 2) { if (meshDims[D-3] != MESH_DIMS.z) match = false; }
      } else {
         if (meshDims[0] != MESH_DIMS.x) match = false;
         if (D > 1) { if (meshDims[1] != MESH_DIMS.y) match = false; }
         if (D > 2) { if (meshDims[2] != MESH_DIMS.z) match = false; }
      }

      if (match) {
         // Do nothing. Reuse previous execution configuration.
         return;
      }

      // Store meshDims in MESH_DIMS in dim3 format
      MESH_DIMS.y = 1;
      MESH_DIMS.z = 1;
      if (invert) {
         MESH_DIMS.x = meshDims[D-1];
         if (D > 1) MESH_DIMS.y = meshDims[D-2];
         if (D > 2) MESH_DIMS.z = meshDims[D-3];
      } else {
         MESH_DIMS.x = meshDims[0];
         if (D > 1) MESH_DIMS.y = meshDims[1];
         if (D > 2) MESH_DIMS.z = meshDims[2];
      }
      
      // Assign the total block size parameter if not provided
      bool manualBlockSize = false;
      if (blockSize < 0) { 
         blockSize = MAX_THREADS_PER_BLOCK;
      } else {
         manualBlockSize = true;
         if ((blockSize & (blockSize - 1)) != 0) { // blockSize not power of 2
            UTIL_THROW("Manual block size entry must be a power of 2.");
         }
         if (blockSize < WARP_SIZE) {
            Log::file() << "\nRequested threads per block: " << blockSize
                        << "\nWarp size: " << WARP_SIZE << std::endl;
            UTIL_THROW("Threads per block cannot be smaller than warp size.");
         }
      }
      THREADS_PER_BLOCK = blockSize;

      // Set BLOCK_DIMS.x to maximize opportunities for coalescence.
      // (Note: we restrict block dimensions to powers of 2.)
      if (MESH_DIMS.x % WARP_SIZE == 0) {
         // If MESH_DIMS.x is a multiple of WARP_SIZE, BLOCK_DIMS.x can be 
         // set to WARP_SIZE, resulting in no unused threads along x and
         // optimal coalescence.
         BLOCK_DIMS.x = WARP_SIZE;
      } else if ((MESH_DIMS.x & (MESH_DIMS.x - 1)) == 0) {
         // MESH_DIMS.x is a small power of two, so we can set BLOCK_DIMS.x
         // equal to MESH_DIMS.x and have no unused threads along x, 
         // resulting in optimal coalescence.
         BLOCK_DIMS.x = MESH_DIMS.x;
      } else {
         // MESH_DIMS.x is not a power of 2 or multiple of WARP_SIZE.
         // Choose block dimensions to be either a warp or half-warp in the
         // x dimension depending on the number of unused threads.
         if ((MESH_DIMS.x / (WARP_SIZE/2)) % 2 == 0) {
            // half warp will have smaller number of unused threads
            BLOCK_DIMS.x = WARP_SIZE / 2;
         } else {
            // half and full warp will have same number of unused threads
            BLOCK_DIMS.x = WARP_SIZE;
         }
      }

      // Set BLOCK_DIMS.y to lowest power of 2 that is >= MESH_DIMS.y
      BLOCK_DIMS.y = 1;
      while (BLOCK_DIMS.y < MESH_DIMS.y) {
         BLOCK_DIMS.y *= 2;
      }
      GRID_DIMS.y = 1; // set GRID_DIMS.y

      // Set BLOCK_DIMS.z to lowest power of 2 that is >= MESH_DIMS.z
      BLOCK_DIMS.z = 1;
      while (BLOCK_DIMS.z < MESH_DIMS.z) {
         BLOCK_DIMS.z *= 2;
      }
      
      // BLOCK_DIMS.z cannot exceed 64. Adjust if needed.
      if (BLOCK_DIMS.z > 64) {
         BLOCK_DIMS.z = 64;
      }
      
      // Set GRID_DIMS.z
      GRID_DIMS.z = MESH_DIMS.z / BLOCK_DIMS.z;
      if (MESH_DIMS.z % BLOCK_DIMS.z) GRID_DIMS.z++;

      // Shrink BLOCK_DIMS until block size size equals THREADS_PER_BLOCK
      blockSize = BLOCK_DIMS.x * BLOCK_DIMS.y * BLOCK_DIMS.z;
      if (blockSize >= THREADS_PER_BLOCK) {
         
         while (blockSize > THREADS_PER_BLOCK) {
            // Determine how many unused gridpoints lie along y and z
            int yOvershoot = BLOCK_DIMS.y * GRID_DIMS.y - MESH_DIMS.y;
            int zOvershoot = BLOCK_DIMS.z * GRID_DIMS.z - MESH_DIMS.z;

            // Halve the element of BLOCK_DIMS with largest overshoot
            if ((yOvershoot > zOvershoot) || (BLOCK_DIMS.z == 1)) {
               UTIL_CHECK(BLOCK_DIMS.y > 1);
               BLOCK_DIMS.y /= 2;
               GRID_DIMS.y = MESH_DIMS.y / BLOCK_DIMS.y;
               if (MESH_DIMS.y % BLOCK_DIMS.y) GRID_DIMS.y++;
            } else {
               BLOCK_DIMS.z /= 2;
               GRID_DIMS.z = MESH_DIMS.z / BLOCK_DIMS.z;
               if (MESH_DIMS.z % BLOCK_DIMS.z) GRID_DIMS.z++;
            }

            blockSize = BLOCK_DIMS.x * BLOCK_DIMS.y * BLOCK_DIMS.z;
         }
         UTIL_CHECK(blockSize == THREADS_PER_BLOCK);
         UTIL_CHECK(blockSize % WARP_SIZE == 0);

      } else {
         // blockSize is less than THREADS_PER_BLOCK.
         // BLOCK_DIMS.x can be increased.
         while ((BLOCK_DIMS.x < MESH_DIMS.x) && (blockSize < THREADS_PER_BLOCK))
         {
            BLOCK_DIMS.x *= 2;
            blockSize *= 2;
         }
         if (blockSize < WARP_SIZE) {
            // Make sure blockSize is at least one warp
            UTIL_CHECK(blockSize > 0);
            while (blockSize < WARP_SIZE) {
               BLOCK_DIMS.x *= 2;
               blockSize *= 2;
            }
         }

         THREADS_PER_BLOCK = blockSize;
      }

      // Set GRID_DIMS.x
      GRID_DIMS.x = MESH_DIMS.x / BLOCK_DIMS.x;
      if (MESH_DIMS.x % BLOCK_DIMS.x) GRID_DIMS.x++;

      // If block is smaller than manually requested size, print warning
      if ((manualBlockSize) && (blockSize < THREADS_PER_BLOCK)) {
         Log::file() << "WARNING: The number of threads per block (" 
                  << blockSize
                  << ") will be smaller than \nthe requested size of "
                  << THREADS_PER_BLOCK << "." << std::endl;
      }

      // Set UNUSED_THREADS
      if ((MESH_DIMS.x % BLOCK_DIMS.x) || (MESH_DIMS.y % BLOCK_DIMS.y) || 
                                          (MESH_DIMS.z % BLOCK_DIMS.z)) {
         UNUSED_THREADS = true;
      } else {
         UNUSED_THREADS = false;
      }

      checkConfig();
   }

   // Instantiate setConfig methods for D = 1, 2, 3
   template void setConfig<1>(IntVec<1> const & meshDims, 
                              bool invert, int blockSize);
   template void setConfig<2>(IntVec<2> const & meshDims, 
                              bool invert, int blockSize);
   template void setConfig<3>(IntVec<3> const & meshDims, 
                              bool invert, int blockSize);

   void checkConfig()
   {
      // Check that max threads per block is a power of two. 
      if ((MAX_THREADS_PER_BLOCK & (MAX_THREADS_PER_BLOCK - 1)) != 0) {
         Log::file() << "\nMax threads per block: " << MAX_THREADS_PER_BLOCK
                     << std::endl;
         UTIL_THROW("Max threads per block must be a power of two.");
      }

      // Check that max threads per block is multiple of WARP_SIZE.
      if (MAX_THREADS_PER_BLOCK % WARP_SIZE != 0)
      {
         Log::file() << "\nMax threads per block: " << MAX_THREADS_PER_BLOCK
                     << "\nWarp size: " << WARP_SIZE << std::endl;
         UTIL_THROW("Max threads per block must be a multiple of warp size.");
      }

      // Check that threads per block is multiple of WARP_SIZE.
      if (THREADS_PER_BLOCK % WARP_SIZE != 0)
      {
         Log::file() << "\nThreads per block: " << THREADS_PER_BLOCK
                     << "\nWarp size: " << WARP_SIZE << std::endl;
         UTIL_THROW("Threads per block must be a multiple of warp size.");
      }

      // Check that block dimensions are each 1 or greater
      if ((BLOCK_DIMS.x < 1) || (BLOCK_DIMS.y < 1) || (BLOCK_DIMS.z < 1)) {
         Log::file() << "\nBlock dimensions: " << BLOCK_DIMS.x << " "
                     << BLOCK_DIMS.y << " " << BLOCK_DIMS.z << std::endl;
         UTIL_THROW("Block dimensions must each be a power of two.");
      }
      
      // Check that BLOCK_DIMS multiply to THREADS_PER_BLOCK
      if (BLOCK_DIMS.x * BLOCK_DIMS.y * BLOCK_DIMS.z != THREADS_PER_BLOCK) {
         UTIL_THROW("THREADS_PER_BLOCK not properly set.");
      }

      // Check that the thread grid is at least as big as the requested mesh
      if ((BLOCK_DIMS.x * GRID_DIMS.x < MESH_DIMS.x) ||
          (BLOCK_DIMS.y * GRID_DIMS.y < MESH_DIMS.y) ||
          (BLOCK_DIMS.z * GRID_DIMS.z < MESH_DIMS.z)) {
         Log::file() << "\nBlock dimensions: " << BLOCK_DIMS.x << " "
                     << BLOCK_DIMS.y << " " << BLOCK_DIMS.z << std::endl;
         Log::file() << "\nGrid dimensions: " << GRID_DIMS.x << " "
                     << GRID_DIMS.y << " " << GRID_DIMS.z << std::endl;
         Log::file() << "\nMesh dimensions: " << MESH_DIMS.x << " "
                     << MESH_DIMS.y << " " << MESH_DIMS.z << std::endl;
         UTIL_THROW("Thread grid smaller than the requested mesh.");
      }

      // Check that the maximum number of threads per multiprocessor is an 
      // integer multiple of the threads per block. This is not required 
      // for validity, but performance will be suboptimal if not the case, 
      // as it will limit the total number of threads that can be 
      // scheduled at any given time.
      if (MAX_THREADS_PER_SM % MAX_THREADS_PER_BLOCK != 0) {
         Log::file() << "WARNING: The number of threads per block (" 
                     << MAX_THREADS_PER_BLOCK 
                     << ") is not an even divisor of the maximum number"
                     << " of threads per streaming multiprocessor ("
                     << MAX_THREADS_PER_SM
                     << "). Performance will be suboptimal." 
                     << std::endl;
      }
   }

   /**
   * Manually set the block size that should be used by default.
   */
   void setThreadsPerBlock(int blockSize)
   {  
      // If MAX_THREADS_PER_BLOCK hasn't been set at all, initialize
      if (MAX_THREADS_PER_BLOCK == -1) {
         init();
      }

      if ((blockSize & (blockSize - 1)) != 0) { // blockSize not power of 2
         UTIL_THROW("Manual block size entry must be a power of 2.");
      }

      if (blockSize < WARP_SIZE) { 
         Log::file() << "\nRequested threads per block: " << blockSize
                     << "\nWarp size: " << WARP_SIZE << std::endl;
         UTIL_THROW("Threads per block cannot be smaller than warp size.");
      }

      MAX_THREADS_PER_BLOCK = blockSize; 
   }

   // Accessors

   dim3 gridDims()
   {  return GRID_DIMS; }

   dim3 blockDims()
   {  return BLOCK_DIMS; }

   dim3 meshDims()
   { return MESH_DIMS; }

   int warpSize()
   {  return WARP_SIZE; }

   bool hasUnusedThreads()
   {  return UNUSED_THREADS; }

} // namespace ThreadMesh
} // namespace Pscf