#ifndef PSCF_CUDA_THREAD_GRID_TEST_H
#define PSCF_CUDA_THREAD_GRID_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/ThreadArray.h>
#include <pscf/cuda/ThreadMesh.h>
#include <cuda_runtime.h>

using namespace Util;
using namespace Pscf;

class CudaThreadGridTest : public UnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   void testThreadArray()
   {
      printMethod(TEST_FUNC);

      cudaDeviceProp dprop;
      // Get properties, assuming one GPU.
      cudaGetDeviceProperties(&dprop, 0);

      int nBlocks, nThreads;
      
      // Test with odd number of requested threads
      int nThreadsLogical(10001);
      ThreadArray::setThreadsLogical(nThreadsLogical, nBlocks, nThreads);
      ThreadArray::checkExecutionConfig();

      if (verbose()) {
         Log::file() << "\n  nThreadsLogical: " << nThreadsLogical
                     << "\n  nBlocks:         " << nBlocks
                     << "\n  nThreads:        " << nThreads << std::endl;
      }

      TEST_ASSERT(nThreads <= dprop.maxThreadsPerBlock);
      TEST_ASSERT(nThreads >= dprop.warpSize);
      TEST_ASSERT(nBlocks * nThreads > nThreadsLogical);
      TEST_ASSERT((nBlocks - 1) * nThreads < nThreadsLogical);
      TEST_ASSERT(ThreadArray::hasUnusedThreads());
      TEST_ASSERT(ThreadArray::nThreads() == nThreads);
      TEST_ASSERT(ThreadArray::nBlocks() == nBlocks);
      TEST_ASSERT(ThreadArray::nThreadsLogical() == nThreadsLogical);

      // Test with a power of two as number of requested threads
      nThreadsLogical = 16384;
      ThreadArray::setThreadsLogical(nThreadsLogical, nBlocks, nThreads);
      ThreadArray::checkExecutionConfig();

      if (verbose()) {
         Log::file() << "\n  nThreadsLogical: " << nThreadsLogical
                     << "\n  nBlocks:         " << nBlocks
                     << "\n  nThreads:        " << nThreads << std::endl;
      }

      TEST_ASSERT(nThreads <= dprop.maxThreadsPerBlock);
      TEST_ASSERT(nThreads >= dprop.warpSize);
      TEST_ASSERT(nBlocks * nThreads == nThreadsLogical);
      TEST_ASSERT(!ThreadArray::hasUnusedThreads());
      TEST_ASSERT(ThreadArray::nThreads() == nThreads);
      TEST_ASSERT(ThreadArray::nBlocks() == nBlocks);
      TEST_ASSERT(ThreadArray::nThreadsLogical() == nThreadsLogical);

      // Test with manual threads per block
      ThreadArray::setThreadsPerBlock(32);
      ThreadArray::setThreadsLogical(nThreadsLogical, nBlocks, nThreads);
      ThreadArray::checkExecutionConfig();

      if (verbose()) {
         Log::file() << "\n  nThreadsLogical: " << nThreadsLogical
                     << "\n  nBlocks:         " << nBlocks
                     << "\n  nThreads:        " << nThreads << std::endl;
      }

      TEST_ASSERT(nThreads == 32);
      TEST_ASSERT(nBlocks * nThreads == nThreadsLogical);
      TEST_ASSERT(!ThreadArray::hasUnusedThreads());
      TEST_ASSERT(ThreadArray::nThreads() == nThreads);
      TEST_ASSERT(ThreadArray::nBlocks() == nBlocks);
      TEST_ASSERT(ThreadArray::nThreadsLogical() == nThreadsLogical);
   }

   void testThreadMesh1D()
   {
      printMethod(TEST_FUNC);

      cudaDeviceProp dprop;
      // Get properties, assuming one GPU.
      cudaGetDeviceProperties(&dprop, 0);

      IntVec<1> mesh;
      dim3 blockDims, gridDims, meshDims;

      // Test with odd number of requested threads
      mesh[0] = dprop.warpSize * 3 + 3;
      ThreadMesh::setConfig(mesh, false);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x 
                     << "\n  Block dims: " << blockDims.x
                     << "\n  Grid dims:  " << gridDims.x << std::endl;
      }

      TEST_ASSERT(blockDims.x >= dprop.warpSize); 
      TEST_ASSERT(blockDims.x * gridDims.x > mesh[0]); 
      TEST_ASSERT(blockDims.x * (gridDims.x-1) < mesh[0]); 
      TEST_ASSERT(meshDims.x == mesh[0]);
      TEST_ASSERT(ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.y * blockDims.z == 1);
      TEST_ASSERT(gridDims.y * gridDims.z == 1);
      TEST_ASSERT(meshDims.y * meshDims.z == 1);

      // Test with a power of two as number of requested threads
      mesh[0] = dprop.warpSize * 128;
      ThreadMesh::setConfig(mesh, true);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x 
                     << "\n  Block dims: " << blockDims.x
                     << "\n  Grid dims:  " << gridDims.x << std::endl;
      }

      TEST_ASSERT(blockDims.x >= dprop.warpSize); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[0]); 
      TEST_ASSERT(meshDims.x == mesh[0]);
      TEST_ASSERT(!ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.y * blockDims.z == 1);
      TEST_ASSERT(gridDims.y * gridDims.z == 1);
      TEST_ASSERT(meshDims.y * meshDims.z == 1);

      // Test with manual threads per block
      mesh[0] = dprop.warpSize * 128;
      ThreadMesh::setConfig(mesh, true, 64);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x 
                     << "\n  Block dims: " << blockDims.x
                     << "\n  Grid dims:  " << gridDims.x << std::endl;
      }

      TEST_ASSERT(blockDims.x == 64); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[0]); 
      TEST_ASSERT(meshDims.x == mesh[0]);
      TEST_ASSERT(!ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.y * blockDims.z == 1);
      TEST_ASSERT(gridDims.y * gridDims.z == 1);
      TEST_ASSERT(meshDims.y * meshDims.z == 1);

      
   }

   void testThreadMesh2D()
   {
      printMethod(TEST_FUNC);

      cudaDeviceProp dprop;
      // Get properties, assuming one GPU.
      cudaGetDeviceProperties(&dprop, 0);

      IntVec<2> mesh;
      dim3 blockDims, gridDims, meshDims;

      // Test with odd number of requested threads in x, power of 2 in y
      mesh[0] = dprop.warpSize + 3;
      mesh[1] = dprop.warpSize * 2;
      ThreadMesh::setConfig(mesh, false);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " << meshDims.y
                     << "\n  Block dims: " << blockDims.x << " " << blockDims.y
                     << "\n  Grid dims:  " << gridDims.x << " " << gridDims.y 
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == dprop.warpSize / 2); 
      TEST_ASSERT(blockDims.x * gridDims.x >= mesh[0]); 
      TEST_ASSERT(blockDims.x * (gridDims.x-1) < mesh[0]); 
      TEST_ASSERT(blockDims.y * gridDims.y == mesh[1]); 
      TEST_ASSERT(blockDims.x * blockDims.y >= dprop.warpSize);
      TEST_ASSERT(meshDims.x == mesh[0]);
      TEST_ASSERT(meshDims.y == mesh[1]);
      TEST_ASSERT(ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.z == 1);
      TEST_ASSERT(gridDims.z == 1);
      TEST_ASSERT(meshDims.z == 1);

      // Test previous configuration but inverted
      ThreadMesh::setConfig(mesh, true);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " << meshDims.y
                     << "\n  Block dims: " << blockDims.x << " " << blockDims.y
                     << "\n  Grid dims:  " << gridDims.x << " " << gridDims.y 
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == mesh[1]); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[1]); 
      TEST_ASSERT(blockDims.y * gridDims.y >= mesh[0]); 
      TEST_ASSERT(blockDims.y * (gridDims.y-1) < mesh[0]); 
      TEST_ASSERT(meshDims.x == mesh[1]);
      TEST_ASSERT(meshDims.y == mesh[0]);
      TEST_ASSERT(ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.z == 1);
      TEST_ASSERT(gridDims.z == 1);
      TEST_ASSERT(meshDims.z == 1);

      // Test with manual threads per block
      ThreadMesh::setConfig(mesh, true, dprop.warpSize);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " << meshDims.y
                     << "\n  Block dims: " << blockDims.x << " " << blockDims.y
                     << "\n  Grid dims:  " << gridDims.x << " " << gridDims.y 
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == dprop.warpSize); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[1]); 
      TEST_ASSERT(blockDims.y == 1); 
      TEST_ASSERT(blockDims.y * gridDims.y == mesh[0]); 

      TEST_ASSERT(meshDims.x == mesh[1]);
      TEST_ASSERT(meshDims.y == mesh[0]);
      TEST_ASSERT(!ThreadMesh::hasUnusedThreads());

      TEST_ASSERT(blockDims.z == 1);
      TEST_ASSERT(gridDims.z == 1);
      TEST_ASSERT(meshDims.z == 1);
   }

   void testThreadMesh3D()
   {
      printMethod(TEST_FUNC);
      
      cudaDeviceProp dprop;
      // Get properties, assuming one GPU.
      cudaGetDeviceProperties(&dprop, 0);

      IntVec<3> mesh;
      dim3 blockDims, gridDims, meshDims;

      // Test with odd number of requested threads in x, power of 2 in y/z
      mesh[0] = dprop.warpSize + 3;
      mesh[1] = dprop.warpSize;
      mesh[2] = dprop.warpSize * 2;
      ThreadMesh::setConfig(mesh, false);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " 
                     << meshDims.y << " " << meshDims.z
                     << "\n  Block dims: " << blockDims.x << " " 
                     << blockDims.y << " " << blockDims.z
                     << "\n  Grid dims:  " << gridDims.x << " " 
                     << gridDims.y << " " << gridDims.z
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == dprop.warpSize / 2); 
      TEST_ASSERT(blockDims.x * gridDims.x >= mesh[0]); 
      TEST_ASSERT(blockDims.x * (gridDims.x-1) < mesh[0]); 
      TEST_ASSERT(blockDims.y * gridDims.y == mesh[1]); 
      TEST_ASSERT(blockDims.z * gridDims.z == mesh[2]); 
      TEST_ASSERT(blockDims.x * blockDims.y * blockDims.z >= dprop.warpSize);
      TEST_ASSERT(meshDims.x == mesh[0]);
      TEST_ASSERT(meshDims.y == mesh[1]);
      TEST_ASSERT(meshDims.z == mesh[2]);
      TEST_ASSERT(ThreadMesh::hasUnusedThreads());

      // Test previous configuration but inverted
      ThreadMesh::setConfig(mesh, true);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " 
                     << meshDims.y << " " << meshDims.z
                     << "\n  Block dims: " << blockDims.x << " " 
                     << blockDims.y << " " << blockDims.z
                     << "\n  Grid dims:  " << gridDims.x << " " 
                     << gridDims.y << " " << gridDims.z
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == mesh[2]); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[2]); 
      TEST_ASSERT(blockDims.y * gridDims.y == mesh[1]); 
      TEST_ASSERT(blockDims.z * gridDims.z >= mesh[0]); 
      TEST_ASSERT(blockDims.z * (gridDims.z-1) < mesh[0]); 
      TEST_ASSERT(meshDims.x == mesh[2]);
      TEST_ASSERT(meshDims.y == mesh[1]);
      TEST_ASSERT(meshDims.z == mesh[0]);
      TEST_ASSERT(!ThreadMesh::hasUnusedThreads());

      // Test with manual threads per block
      mesh[2] = dprop.warpSize * 4; // make larger
      ThreadMesh::setConfig(mesh, true, mesh[2]/2);

      blockDims = ThreadMesh::blockDims();
      gridDims = ThreadMesh::gridDims();
      meshDims = ThreadMesh::meshDims();

      if (verbose()) {
         Log::file() << "\n  Mesh dims:  " << meshDims.x << " " 
                     << meshDims.y << " " << meshDims.z
                     << "\n  Block dims: " << blockDims.x << " " 
                     << blockDims.y << " " << blockDims.z
                     << "\n  Grid dims:  " << gridDims.x << " " 
                     << gridDims.y << " " << gridDims.z
                     << std::endl;
      }

      TEST_ASSERT(blockDims.x == mesh[2]/2); 
      TEST_ASSERT(blockDims.x * gridDims.x == mesh[2]); 
      TEST_ASSERT(blockDims.y == 1); 
      TEST_ASSERT(blockDims.z == 1); 
      TEST_ASSERT(gridDims.y == meshDims.y); 
      TEST_ASSERT(gridDims.z == meshDims.z); 
      TEST_ASSERT(meshDims.x == mesh[2]);
      TEST_ASSERT(meshDims.y == mesh[1]);
      TEST_ASSERT(meshDims.z == mesh[0]);
      TEST_ASSERT(!ThreadMesh::hasUnusedThreads());
   }

};

TEST_BEGIN(CudaThreadGridTest)
TEST_ADD(CudaThreadGridTest, testThreadArray)
TEST_ADD(CudaThreadGridTest, testThreadMesh1D)
TEST_ADD(CudaThreadGridTest, testThreadMesh2D)
TEST_ADD(CudaThreadGridTest, testThreadMesh3D)
TEST_END(CudaThreadGridTest)

#endif