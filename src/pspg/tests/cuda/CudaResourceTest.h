#ifndef PSPG_CUDA_RESOURCE_TEST_H
#define PSPG_CUDA_RESOURCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/field/RDField.h>
#include <util/math/Constants.h>
#include <pspg/GpuResources.h>

#include <cstdlib>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class CudaResourceTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReductionMaxSmall()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 1; // parallel reduction into 1 block.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = 32;
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 100);
      }
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Host find max
      maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
         }
      }

      // Launch kernel and get output
      reductionMax<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(d_max, d_num, n);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMaxLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      NUMBER_OF_BLOCKS = 8; // parallel reduction into 32 blocks.
      THREADS_PER_BLOCK = 32;

      // Create device and host arrays
      const int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK;
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, NUMBER_OF_BLOCKS*sizeof(cudaReal));
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      int maxCheckIdx;
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
            maxCheckIdx = i;
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMax<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMax<<<1, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(cudaReal)>>>(d_max, d_temp, NUMBER_OF_BLOCKS);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

};

TEST_BEGIN(CudaResourceTest)
TEST_ADD(CudaResourceTest, testReductionMaxSmall)
TEST_ADD(CudaResourceTest, testReductionMaxLarge)
TEST_END(CudaResourceTest)

#endif
