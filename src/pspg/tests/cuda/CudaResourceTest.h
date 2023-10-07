#ifndef PSPG_CUDA_RESOURCE_TEST_H
#define PSPG_CUDA_RESOURCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/gpu/RField.h>
#include <util/math/Constants.h>
#include <pspg/math/GpuResources.h>

#include <cstdlib>
#include <cmath>

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

   void testReductionSum() 
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 32;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks
      const int n = 32*nThreads*2;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Create device and host arrays
      cudaReal sum = 0;
      cudaReal sumCheck = 0;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, nBlocks*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         sumCheck+=num[i];
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionSum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionSum<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>(d_temp, d_temp, nBlocks);
      cudaMemcpy(&sum, d_temp, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(sum == sumCheck);
   }

   void testReductionMaxSmall()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 32;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks
      const int n = 1*nThreads*2;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Create device and host arrays
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
      reductionMax<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_max, d_num, n);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMaxLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 32;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks
      const int n = 8*nThreads*2;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Create device and host arrays
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, nBlocks*sizeof(cudaReal));
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMax<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMax<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>(d_max, d_temp, nBlocks);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMaxAbsLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 32;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks
      const int n = 8*nThreads*2;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Create device and host arrays
      cudaReal max = -1;
      cudaReal maxCheck = -10;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_max;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, nBlocks*sizeof(cudaReal));
      cudaMalloc((void**) &d_max, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (fabs(num[i]) > maxCheck) {
            maxCheck = fabs(num[i]);
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMaxAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMaxAbs<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>(d_max, d_temp, nBlocks);
      cudaMemcpy(&max, d_max, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(max == maxCheck);

   }

   void testReductionMinLarge()
   {
      printMethod(TEST_FUNC);

      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 32;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks
      const int n = 8*nThreads*2;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Create device and host arrays
      cudaReal min = 100000;
      cudaReal minCheck = 100000;
      cudaReal* num = new cudaReal[n];
      cudaReal* d_temp;
      cudaReal* d_min;
      cudaReal* d_num;
      cudaMalloc((void**) &d_num, n*sizeof(cudaReal));
      cudaMalloc((void**) &d_temp, nBlocks*sizeof(cudaReal));
      cudaMalloc((void**) &d_min, 1*sizeof(cudaReal));

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         if (num[i] < minCheck) {
            minCheck = num[i];
         }
      }

      // Launch kernel twice and get output
      cudaMemcpy(d_num, num, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      reductionMin<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>(d_temp, d_num, n);
      reductionMin<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>(d_min, d_temp, nBlocks);
      cudaMemcpy(&min, d_min, 1*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      TEST_ASSERT(min == minCheck);

   }

   void testGpuInnerProduct()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_a, d_b;
      d_a.allocate(n);
      d_b.allocate(n);

      // Host arrays
      cudaReal* a = new cudaReal[n];
      cudaReal* b = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++ ) {
         a[i] = (cudaReal)(std::rand() % 10000);
         b[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_a.cField(), a, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_b.cField(), b, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += a[i]*b[i];
      }

      // Inner product on device
      cudaReal prod = gpuInnerProduct(d_a.cField(),d_b.cField(),n);

      TEST_ASSERT(prodCheck==prod);
   }

   void testGpuSum()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += data[i];
      }

      // Inner product on device
      cudaReal prod = gpuSum(d_data.cField(),n);

      TEST_ASSERT(prodCheck==prod);
   }

   void testGpuMax()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (data[i] > maxCheck) maxCheck = data[i];
      }

      // Inner product on device
      cudaReal max = gpuMax(d_data.cField(),n);

      TEST_ASSERT(max==maxCheck);
   }

   void testGpuMaxAbs()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      cudaMemcpy(d_data.cField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) > maxCheck) maxCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal max = gpuMaxAbs(d_data.cField(),n);

      TEST_ASSERT(max==maxCheck);
   }

   void testGpuMin()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      cudaMemcpy(d_data.cField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal minCheck = 1000000;
      for (int i = 0; i < n; i++) {
         if (data[i] < minCheck) minCheck = data[i];
      }

      // Inner product on device
      cudaReal min = gpuMin(d_data.cField(),n);

      TEST_ASSERT(min==minCheck);
   }

   void testGpuMinAbs()
   {
      printMethod(TEST_FUNC);
      
      // GPU Resources
      ThreadGrid::init();
      const int nThreads = 128;
      ThreadGrid::setThreadsPerBlock(nThreads);

      // Data size and number of blocks.
      // Non-power-of-two to check performance in weird situations
      const int n = 14*nThreads + 77;
      int nBlocks;
      ThreadGrid::setThreadsLogical(n/2, nBlocks);

      // Device arrays
      Field<cudaReal> d_data;
      d_data.allocate(n);

      // Host arrays
      cudaReal* data = new cudaReal[n];
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      cudaMemcpy(d_data.cField(), data, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal minCheck = 1E300;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) < minCheck) minCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal min = gpuMinAbs(d_data.cField(),n);

      TEST_ASSERT(min==minCheck);
   }

};

TEST_BEGIN(CudaResourceTest)
TEST_ADD(CudaResourceTest, testReductionSum)
TEST_ADD(CudaResourceTest, testReductionMaxSmall)
TEST_ADD(CudaResourceTest, testReductionMaxLarge)
TEST_ADD(CudaResourceTest, testReductionMaxAbsLarge)
TEST_ADD(CudaResourceTest, testReductionMinLarge)
TEST_ADD(CudaResourceTest, testGpuInnerProduct)
TEST_ADD(CudaResourceTest, testGpuSum)
TEST_ADD(CudaResourceTest, testGpuMax)
TEST_ADD(CudaResourceTest, testGpuMaxAbs)
TEST_ADD(CudaResourceTest, testGpuMin)
TEST_ADD(CudaResourceTest, testGpuMinAbs)
TEST_END(CudaResourceTest)

#endif
