#ifndef PRDC_CUDA_RESOURCE_TEST_H
#define PRDC_CUDA_RESOURCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/GpuResources.h>
#include <util/math/Constants.h>

#include <cstdlib>
#include <cmath>

using namespace Util;
using namespace Pscf;

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
      cudaReal sumCheck = 0;
      HostDArray<cudaReal> num(n), sum(1);
      DeviceArray<cudaReal> d_num(n), d_temp(nBlocks), d_sum(1);

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 10000);
      }

      // Host find max
      for (int i = 0; i < n; i++) {
         sumCheck+=num[i];
      }

      // Launch kernel twice and get output
      d_num = num;
      reductionSum<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (d_temp.cArray(), d_num.cArray(), n);
      reductionSum<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>
                              (d_sum.cArray(), d_temp.cArray(), nBlocks);
      sum = d_sum; // copy first element of d_sum from device

      TEST_ASSERT(sum[0] == sumCheck);
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
      cudaReal maxCheck = -10;
      HostDArray<cudaReal> num(n), max(1);
      DeviceArray<cudaReal> d_num(n), d_max(1);

      // Test data
      for (int i = 0; i < n; i++) {
         num[i] = (cudaReal)(std::rand() % 100);
      }

      // Host find max
      maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (num[i] > maxCheck) {
            maxCheck = num[i];
         }
      }

      // Launch kernel and get output
      d_num = num;
      reductionMax<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                 (d_max.cArray(), d_num.cArray(), n);
      max = d_max; // copy first element of d_max from device

      TEST_ASSERT(max[0] == maxCheck);

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
      cudaReal maxCheck = -10;
      HostDArray<cudaReal> num(n), max(1);
      DeviceArray<cudaReal> d_num(n), d_temp(nBlocks), d_max(1);

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
      d_num = num;
      reductionMax<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                 (d_temp.cArray(), d_num.cArray(), n);
      reductionMax<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>
                           (d_max.cArray(), d_temp.cArray(), nBlocks);
      max = d_max;

      TEST_ASSERT(max[0] == maxCheck);

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
      cudaReal maxCheck = -10;
      HostDArray<cudaReal> num(n), max(1);
      DeviceArray<cudaReal> d_num(n), d_temp(nBlocks), d_max(1);

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
      d_num = num;
      reductionMaxAbs<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                    (d_temp.cArray(), d_num.cArray(), n);
      reductionMaxAbs<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>
                              (d_max.cArray(), d_temp.cArray(), nBlocks);
      max = d_max;

      TEST_ASSERT(max[0] == maxCheck);

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
      cudaReal minCheck = 100000;
      HostDArray<cudaReal> num(n), min(1);
      DeviceArray<cudaReal> d_num(n), d_temp(nBlocks), d_min(1);

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
      d_num = num;
      reductionMin<<<nBlocks, nThreads, nThreads*sizeof(cudaReal)>>>
                                 (d_temp.cArray(), d_num.cArray(), n);
      reductionMin<<<1, nBlocks/2, nBlocks/2*sizeof(cudaReal)>>>
                           (d_min.cArray(), d_temp.cArray(), nBlocks);
      min = d_min; 

      TEST_ASSERT(min[0] == minCheck);

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
      DeviceArray<cudaReal> d_a(n), d_b(n);

      // Host arrays
      HostDArray<cudaReal> a(n), b(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++ ) {
         a[i] = (cudaReal)(std::rand() % 10000);
         b[i] = (cudaReal)(std::rand() % 10000);
      }
      d_a = a; 
      d_b = b;

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += a[i]*b[i];
      }

      // Inner product on device
      cudaReal prod = gpuInnerProduct(d_a.cArray(),d_b.cArray(),n);

      TEST_ASSERT(prodCheck == prod);
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
      DeviceArray<cudaReal> d_data(n);

      // Host arrays
      HostDArray<cudaReal> data(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      d_data = data;

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += data[i];
      }

      // Inner product on device
      cudaReal prod = gpuSum(d_data.cArray(),n);

      TEST_ASSERT(prodCheck == prod);
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
      DeviceArray<cudaReal> d_data(n);

      // Host arrays
      HostDArray<cudaReal> data(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      d_data = data;

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (data[i] > maxCheck) maxCheck = data[i];
      }

      // Inner product on device
      cudaReal max = gpuMax(d_data.cArray(),n);

      TEST_ASSERT(max == maxCheck);
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
      DeviceArray<cudaReal> d_data(n);

      // Host arrays
      HostDArray<cudaReal> data(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      d_data = data;

      // Inner product on host
      cudaReal maxCheck = 0;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) > maxCheck) maxCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal max = gpuMaxAbs(d_data.cArray(),n);

      TEST_ASSERT(max == maxCheck);
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
      DeviceArray<cudaReal> d_data(n);

      // Host arrays
      HostDArray<cudaReal> data(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000);
      }
      d_data = data;

      // Inner product on host
      cudaReal minCheck = 1000000;
      for (int i = 0; i < n; i++) {
         if (data[i] < minCheck) minCheck = data[i];
      }

      // Inner product on device
      cudaReal min = gpuMin(d_data.cArray(),n);

      TEST_ASSERT(min == minCheck);
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
      DeviceArray<cudaReal> d_data(n);

      // Host arrays
      HostDArray<cudaReal> data(n);
      
      // Create random data, store on host and device
      for (int i = 0; i < n; i++) {
         data[i] = (cudaReal)(std::rand() % 10000 - 6000);
      }
      d_data = data;

      // Inner product on host
      cudaReal minCheck = 1E300;
      for (int i = 0; i < n; i++) {
         if (fabs(data[i]) < minCheck) minCheck = fabs(data[i]);
      }

      // Inner product on device
      cudaReal min = gpuMinAbs(d_data.cArray(),n);

      TEST_ASSERT(min == minCheck);
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
