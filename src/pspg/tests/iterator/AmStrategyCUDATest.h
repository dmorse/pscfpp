#ifndef PSPG_AM_STRATEGY_CUDA_TEST_H
#define PSPG_AM_STRATEGY_CUDA_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspg/System.h>
#include <pspg/iterator/AmStrategyCUDA.h>

#include <pspg/field/RDField.h>
#include <util/math/Constants.h>
#include <pspg/math/GpuResources.h>

#include <fstream>
#include <iomanip>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspg;

class AmStrategyCUDATest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testInnerProduct()
   {
      printMethod(TEST_FUNC);
      // GPU Resources
      NUMBER_OF_BLOCKS = 16;
      THREADS_PER_BLOCK = 32;
      
      AmStrategyCUDA strat;
      // Data size
      // Non-power-of-two to check performance in weird situations
      int n = NUMBER_OF_BLOCKS*THREADS_PER_BLOCK*2+7;

      // Device arrays
      FieldCUDA d_a, d_b;
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
      cudaMemcpy(d_a.cDField(), a, n*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(d_b.cDField(), b, n*sizeof(cudaReal), cudaMemcpyHostToDevice);

      // Inner product on host
      cudaReal prodCheck = 0;
      for (int i = 0; i < n; i++) {
         prodCheck += a[i]*b[i];
      }

      // Inner product on device
      cudaReal prod = strat.innerProduct(d_a,d_b);

      TEST_ASSERT(prodCheck==prod);


   }
   

};

TEST_BEGIN(AmStrategyCUDATest)
// This test is currently disabled because the innerProduct member should be private.
TEST_ADD(AmStrategyCUDATest, testInnerProduct)
TEST_END(AmStrategyCUDATest)
#endif
