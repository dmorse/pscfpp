#ifndef PRDC_CUDA_MEM_TEST_H
#define PRDC_CUDA_MEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/RField.h>
#include <pscf/cuda/GpuResources.h>
#include <util/math/Constants.h>

#include <fstream>
#include <iomanip>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaMemTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testCopyRoundTrip()
   {
      printMethod(TEST_FUNC);
      
      int nx = 10;
      // device array
      Cuda::RField<1> d_in;
      d_in.allocate(10);
      // host arrays
      cudaReal* in = new cudaReal[nx];
      cudaReal* out = new cudaReal[nx];

      // random data
      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         in[i] = cos(twoPi*double(i)/double(nx));
      }

      // copy round trip
      cudaMemcpy(d_in.cField(), in, nx*sizeof(cudaReal), cudaMemcpyHostToDevice);
      cudaMemcpy(out, d_in.cField(), nx*sizeof(cudaReal), cudaMemcpyDeviceToHost);

      // compare
      double maxDiff = 0, currDiff = 0;
      for (int i = 0; i < nx; ++i ) {
         currDiff = std::abs( in[i] - out[i] );
         if ( currDiff > maxDiff ) {
            maxDiff = currDiff;
         }
         TEST_ASSERT( in[i] == out[i] );
      }
      // std::cout << std::setprecision(16) << maxDiff << std::endl;
   }

};

TEST_BEGIN(CudaMemTest)
TEST_ADD(CudaMemTest, testCopyRoundTrip)
TEST_END(CudaMemTest)

#endif
