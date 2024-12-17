#ifndef PSCF_CUDA_RANDOM_TEST_H
#define PSCF_CUDA_RANDOM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/CudaRandom.h>
#include <pscf/cuda/GpuTypes.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class CudaRandomTest : public UnitTest 
{

public:

   DeviceArray<cudaReal> d_;
   HostDArray<cudaReal> h_;
   int n_;

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}
  
   void allocate(int n)
   {
      n_ = n;
      d_.allocate(n);
      h_.allocate(n);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);
   }

   void testUniform()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);
      int n = 100000;
      allocate(n);
      random.uniform(d_.cArray(), n);

      h_ = d_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      cudaReal mean = 0.5;
      cudaReal stddev = 1.0/sqrt(12.0);
      cudaReal ave = 0.0;
      cudaReal var = 0.0;
      cudaReal val = 0.0;
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(h_[i] > 0.0);
         TEST_ASSERT(h_[i] <= 1.0);
         val = h_[i] - mean;
         ave += val;
         var += val*val;
         if (verbose() > 1) {
            std::cout << Int(i,5) << "  " 
                      << Dbl(h_[i]) << std::endl;
         }
      }
      ave = ave/cudaReal(n);
      var = var/cudaReal(n);
      ave = ave;
      var = sqrt(var) - stddev;
      ave = ave/stddev;
      var = var/stddev;
      if (verbose() > 0) {
         std::cout << "Average  " << ave << std::endl;
         std::cout << "StdDev   " << var << std::endl;
      }
      TEST_ASSERT(fabs(ave) < 0.1);
      TEST_ASSERT(fabs(var) < 0.1);

   }

   void testNormal()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);
      int n = 100000;
      allocate(n);

      cudaReal mean = 1.0;
      cudaReal stddev = 0.5; 
      random.normal(d_.cArray(), n, stddev, mean);

      h_ = d_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      cudaReal ave = 0.0;
      cudaReal var = 0.0;
      cudaReal val = 0.0;
      for (int i = 0; i < n; ++i) {
         val = h_[i] - mean;
         ave += val;
         var += val*val;
      }
      ave = ave/cudaReal(n);
      var = var/cudaReal(n);
      ave = ave;
      var = sqrt(var) - stddev;
      ave = ave/stddev;
      var = var/stddev;
      if (verbose() > 0) {
         std::cout << "Average  " << ave << std::endl;
         std::cout << "StdDev   " << var << std::endl;
      }
      TEST_ASSERT(fabs(ave) < 0.1);
      TEST_ASSERT(fabs(var) < 0.1);
   }
};

TEST_BEGIN(CudaRandomTest)
TEST_ADD(CudaRandomTest, testConstructor)
TEST_ADD(CudaRandomTest, testUniform)
TEST_ADD(CudaRandomTest, testNormal)
TEST_END(CudaRandomTest)

#endif
