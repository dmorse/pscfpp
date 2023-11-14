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

   cudaReal* devPtr_;
   cudaReal* hostPtr_;
   int n_;

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}
  
   void allocate(int n)
   {
      n_ = n;
      gpuErrchk(cudaMalloc((void**) &devPtr_, n * sizeof(cudaReal)));
      hostPtr_ = new cudaReal[n];
   }
  
   void copyDevToHost()
   {
      cudaMemcpy((void *)hostPtr_, (void const *) devPtr_, 
                 n_*sizeof(cudaReal), cudaMemcpyDeviceToHost);
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
      random.uniform(devPtr_, n);

      copyDevToHost();

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
         TEST_ASSERT(hostPtr_[i] > 0.0);
         TEST_ASSERT(hostPtr_[i] <= 1.0);
         val = hostPtr_[i] - mean;
         ave += val;
         var += val*val;
         if (verbose() > 1) {
            std::cout << Int(i,5) << "  " 
                      << Dbl(hostPtr_[i]) << std::endl;
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
      random.normal(devPtr_, n, stddev, mean);

      copyDevToHost();

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      cudaReal ave = 0.0;
      cudaReal var = 0.0;
      cudaReal val = 0.0;
      for (int i = 0; i < n; ++i) {
         val = hostPtr_[i] - mean;
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
