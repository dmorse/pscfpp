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
   {}

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
      int n = 400;
      allocate(n);
      random.uniform(devPtr_, n);

      copyDevToHost();

      std::cout << std::endl;
      for (int i = 0; i < n; ++i) {
         std::cout << Int(i,5) << "  " 
                   << Dbl(hostPtr_[i]) << std::endl;
      }
   }

};

TEST_BEGIN(CudaRandomTest)
TEST_ADD(CudaRandomTest, testConstructor)
TEST_ADD(CudaRandomTest, testUniform)
TEST_END(CudaRandomTest)

#endif
