#ifndef PSCF_CUDA_RANDOM_TEST_H
#define PSCF_CUDA_RANDOM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/CudaRandom.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class CudaRandomTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);
   }

   void testUniformDouble()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);

      int n = 100000;
      DeviceArray<float> d_(n);
      HostDArray<float> h_(n);

      random.uniform(d_);

      h_ = d_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      double mean = 0.5;
      double stddev = 1.0/sqrt(12.0);
      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
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
      ave = ave/double(n);
      var = var/double(n);
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

   void testUniformFloat()
   {
      printMethod(TEST_FUNC);

      CudaRandom random;
      random.setSeed(6712983651284);

      int n = 100000;
      DeviceArray<float> df_(n);
      HostDArray<float> hf_(n);

      random.uniform(df_);

      hf_ = df_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      float mean = 0.5;
      float stddev = 1.0/sqrt(12.0);
      float ave = 0.0;
      float var = 0.0;
      float val = 0.0;
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(hf_[i] > 0.0);
         TEST_ASSERT(hf_[i] <= 1.0);
         val = hf_[i] - mean;
         ave += val;
         var += val*val;
         if (verbose() > 1) {
            std::cout << Int(i,5) << "  " 
                      << Dbl(hf_[i]) << std::endl;
         }
      }
      ave = ave/float(n);
      var = var/float(n);
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

   void testNormalDouble()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);
      
      int n = 100000;
      DeviceArray<float> d_(n);
      HostDArray<float> h_(n);

      double mean = 1.0;
      double stddev = 0.5; 
      random.normal(d_, stddev, mean);

      h_ = d_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
         val = h_[i] - mean;
         ave += val;
         var += val*val;
      }
      ave = ave/double(n);
      var = var/double(n);
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

   void testNormalFloat()
   {
      printMethod(TEST_FUNC);
      CudaRandom random;
      random.setSeed(6712983651284);

      int n = 100000;
      DeviceArray<float> df_(n);
      HostDArray<float> hf_(n);

      float mean = 1.0;
      float stddev = 0.5; 
      random.normal(df_, stddev, mean);

      hf_ = df_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      float ave = 0.0;
      float var = 0.0;
      float val = 0.0;
      for (int i = 0; i < n; ++i) {
         val = hf_[i] - mean;
         ave += val;
         var += val*val;
      }
      ave = ave/float(n);
      var = var/float(n);
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
TEST_ADD(CudaRandomTest, testUniformDouble)
TEST_ADD(CudaRandomTest, testUniformFloat)
TEST_ADD(CudaRandomTest, testNormalDouble)
TEST_ADD(CudaRandomTest, testNormalFloat)
TEST_END(CudaRandomTest)

#endif
