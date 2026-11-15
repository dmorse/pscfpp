#ifndef PSCF_CUDA_RANDOM_TEST_H
#define PSCF_CUDA_RANDOM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cuda/CUT.h>
#include <pscf/backend/cuda/VecRandom.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/misc/Exception.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class CudaVecRandomTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      VecRandom<CUT> random;
      random.setSeed(6712983651284);
   }

   void testUniformDouble()
   {
      printMethod(TEST_FUNC);
      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      int n = 100000;
      BackendArray<float,CUT> d_(n);
      HostSendArray<float,CUT> h_(n);

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

      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      int n = 100000;
      BackendArray<float,CUT> df_(n);
      HostSendArray<float,CUT> hf_(n);

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

   void testUniformRangeDouble()
   {
      printMethod(TEST_FUNC);
      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      int n = 100000;
      BackendArray<float,CUT> d_(n);
      HostSendArray<float,CUT> h_(n);

      const double min = -1.34;
      const double max =  3.57;
      random.uniform(d_, min, max);

      h_ = d_;

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      const double mean = 0.5*(min + max);
      const double diff = max - min;
      const double stddev = diff/sqrt(12.0);

      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(h_[i] >= min);
         TEST_ASSERT(h_[i] < max);
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

   void testNormalDouble()
   {
      printMethod(TEST_FUNC);
      VecRandom<CUT> random;
      random.setSeed(6712983651284);
      
      int n = 100000;
      BackendArray<float,CUT> d_(n);
      HostSendArray<float,CUT> h_(n);

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
      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      int n = 100000;
      BackendArray<float,CUT> df_(n);
      HostSendArray<float,CUT> hf_(n);

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

   /*
   * Normal random numbers for arrays with an odd number of elements.
   *
   * Every element must be set, and the values of a large array must
   * have the requested mean and standard deviation.
   */
   template <typename Data>
   void checkNormalOdd()
   {
      VecRandom<CUT> random;
      random.setSeed(6712983651284);

      const Data mean = 1.0;
      const Data stddev = 0.5;
      const Data unset = 1.0E+30;

      int sizes[3] = {1, 3, 100001};
      for (int s = 0; s < 3; ++s) {
         int n = sizes[s];
         HostSendArray<Data,CUT> h_(n);
         for (int i = 0; i < n; ++i) {
            h_[i] = unset;
         }
         BackendArray<Data,CUT> d_(n);
         d_ = h_;

         bool thrown = false;
         try {
            random.normal(d_, stddev, mean);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(!thrown);
         h_ = d_;

         double ave = 0.0;
         double var = 0.0;
         double val = 0.0;
         for (int i = 0; i < n; ++i) {
            val = h_[i] - mean;
            TEST_ASSERT(fabs(val) < 10.0*stddev);
            ave += val;
            var += val*val;
         }
         if (n > 1000) {
            ave = ave/double(n);
            var = var/double(n);
            var = sqrt(var) - stddev;
            ave = ave/stddev;
            var = var/stddev;
            TEST_ASSERT(fabs(ave) < 0.1);
            TEST_ASSERT(fabs(var) < 0.1);
         }
      }
   }

   void testNormalOdd()
   {
      printMethod(TEST_FUNC);
      checkNormalOdd<float>();
      checkNormalOdd<double>();
   }

};

TEST_BEGIN(CudaVecRandomTest)
TEST_ADD(CudaVecRandomTest, testConstructor)
TEST_ADD(CudaVecRandomTest, testUniformDouble)
TEST_ADD(CudaVecRandomTest, testUniformFloat)
TEST_ADD(CudaVecRandomTest, testUniformRangeDouble)
TEST_ADD(CudaVecRandomTest, testNormalDouble)
TEST_ADD(CudaVecRandomTest, testNormalFloat)
TEST_ADD(CudaVecRandomTest, testNormalOdd)
TEST_END(CudaVecRandomTest)

#endif
