#ifndef PSCF_CPU_VEC_RANDOM_TEST_H
#define PSCF_CPU_VEC_RANDOM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/CpuVecRandom.h>
#include <util/random/Random.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class CpuVecRandomTest : public UnitTest 
{
private:

   Util::Random sRandom_;
   Pscf::CpuVecRandom vRandom_;
   DArray<double> data_;

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {}

   void initialize(int n, long shift = 0)
   {
      long seed = 6712983651284;
      seed += shift;
      TEST_ASSERT(seed >= 0);
      sRandom_.setSeed(seed);
      vRandom_.associate(sRandom_);
      data_.allocate(n);
   }

   void testConstructor()
   {  printMethod(TEST_FUNC); }

   void testUniform()
   {
      printMethod(TEST_FUNC);

      int n = 100000;
      initialize(n, 4);

      vRandom_.uniform(data_);

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      const double mean = 0.5;
      const double stddev = 1.0/sqrt(12.0);
      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(data_[i] >= 0.0);
         TEST_ASSERT(data_[i] < 1.0);
         val = data_[i] - mean;
         ave += val;
         var += val*val;
         if (verbose() > 1) {
            std::cout << Int(i,5) << "  " 
                      << Dbl(data_[i]) << std::endl;
         }
      }
      ave = ave/double(n);
      var = var/double(n);
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

   void testUniformRange()
   {
      printMethod(TEST_FUNC);

      int n = 100000;
      initialize(n, -30);

      const double min = -1.3;
      const double max = +3.0;
      vRandom_.uniform(data_, min, max);

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      const double mean = 0.5 * (min + max);
      const double diff = max - min;
      const double stddev = diff/sqrt(12.0);
      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
         TEST_ASSERT(data_[i] >= min);
         TEST_ASSERT(data_[i] < max);
         val = data_[i] - mean;
         ave += val;
         var += val*val;
         if (verbose() > 1) {
            std::cout << Int(i,5) << "  " 
                      << Dbl(data_[i]) << std::endl;
         }
      }
      ave = ave/double(n);
      var = var/double(n);
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
      
      int n = 100000;
      initialize(n, 23);

      double stddev = 0.5; 
      double mean = 1.0;
      vRandom_.normal(data_, stddev, mean);

      //setVerbose(1);
      if (verbose() > 0) {
         std::cout << std::endl;
      }

      double ave = 0.0;
      double var = 0.0;
      double val = 0.0;
      for (int i = 0; i < n; ++i) {
         val = data_[i] - mean;
         ave += val;
         var += val*val;
      }
      ave = ave/double(n);
      var = var/double(n);
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

TEST_BEGIN(CpuVecRandomTest)
TEST_ADD(CpuVecRandomTest, testConstructor)
TEST_ADD(CpuVecRandomTest, testUniform)
TEST_ADD(CpuVecRandomTest, testUniformRange)
TEST_ADD(CpuVecRandomTest, testNormal)
TEST_END(CpuVecRandomTest)

#endif
