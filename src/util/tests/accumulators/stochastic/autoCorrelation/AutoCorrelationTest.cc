#ifndef AUTO_CORRELATION_TEST_H
#define AUTO_CORRELATION_TEST_H

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/random/Random.h>
#include <util/random/Ar1Process.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/accumulators/AutoCorrelation.tpp>     // template 
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class AutoCorrelationTest 
{

public:

   AutoCorrelationTest()
    : verbose_(2)
   {}

   void setUp(const char* functionName)
   { 
      std::cout << std::string(functionName) << " :" << std::endl << std::endl;
      std::ifstream file("in/Random");
      random().readParam(file);
   }

   void tearDown()
   { 
      std::cout << "----------------------------------------------------" 
                << std::endl << std::endl;
   }


   void testAutoCorrDouble(int maxStageId) 
   {
      setUp("testAutoCorrDouble");

      AutoCorrelation<double, double> autocorr;
      autocorr.setParam(64, maxStageId);

      Ar1Process process(random_);
      double tau = 10.0;
      process.init(tau);

      int nSample = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = process();
         autocorr.sample(x);
      }
      autocorr.output(std::cout);

      std::cout << "maxDelay = " 
                << autocorr.maxDelay()
                << std::endl;

      tearDown();
   }

   void testAutoCorrVector(int maxStageId) 
   {
      setUp("testAutoVector");

      AutoCorrelation<Vector, double> autocorr;
      autocorr.setParam(64, maxStageId);

      Ar1Process process0(random_);
      Ar1Process process1(random_);
      Ar1Process process2(random_);
      double tau = 10.0;
      process0.init(tau);
      process1.init(tau);
      process2.init(tau);

      int    nSample = 100000;
      Vector x;
      for (int i=0; i < nSample; i++) {
         x[0] = process0();
         x[1] = process1();
         x[2] = process2();
         autocorr.sample(x);
      }
      autocorr.output(std::cout);
      std::cout << "maxDelay = " 
                << autocorr.maxDelay()
                << std::endl;

      tearDown();
   }

   void testAutoCorrSerialize(int maxStageId) 
   {
      setUp("testAutoCorrSerialize");

      AutoCorrelation<double, double> autocorr;
      autocorr.setParam(64, maxStageId);

      Ar1Process process(random_);
      double tau = 10.0;
      process.init(tau);

      int nSample = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = process();
         autocorr.sample(x);
      }

      autocorr.output(std::cout);
      std::cout << "End Original ----------------------------------------" 
                << std::endl << std::endl;

      // Save to binary file
      BinaryFileOArchive u("binary");
      u << autocorr;
      u.file().close();

      // Load from binary file
      BinaryFileIArchive v("binary");
      AutoCorrelation<double, double> autocorr2;
      v >> autocorr2;

      autocorr2.output(std::cout);
      std::cout << "End Clone -------------------------------------------" 
                << std::endl << std::endl;

      std::cout << "maxDelay = " 
                << autocorr2.maxDelay()
                << std::endl;

      tearDown();
   }

   Random& random() 
   { return random_; }

private:

   Random random_;
   int    verbose_;

};

int main()
{
   AutoCorrelationTest test;

   test.testAutoCorrDouble(0);
   test.testAutoCorrDouble(3);
   test.testAutoCorrSerialize(2);
   test.testAutoCorrVector(3);
}

#endif
