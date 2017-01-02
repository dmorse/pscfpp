#ifndef AVERAGE_STAGE_TEST_H
#define AVERAGE_STAGE_TEST_H

#include <util/random/Random.h>
#include <util/accumulators/AverageStage.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class AverageStageTest 
{

public:

   AverageStageTest()
    : randomPtr_(0),
      verbose_(2)
   {}

   void setUp(const char* functionName)
   { 
      std::cout << std::string(functionName) << " :" << std::endl << std::endl;

      randomPtr_ = new Random;

      std::ifstream file("in/Random");
      random().readParam(file);
   }

   void tearDown()
   { 
      delete randomPtr_;
      std::cout << "----------------------------------------------------" 
                << std::endl << std::endl;
   }


   void testAverageStage() 
   {
      setUp("testAverageStage");

      AverageStage average;

      double tau = 64.0;
      double A   = 1/tau;
      double B   = sqrt(2.0*A);

      int    nSample = 10000000;
      double x = random().gaussian();
      for (int i=0; i < nSample; i++) {
         x += -A*x + B*random().gaussian();
         average.sample(x);
      }
      std::cout << "Average  = " << average.average()  << std::endl;
      std::cout << "Variance = " << average.variance() << std::endl;
      //average.outputError(std::cout);

      tearDown();
   }

   Random& random() 
   { return *randomPtr_; }

private:

   Random* randomPtr_;
   int     verbose_;

};


int main()
{
   AverageStageTest test;

   test.testAverageStage();
}

#endif
