#ifndef AVERAGE_TEST_H
#define AVERAGE_TEST_H

#include <util/random/Random.h>
#include <util/random/Ar1Process.h>
#include <util/accumulators/Average.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class AverageTest 
{

public:

   AverageTest()
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


   void testAverage() 
   {
      setUp("testAverage");

      Average                  average;

      // Set parameters
      //average.setParam(1, 200, true);

      // Set parameters
      std::ifstream file("in/Average");
      average.readParam(file);

      Ar1Process process(random());
      double tau = 10.0;
      process.init(tau);

      long   nSample = 1000;
      double x;
      for (long i = 0; i < nSample; ++i) {
         x = process();
         std::cout << x << std::endl;
         average.sample(x);
      }
      std::cout << "Average = " << average.average() 
                <<  " +- " << average.blockingError() 
                << std::endl;
      average.output(std::cout);

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
   AverageTest test;

   test.testAverage();
}

#endif
