#ifndef TENSOR_AVERAGE_TEST_H
#define TENSOR_AVERAGE_TEST_H

#include <util/random/Random.h>
#include <util/random/Ar1Process.h>
#include <util/accumulators/TensorAverage.h>
#include <util/space/Tensor.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class TensorAverageTest 
{

public:

   TensorAverageTest()
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

      TensorAverage accumulator;

      // Set parameters
      //accumulator.setParam(1, 200, true);

      // Set parameters
      std::ifstream file("in/TensorAverage");
      accumulator.readParam(file);

      // Create array of auto-regressive AR1 processes
      const int N = 9;
      double tau = 10.0;
      FArray<Ar1Process, N> processes;
      for (int m = 0; m < N; ++m) {
         processes[m].setRNG(random());
         processes[m].init(tau);
      }

      int nSample = 1000;
      Tensor value;
      int i, j, k;
      for (int m = 0; m < nSample; ++m) {
         k = 0;
         for (i = 0; i < Dimension; ++i) {
            for (j = 0; j < Dimension; ++j) {
               value(i, j) = processes[k]();
               ++k;
            }
         }
         //std::cout << x << std::endl;
         accumulator.sample(value);
      }

      double ave;
      double err;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            ave = accumulator(i, j).average();
            err = accumulator(i, j).blockingError();
            std::cout << "average(" << i << ", " << j << ") = "
                      << ave <<  " +- " << err << "\n";
            accumulator(i, j).output(std::cout);
         }
      }

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
   TensorAverageTest test;

   test.testAverage();
}

#endif
