#ifndef SYMM_TENSOR_AVERAGE_TEST_H
#define SYMM_TENSOR_AVERAGE_TEST_H

#include <util/random/Random.h>
#include <util/random/Ar1Process.h>
#include <util/accumulators/SymmTensorAverage.h>
#include <util/space/Tensor.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class SymmTensorAverageTest 
{

public:

   SymmTensorAverageTest()
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

      SymmTensorAverage accumulator;

      // Set parameters
      //accumulator.setParam(1, 200, true);

      // Set parameters
      std::ifstream file("in/SymmTensorAverage");
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
            for (j = 0; j <= i; ++j) {
               value(i, j) = processes[k]();
               if (j != i) {
                  value(j, i) = value(i, j);
               }
               //std::cout << value(i, j);
               ++k;
            }
         }
         //std::cout << "\n";
         accumulator.sample(value);
      }

      double ave;
      double err;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j <= i; ++j) {
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
   SymmTensorAverageTest test;

   test.testAverage();
}

#endif
