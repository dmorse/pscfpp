#ifndef MEAN_SQ_DISP_TEST_H
#define MEAN_SQ_DISP_TEST_H

#include <util/random/Random.h>
#include <util/space/Vector.h>
#include <util/containers/DArray.h>
#include <util/accumulators/MeanSqDispArray.h>  // template implementation
#include <util/accumulators/Average.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class MeanSqDispTest 
{

public:

   MeanSqDispTest()
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

   void testMeanSqDispArrayDouble() 
   {
      setUp("testMeanSqDispArrayDouble");

      int    nEnsemble = 1000;
      int    nSample   = 1000;
      int    capacity  = 201;
      MeanSqDispArray<double> msd;
      msd.setParam(nEnsemble, capacity);

      DArray<double> x;
      x.allocate(nEnsemble);

      double D   = 0.10;
      double B   = sqrt(2.0*D);

      int    i, j;
      for (j = 0; j < nEnsemble; ++j) {
         x[j] = random().gaussian();
      }
      for (i = 0; i < nSample; ++i) {
         for (j = 0; j < nEnsemble; ++j) {
            x[j] += B*random().gaussian();
         }
         msd.sample(x);
      }

      msd.output(std::cout);

      tearDown();
   }

   void testMeanSqDispArrayVector() 
   {
      setUp("testMeanSqDispArrayVector");

      int    nEnsemble = 1000;
      int    nSample   = 1000;
      int    capacity  = 201;
      MeanSqDispArray<Vector> msd;
      msd.setParam(nEnsemble, capacity);

      DArray<Vector> x;
      x.allocate(nEnsemble);

      double D   = 0.1;
      double B   = sqrt(2.0*D);

      int    i, j;
      for (j = 0; j < nEnsemble; ++j) {
         x[j][0] = random().gaussian();
         x[j][1] = random().gaussian();
         x[j][2] = random().gaussian();
      }
      for (i = 0; i < nSample; ++i) {
         for (j = 0; j < nEnsemble; ++j) {
            //x[j][0] += -A*x[j][0] + B*random().gaussian();
            //x[j][1] += -A*x[j][1] + B*random().gaussian();
            //x[j][2] += -A*x[j][2] + B*random().gaussian();
            x[j][0] += B*random().gaussian();
            x[j][1] += B*random().gaussian();
            x[j][2] += B*random().gaussian();
         }
         msd.sample(x);
      }

      msd.output(std::cout);

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
   MeanSqDispTest test;

   test.testMeanSqDispArrayDouble();
   test.testMeanSqDispArrayVector();
}

#endif
