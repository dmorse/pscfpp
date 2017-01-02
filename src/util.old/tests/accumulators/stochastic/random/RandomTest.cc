#ifndef RANDOM_TEST_H
#define RANDOM_TEST_H

#include <util/random/Random.h>
#include <util/accumulators/Average.h>
#include <util/accumulators/Distribution.h>
#include <util/accumulators/IntDistribution.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class RandomTest 
{

public:

   RandomTest()
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


   void testReadParam() 
   {
      setUp("testReadParam");

      // Verbose output
      if (verbose_ > 0) {
         printf("idum: %ld\n", random().seed() );
      }

      tearDown();
   }

   void testGetFloat() 
   {
      setUp("testGetFloat");

      Average average;
      int nSamplePerBlock   = 1;
      average.setNSamplePerBlock(nSamplePerBlock);

      Distribution distribution;
      double min = 0.0;
      double max = 1.0;
      int nBin = 10;
      distribution.setParam(min, max, nBin);

      const int nSample  = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = random().uniform(0.0,1.0);
         average.sample(x);
         distribution.sample(x);
      }

      average.output(std::cout);
      distribution.output(std::cout);

      tearDown();
   }

   void testGetInteger() 
   {
      setUp("testGetInteger");

      const int nSample  = 100000;

      Average average;
      int nSamplePerBlock   = 1;
      average.setNSamplePerBlock(nSamplePerBlock);


      int  nBin   = 11;
      IntDistribution distribution;
      int  min = 0;
      int  max = min + nBin - 1;
      distribution.setParam(min, max);

      long x;
      for (int i=0; i < nSample; i++) {
         x = random().uniformInt(0, nBin);
         average.sample(double(x));
         distribution.sample(x);
      }

      average.output(std::cout);
      distribution.output(std::cout);

      tearDown();
   }


   void testGaussian() 
   {
      setUp("testGaussian");

      Average average;
      int nSamplePerBlock   = 1;
      average.setNSamplePerBlock(nSamplePerBlock);

      Distribution distribution;
      double min = -3.0;
      double max =  3.0;
      int nBin   =  60;
      distribution.setParam(min, max, nBin);

      const int nSample  = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = random().gaussian();
         average.sample(x);
         distribution.sample(x);
      }

      average.output(std::cout);
      distribution.output(std::cout);

      tearDown();
   }


   void testUnitVector() 
   {
      setUp("testUnitVector");

      Average average;
      int nSamplePerBlock   = 1;
      average.setNSamplePerBlock(nSamplePerBlock);

      Distribution distribution;
      double min = -1.1;
      double max =  1.1;
      int nBin   =  22;
      distribution.setParam(min, max, nBin);

      Vector v;
      double vsq;
      int nSample = 100000;
      for (int i=0; i < nSample; i++ ) {
         random().unitVector(v);
         vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
         average.sample(v[1]);
         distribution.sample(v[1]);
      }
      average.output(std::cout);
      distribution.output(std::cout);

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
   RandomTest test;

   test.testReadParam();
   test.testGetFloat();
   test.testGetInteger();
   test.testGaussian();
   test.testUnitVector();
}

#endif
