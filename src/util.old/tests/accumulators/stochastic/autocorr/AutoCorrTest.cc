#ifndef AUTO_CORR_TEST_H
#define AUTO_CORR_TEST_H

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/random/Random.h>
#include <util/random/Ar1Process.h>
#include <util/space/Vector.h>
#include <util/space/Tensor.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
//#include <util/containers/PackedData.h>
#include <util/accumulators/AutoCorr.h>        // template 
//#include <util/accumulators/AutoCorr_serial.h> // template 
#include <util/accumulators/AutoCorrArray.h>   // template 
#include <util/accumulators/Average.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

class AutoCorrTest 
{

public:

   AutoCorrTest()
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


   void testAutoCorrDouble() 
   {
      setUp("testAutoCorrDouble");

      AutoCorr<double, double> autocorr;
      Average                  average;
      autocorr.setParam(100);
      average.setNSamplePerBlock(1);

      Ar1Process process(random_);
      double tau = 10.0;
      process.init(tau);

      int    nSample = 1000000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = process();
         autocorr.sample(x);
         average.sample(x);
      }
      autocorr.output(std::cout);
      average.output(std::cout);
      tearDown();
   }

   void testAutoCorrVector() 
   {
      setUp("testAutoVector");

      AutoCorr<Vector, double> autocorr;
      autocorr.setParam(100);

      Ar1Process process0(random_);
      Ar1Process process1(random_);
      Ar1Process process2(random_);
      double tau = 10.0;
      process0.init(tau);
      process1.init(tau);
      process2.init(tau);

      int    nSample =  100000;
      Vector x;
      for (int i=0; i < nSample; i++) {
         x[0] = process0();
         x[1] = process1();
         x[2] = process2();
         autocorr.sample(x);
      }
      autocorr.output(std::cout);

      tearDown();
   }

   void testAutoCorrArrayDouble() 
   {
      setUp("testAutoCorrArrayDouble");

      int    nEnsemble = 5;
      int    nSample   = 50000;
      int    capacity  = 200;
      AutoCorrArray<double, double> autocorr;
      autocorr.setParam(nEnsemble, capacity);

      DArray<double>     x;
      x.allocate(nEnsemble);

      DArray<Ar1Process> processes;
      processes.allocate(nEnsemble);

      double tau = 20.0;
      int  i, j;

      for (int j = 0; j < nEnsemble; j++) {
         processes[j].setRNG(random());
         processes[j].init(tau);
      }
     
      for (i = 0; i < nSample; i++) {
         for (j = 0; j < nEnsemble; j++) {
            x[j] = processes[j]();;
         }
         autocorr.sample(x);
      }
      autocorr.output(std::cout);

      tearDown();
   }

   void testAutoCorrArrayVector() 
   {
      setUp("testAutoCorrArrayVector");

      int    nEnsemble = 5;
      int    nSample   = 50000;
      int    capacity  = 200;
      AutoCorrArray<Vector, double> autocorr;
      autocorr.setParam(nEnsemble, capacity);

      DArray<Vector> x;
      x.allocate(nEnsemble);

      DArray< FArray<Ar1Process, 3> > processes;
      processes.allocate(nEnsemble);

      double tau = 20.0;
      int  i, j;

      for (i = 0; i < nEnsemble; ++i) {
         for (j = 0; j < 3; ++j) {
            processes[i][j].setRNG(random());
            processes[i][j].init(tau);
         }
      }
    
      for (i = 0; i < nSample; ++i) {
         for (j = 0; j < nEnsemble; ++j) {
            x[j][0] = processes[j][0]();
            x[j][1] = processes[j][1]();
            x[j][2] = processes[j][2]();
         }
         autocorr.sample(x);
      }
      autocorr.output(std::cout);

      tearDown();
   }

   void testAutoCorrArrayTensor() 
   {
      setUp("testAutoCorrArrayTensor");

      int   nEnsemble = 5;
      int   nSample   = 100;
      int   capacity  = 4;
      AutoCorrArray<Tensor, double> autocorr;
      autocorr.setParam(nEnsemble, capacity);

      DArray<Tensor> x;
      Vector u;

      x.allocate(nEnsemble);

      double trace;
      int i, j, k;
      for (i = 0; i < nSample; ++i) {
         for (j = 0; j < nEnsemble; ++j) {
            random().unitVector(u);
            x[j].zero();
            x[j].dyad(u, u);
            trace = x[j].trace()/double(Dimension);
            for (k = 0; k < Dimension; ++k) {
               x[j](k, k) -= trace;
            }
         }
         autocorr.sample(x);
      }
      autocorr.output(std::cout);

      tearDown();
   }

   #if 0
   void testAutoCorrPack() 
   {
      setUp("testAutoCorrPack");

      AutoCorr<double, double> autocorr;
      Average average;
      autocorr.setParam(100);
      average.setNSamplePerBlock(1);

      Ar1Process process(random_);
      double tau = 10.0;
      process.init(tau);

      int    nSample = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = process();
         autocorr.sample(x);
         average.sample(x);
      }

      PackedData u;
      u.allocate(autocorr.packedSize());
      u.beginPacking();
      PackedData::Byte* current = u.begin();
      PackedData::Byte* end     = u.begin() + autocorr.packedSize();
      autocorr.pack(current, end);

      AutoCorr<double, double> clone;
      clone.setParam(100);

      u.beginUnpacking();
      current = u.begin();
      clone.unpack(current, end);

      clone.output(std::cout);
      tearDown();
   }
   #endif

   #if 0
   void testAutoCorrSerialize() 
   {
      setUp("testAutoCorrSerialize");

      AutoCorr<double, double> autocorr;
      Average average;
      autocorr.setParam(100);
      average.setNSamplePerBlock(1);

      Ar1Process process(random_);
      double tau = 10.0;
      process.init(tau);

      int    nSample = 100000;
      double x;
      for (int i=0; i < nSample; i++) {
         x = process();
         autocorr.sample(x);
         average.sample(x);
      }

      MemoryOArchive u;
      MemoryIArchive v;
      u.allocate(autocorr.packedSize());
      AutoCorr<double, double> clone;
      u << autocorr;
      v = u;
      v >> clone;

      clone.output(std::cout);
      tearDown();
   }
   #endif

   Random& random() 
   { return random_; }

private:

   Random random_;
   int    verbose_;

};


int main()
{
   AutoCorrTest test;

   test.testAutoCorrDouble();
   test.testAutoCorrVector();
   test.testAutoCorrArrayDouble();
   test.testAutoCorrArrayVector();
   test.testAutoCorrArrayTensor();
   //test.testAutoCorrPack();
   //test.testAutoCorrSerialize();
}

#endif
