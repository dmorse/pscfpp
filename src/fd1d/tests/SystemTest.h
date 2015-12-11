#ifndef FD1D_SYSTEM_TEST_H
#define FD1D_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/System.h>
//#include <pfts/Block.h>
//#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pfts;
using namespace Pfts::Fd1d;

class SystemTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      System sys;
   }

   void testReadParameters()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/System", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);
   }

   void testSolve()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/System", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      double nx = (double)sys.nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] = cos(Constants::Pi*(double(i)+0.5)/nx);
         sys.wField(1)[i] = -cos(Constants::Pi*(double(i)+0.5)/nx);
      }
      sys.compute();

      // Test if same Q is obtained from different methods
      std::cout << sys.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << sys.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << sys.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << sys.polymer(0).propagator(0, 1).computeQ() << "\n";

      // Test spatial integral of block concentration
      double sum0 = 0.0;
      double sum1 = 0.0;
      for (int i = 0; i < nx; ++i) {
         sum0 += sys.cField(0)[i];
         sum1 += sys.cField(1)[i];
      }
      sum0 = sum0/double(nx);
      sum1 = sum1/double(nx);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      
   }
};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testSolve)
TEST_END(SystemTest)

#endif
