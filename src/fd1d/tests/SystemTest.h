#ifndef FD1D_SYSTEM_TEST_H
#define FD1D_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/System.h>
#include <fd1d/Mixture.h>
#include <fd1d/Domain.h>
#include <fd1d/Iterator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Fd1d;

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



   void testSolveMDE()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/System", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(Constants::Pi*(double(i)+0.5)/nx);
         sys.wField(0)[i] = 0.5 + cs;
         sys.wField(1)[i] = 0.7 - cs;
      }
      mix.compute(sys.wFields(), sys.cFields());

      // Test if same Q is obtained from different methods
      std::cout << mix.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(0, 1).computeQ() << "\n";

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

   void testIterator()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/System2", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      double chi = 20.0;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         sys.wField(0)[i] = chi*0.5*(1.0 - cs) - chi*0.5;
         sys.wField(1)[i] = chi*0.5*(1.0 + cs) - chi*0.5;
      }
      mix.compute(sys.wFields(), sys.cFields());

      sys.iterator().solve();

      std::ofstream out;
      openOutputFile("out/w", out);
      sys.writeFields(out, sys.wFields());
      out.close();

      openOutputFile("out/c", out);
      sys.writeFields(out, sys.cFields());
      out.close();

   }

   void testFieldInput()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/System2", in);

      System sys;
      sys.readParam(in);
      in.close();

      std::cout << "\n";
      // sys.writeParam(std::cout);

      openInputFile("in/omega", in);
      sys.readWFields(in);
      in.close();

      sys.iterator().solve();

      std::ofstream out;
      out.open("out/c2");
      sys.writeFields(out, sys.cFields());
      out.close();

   }

   void testReadCommands()
   {
      printMethod(TEST_FUNC);

      System sys;
      std::ifstream in;

      openInputFile("in/System2", in);
      sys.readParam(in);
      in.close();

      std::cout << "\n";
      openInputFile("in/command", in);
      sys.readCommands(in);
      in.close();
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testSolveMDE)
TEST_ADD(SystemTest, testIterator)
TEST_ADD(SystemTest, testFieldInput)
TEST_ADD(SystemTest, testReadCommands)
TEST_END(SystemTest)

#endif
