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
      openInputFile("in/SystemPlanar", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      TEST_ASSERT(sys.domain().mode() == Planar);
   }

   void testSolveMdePlanar()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/SystemPlanar", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*(double(i)+0.5)/double(nx-1));
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
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
   }


   void testSolveMdeSpherical()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/SystemSpherical", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      double nx = (double)domain.nx();
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = -cs;
         sys.wField(1)[i] = +cs;
      }
      mix.compute(sys.wFields(), sys.cFields());

      // Test if same Q is obtained from different methods
      std::cout << mix.polymer(0).propagator(0, 0).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(1, 0).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(1, 1).computeQ() << "\n";
      std::cout << mix.polymer(0).propagator(0, 1).computeQ() << "\n";

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
   }

   void testIteratorPlanar()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/SystemPlanar2", in);

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
         cs = cos(Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = chi*(-0.5*cs + 0.25*cs*cs);
         sys.wField(1)[i] = chi*(+0.5*cs + 0.25*cs*cs);
      }
      double shift = sys.wField(1)[nx-1];
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] -= shift;
         sys.wField(1)[i] -= shift;
      }

      // Compute initial state
      mix.compute(sys.wFields(), sys.cFields());

      std::ofstream out;
      openOutputFile("out/initialPlanar.w", out);
      sys.writeFields(out, sys.wFields());
      out.close();

      openOutputFile("out/initialPlanar.c", out);
      sys.writeFields(out, sys.cFields());
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalPlanar.w", out);
      sys.writeFields(out, sys.wFields());
      out.close();

      openOutputFile("out/finalPlanar.c", out);
      sys.writeFields(out, sys.cFields());
      out.close();

   }

   void testIteratorSpherical()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/SystemSpherical", in);

      System sys;
      sys.readParam(in);

      std::cout << "\n";
      sys.writeParam(std::cout);

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();

      // Create initial chemical potential fields
      double nx = (double)domain.nx();
      double cs;
      double chi = 80.0;
      for (int i = 0; i < nx; ++i) {
         cs = cos(Constants::Pi*double(i)/double(nx-1));
         sys.wField(0)[i] = -chi*cs/2.0;
         sys.wField(1)[i] = +chi*cs/2.0;
      }
      double shift = sys.wField(1)[nx-1];
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] -= shift;
         sys.wField(1)[i] -= shift;
      }

      // Solve MDE for initial fields
      mix.compute(sys.wFields(), sys.cFields());

      std::cout << "Average fraction 0 = " 
                << domain.spatialAverage(sys.cField(0)) << "\n";
      std::cout << "Average fraction 1 = " 
                << domain.spatialAverage(sys.cField(1)) << "\n";

      std::ofstream out;
      openOutputFile("out/initialSpherical.w", out);
      sys.writeFields(out, sys.wFields());
      out.close();

      openOutputFile("out/initialSpherical.c", out);
      sys.writeFields(out, sys.cFields());
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalSpherical.w", out);
      sys.writeFields(out, sys.wFields());
      out.close();

      openOutputFile("out/finalSpherical.c", out);
      sys.writeFields(out, sys.cFields());
      out.close();

   }

   void testFieldInput()
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/SystemPlanar2", in);

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
      std::cout << "\n";

      openInputFile("in/SystemPlanar2", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      std::cout << "Test file prefix = |" 
                << filePrefix() << "|" << std::endl;
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/command", in);
      sys.readCommands(in);
      in.close();
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testSolveMdePlanar)
TEST_ADD(SystemTest, testSolveMdeSpherical)
TEST_ADD(SystemTest, testIteratorPlanar)
TEST_ADD(SystemTest, testIteratorSpherical)
TEST_ADD(SystemTest, testFieldInput)
TEST_ADD(SystemTest, testReadCommands)
TEST_END(SystemTest)

#endif
