#ifndef FD1D_SYSTEM_TEST_H
#define FD1D_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/System.h>
#include <fd1d/domain/Domain.h>
#include <fd1d/solvers/Mixture.h>
#include <fd1d/iterator/Iterator.h>
#include <fd1d/misc/FieldIo.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Fd1d;

class SystemTest : public UnitTest 
{

private:

   std::ofstream logFile_;

public:

   void setUp()
   {}

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
      ParamComponent::setEcho(false);
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
      ParamComponent::setEcho(true);
   }

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      System sys;
   }

   void testReadParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadParameters.log");

      std::ifstream in;
      openInputFile("in/planar1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      TEST_ASSERT(sys.domain().mode() == Planar);
   }

   void testSolveMdePlanar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSolveMdePlanar.log");

      std::ifstream in;
      openInputFile("in/planar1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

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
      Log::file() << mix.polymer(0).propagator(0, 0).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(1, 0).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(1, 1).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(0, 1).computeQ() << "\n";

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";

      TEST_ASSERT(eq(mix.polymer(0).length(), 5.0));
   }


   void testSolveMdeSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSolveMdeSpherical.log");

      std::ifstream in;
      openInputFile("in/spherical1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

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
      Log::file() << mix.polymer(0).propagator(0, 0).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(1, 0).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(1, 1).computeQ() << "\n";
      Log::file() << mix.polymer(0).propagator(0, 1).computeQ() << "\n";

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
   }

   void testIteratorPlanar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanar.log");

      std::ifstream in;
      openInputFile("in/planar2.prm", in);

      System sys;
      sys.readParam(in);
      FieldIo fieldIo(sys);

      Log::file() << "\n";
      sys.writeParam(Log::file());

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
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialPlanar.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalPlanar.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanar.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   void testIteratorSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorSpherical.log");

      std::ifstream in;
      openInputFile("in/spherical1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();
      FieldIo fieldIo(sys);

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

      Log::file() << "Average fraction 0 = " 
                << domain.spatialAverage(sys.cField(0)) << "\n";
      Log::file() << "Average fraction 1 = " 
                << domain.spatialAverage(sys.cField(1)) << "\n";

      std::ofstream out;
      openOutputFile("out/initialSpherical.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialSpherical.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalSpherical.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalSpherical.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   void testFieldInput()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestFieldInput.log");

      std::ifstream in;
      openInputFile("in/planar2.prm", in);

      System sys;
      sys.readParam(in);
      in.close();

      FieldIo fieldIo(sys);

      Log::file() << "\n";
      // sys.writeParam(Log::file());

      openInputFile("in/planar.w", in);
      fieldIo.readFields(sys.wFields(), in);
      in.close();

      sys.iterator().solve();

      std::ofstream out;
      out.open("out/c2");
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   void testReadCommandsPlanar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadCommandsPlanar.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/planar2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      Log::file() << "Test file prefix = |" 
                << filePrefix() << "|" << std::endl;
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/planar.cmd", in);
      sys.readCommands(in);
      in.close();
   }

   void testReadCommandsSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadCommandsSpherical.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      Log::file() << "Test file prefix = |" 
                << filePrefix() << "|" << std::endl;
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/spherical2.cmd", in);
      sys.readCommands(in);
      in.close();
   }

   void testReadCommandsSphericalSweep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadCommandsSphericalSweep.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical3.prm", in);
      sys.readParam(in);
      in.close();
      Log::file() << "Finished reading param file" << std::endl;
      sys.writeParam(Log::file());

      // Set System filemaster prefixes to unit test file prefix
      Log::file() << "Test file prefix = |" 
                  << filePrefix() << "|" << std::endl;
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/spherical3.cmd", in);
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
TEST_ADD(SystemTest, testReadCommandsPlanar)
TEST_ADD(SystemTest, testReadCommandsSpherical)
TEST_ADD(SystemTest, testReadCommandsSphericalSweep)
TEST_END(SystemTest)

#endif
