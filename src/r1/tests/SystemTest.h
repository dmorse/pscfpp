#ifndef R1_SYSTEM_TEST_H
#define R1_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <r1/system/System.h>
#include <r1/field/Domain.h>
#include <r1/solvers/Mixture.h>
#include <r1/solvers/Polymer.h>
#include <r1/iterator/Iterator.h>
#include <r1/field/FieldIo.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::R1;

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
      setVerbose(0);
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
      openInputFile("in/planar_nr1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      TEST_ASSERT(sys.domain().mode() == Planar);
   }

   /*
   * NrIterator for a homopolymer A and a solvent B, with a NaN in the
   * field of B at one grid point, which makes the concentration of B
   * NaN. A residual that is not finite must not be accepted as
   * converged.
   */
   void testIteratorPlanarNrNan()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNrNan.log");

      std::istringstream in(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0"
         " nPolymer 1 nSolvent 1"
         " Polymer{ type linear nBlock 1 blocks 0 1.0 phi 0.5 }"
         " Solvent{ monomerId 1 size 1.0 phi 0.5 }"
         " ds 0.01 }"
         " Interaction{ chi( 0 1 1.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 1.0 nx 11 }"
         " NrIterator{ epsilon 1.0e-8 } }\n");
      System sys;
      sys.readParam(in);
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] = 0.0;
         sys.wField(1)[i] = 0.0;
      }
      sys.wField(1)[nx/2] = std::nan("");
      TEST_ASSERT(sys.iterator().solve() == 1);
   }

   /*
   * COMPUTE for an open homopolymer B and an open solvent A (mu = 0),
   * in fields for which phi of the homopolymer (w_B = -800) or of the
   * solvent (w_A = -800) is above the range of a double: the command
   * must stop with an Exception. For fields w = 0, it must not.
   */
   void testComputeOutOfRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestComputeOutOfRange.log");

      double wA[3] = {0.0, -800.0, 0.0};
      double wB[3] = {-800.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k) {
         std::istringstream in(
            "System{ Mixture{ nMonomer 2 monomers 1.0 1.0"
            " nPolymer 1 nSolvent 1"
            " Polymer{ type linear nBlock 1 blocks 1 1.0 mu 0.0 }"
            " Solvent{ monomerId 0 size 1.0 mu 0.0 }"
            " ds 0.01 }"
            " Interaction{ chi( 0 1 1.0 ) }"
            " Domain{ mode Planar xMin 0.0 xMax 1.0 nx 11 }"
            " NrIterator{ epsilon 1.0e-8 } }\n");
         System sys;
         sys.readParam(in);
         int nx = sys.domain().nx();
         for (int i = 0; i < nx; ++i) {
            sys.wField(0)[i] = wA[k];
            sys.wField(1)[i] = wB[k];
         }

         std::istringstream commands("COMPUTE\nFINISH\n");
         bool thrown = false;
         try {
            sys.readCommands(commands);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown == (k < 2));
      }
   }

   void testSolveMdePlanar()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSolveMdePlanar.log");

      std::ifstream in;
      openInputFile("in/planar_nr1.prm", in);

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
      double q00, q01, q10, q11;
      mix.polymer(0).propagator(0, 0).computeQ(q00);
      mix.polymer(0).propagator(0, 1).computeQ(q01);
      mix.polymer(0).propagator(1, 0).computeQ(q10);
      mix.polymer(0).propagator(1, 1).computeQ(q11);

      setVerbose(1);
      if (verbose() > 0) {
         Log::file() << q00 << "\n";
         Log::file() << q01 << "\n";
         Log::file() << q10 << "\n";
         Log::file() << q11 << "\n";
      }
      UTIL_ASSERT(abs((q01 - q00)/q00) < 1.0E-6);
      UTIL_ASSERT(abs((q10 - q00)/q00) < 1.0E-6);
      UTIL_ASSERT(abs((q11 - q00)/q00) < 1.0E-6);

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
      openInputFile("in/spherical_nr1.prm", in);

      System sys;
      sys.readParam(in);

      TEST_ASSERT( !sys.domain().isShell() );

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
      double q00, q01, q10, q11;
      mix.polymer(0).propagator(0, 0).computeQ(q00);
      mix.polymer(0).propagator(0, 1).computeQ(q01);
      mix.polymer(0).propagator(1, 0).computeQ(q10);
      mix.polymer(0).propagator(1, 1).computeQ(q11);
      if (verbose() > 1) {
         Log::file() << q00 << "\n";
         Log::file() << q01 << "\n";
         Log::file() << q10 << "\n";
         Log::file() << q11 << "\n";
      }
      UTIL_ASSERT(abs((q01 - q00)/q00) < 1.0E-6);
      UTIL_ASSERT(abs((q10 - q00)/q00) < 1.0E-6);
      UTIL_ASSERT(abs((q11 - q00)/q00) < 1.0E-6);

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
   }

   /*
   * Test NR iterator on input field generated within function.
   */
   void testIteratorPlanarNr1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr1.log");

      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);

      System sys;
      sys.readParam(in);
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

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
      openOutputFile("out/initialPlanarNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialPlanarNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalPlanarNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanarNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   /*
   * Test NR iterator on input w field read from file.
   */
   void testIteratorPlanarNr2()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr2.log");

      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);

      System sys;
      sys.readParam(in);
      in.close();

      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

      Log::file() << "\n";
      // sys.writeParam(Log::file());

      openInputFile("in/planar.w", in);
      fieldIo.readFields(sys.wFields(), in);
      in.close();

      sys.iterator().solve();

      std::ofstream out;
      out.open("out/planarNr2.c");
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   /*
   * Test NR iterator controlled by a command file.
   */
   void testIteratorPlanarNr3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarNr3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/planar_nr.cmd", in);
      sys.readCommands(in);
      in.close();
   }

   /*
   * Test NR iterator on input field generated within function.
   */
   void testIteratorSphericalNr1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorSphericalNr1.log");

      std::ifstream in;
      openInputFile("in/spherical_nr1.prm", in);

      System sys;
      sys.readParam(in);

      Log::file() << "\n";
      sys.writeParam(Log::file());

      Mixture& mix = sys.mixture();
      Domain& domain = sys.domain();
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

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
      openOutputFile("out/initialSphericalNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialSphericalNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalSphericalNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalSphericalNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   void testIteratorSphericalNr3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorSphericalNr3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical_nr2.prm", in);
      sys.readParam(in);
      in.close();

      TEST_ASSERT(!sys.domain().isShell());

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/spherical_nr.cmd", in);
      sys.readCommands(in);
      in.close();
   }

   void testIteratorPlanarAm1()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarAm.log");

      std::ifstream in;
      openInputFile("in/planar_am2.prm", in);

      System sys;
      sys.readParam(in);
      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());

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
      openOutputFile("out/initialPlanarAm.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/initialPlanarAm.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

      sys.iterator().solve();

      openOutputFile("out/finalPlanarAm.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanarAm.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   /*
   * Test AM iterator controlled by a command file.
   */
   void testIteratorPlanarAm3()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarAm3.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/planar_am2.prm", in);
      sys.readParam(in);
      in.close();

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/planar_am.cmd", in);
      sys.readCommands(in);
      in.close();
   }

   void testSweepSpherical()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSweepSpherical.log");

      System sys;
      std::ifstream in;
      Log::file() << "\n";

      openInputFile("in/spherical3_nr.prm", in);
      sys.readParam(in);
      in.close();
      sys.writeParam(Log::file());

      TEST_ASSERT( !sys.domain().isShell() );

      // Set System filemaster prefixes to unit test file prefix
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      openInputFile("in/sphericalSweep.cmd", in);
      sys.readCommands(in);
      in.close();
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testIteratorPlanarNrNan)
TEST_ADD(SystemTest, testComputeOutOfRange)
TEST_ADD(SystemTest, testSolveMdePlanar)
TEST_ADD(SystemTest, testSolveMdeSpherical)
TEST_ADD(SystemTest, testIteratorPlanarNr1)
TEST_ADD(SystemTest, testIteratorPlanarNr2)
TEST_ADD(SystemTest, testIteratorPlanarNr3)
TEST_ADD(SystemTest, testIteratorSphericalNr1)
TEST_ADD(SystemTest, testIteratorSphericalNr3)
TEST_ADD(SystemTest, testIteratorPlanarAm1)
TEST_ADD(SystemTest, testIteratorPlanarAm3)
TEST_ADD(SystemTest, testSweepSpherical)
TEST_END(SystemTest)

#endif
