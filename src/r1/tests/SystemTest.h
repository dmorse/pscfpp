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
#include <r1/sweep/LinearSweep.h>
#include <r1/sweep/SweepParameter.h>

#include <util/misc/Exception.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

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

   /*
   * REMESH_W and EXTEND_W with invalid numbers of grid points.
   */
   void testRemeshExtendInvalid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestRemeshExtendInvalid.log");

      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      std::string filename = "out/SystemTestRemeshExtendInvalid.w";
      std::istringstream in1("READ_W in/planar.w\nREMESH_W 51 "
                             + filename + "\nFINISH\n");
      sys.readCommands(in1);
      std::ifstream file;
      openInputFile(filename, file);
      std::ostringstream before;
      before << file.rdbuf();
      file.close();

      // Each command must throw, and leave the existing file unchanged
      std::string commands[3] = {"REMESH_W 1 ", "REMESH_W 0 ",
                                 "EXTEND_W -1 "};
      for (int k = 0; k < 3; ++k) {
         std::istringstream in2(commands[k] + filename + "\nFINISH\n");
         bool thrown = false;
         try {
            sys.readCommands(in2);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown);
         openInputFile(filename, file);
         std::ostringstream after;
         after << file.rdbuf();
         file.close();
         TEST_ASSERT(after.str() == before.str());
      }
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
      TEST_ASSERT(std::abs((q01 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q10 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q11 - q00)/q00) < 1.0E-6);

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
      TEST_ASSERT(std::abs((q01 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q10 - q00)/q00) < 1.0E-6);
      TEST_ASSERT(std::abs((q11 - q00)/q00) < 1.0E-6);

      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(sys.cField(0));
      double sum1 = domain.spatialAverage(sys.cField(1));
      Log::file() << "Volume fraction of block 0 = " << sum0 << "\n";
      Log::file() << "Volume fraction of block 1 = " << sum1 << "\n";
   }

   /*
   * BinaryRelaxIterator with a nonzero diagonal element of chi.
   */
   void testIteratorPlanarBrDiagonalChi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarBrDiagonalChi.log");

      // Diblock with chi(0,0) = 6, chi(0,1) = 20, chi(1,1) = 0
      System sys;
      std::istringstream in(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0 nPolymer 1"
         " nSolvent 0 Polymer{ type linear nBlock 2"
         " blocks 0 0.5 1 0.5 phi 1.0 } ds 0.01 }"
         " Interaction{ chi( 0 0 6.0 0 1 20.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 1.7 nx 101 }"
         " BinaryRelaxIterator{ epsilon 1.0e-8 maxIter 40000"
         " lambdaPlus 2.0 lambdaMinus 2.0 } }\n");
      sys.readParam(in);

      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());
      std::ifstream wIn;
      openInputFile("in/planar.w", wIn);
      fieldIo.readFields(sys.wFields(), wIn);
      wIn.close();

      TEST_ASSERT(sys.iterator().solve() == 0);

      // SCF equations: w1 - w0 = (chi01 - chi00) c0 + (chi11 - chi01) c1
      // and c0 + c1 = 1
      double c0, c1, r;
      double maxExchange = 0.0;
      double maxIncompress = 0.0;
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         c0 = sys.cField(0)[i];
         c1 = sys.cField(1)[i];
         r = sys.wField(1)[i] - sys.wField(0)[i] - 14.0*c0 + 20.0*c1;
         if (fabs(r) > maxExchange) maxExchange = fabs(r);
         r = c0 + c1 - 1.0;
         if (fabs(r) > maxIncompress) maxIncompress = fabs(r);
      }
      Log::file() << "max exchange residual          = "
                  << maxExchange << "\n"
                  << "max incompressibility residual = "
                  << maxIncompress << std::endl;
      TEST_ASSERT(maxExchange < 1.0E-6);
      TEST_ASSERT(maxIncompress < 1.0E-6);
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

      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

      openOutputFile("out/finalPlanarNr1.w", out);
      fieldIo.writeFields(sys.wFields(), out);
      out.close();

      openOutputFile("out/finalPlanarNr1.c", out);
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

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

      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

      std::ofstream out;
      out.open("out/planarNr2.c");
      fieldIo.writeFields(sys.cFields(), out);
      out.close();

   }

   /*
   * BinaryRelaxIterator parameters after a solve that reduced them.
   */
   void testIteratorPlanarBrLambda()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestIteratorPlanarBrLambda.log");

      // Diblock with chi 20, lambdaPlus = lambdaMinus = 6. An early
      // update (the third) increases the error, so that the iterator
      // reduces them.
      System sys;
      std::istringstream in(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0 nPolymer 1"
         " nSolvent 0 Polymer{ type linear nBlock 2"
         " blocks 0 0.5 1 0.5 phi 1.0 } ds 0.01 }"
         " Interaction{ chi( 0 1 20.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 1.7 nx 101 }"
         " BinaryRelaxIterator{ epsilon 1.0e-8 maxIter 40000"
         " lambdaPlus 6.0 lambdaMinus 6.0 } }\n");
      sys.readParam(in);

      FieldIo fieldIo;
      fieldIo.associate(sys.domain(), sys.fileMaster());
      std::ifstream wIn;
      openInputFile("in/planar.w", wIn);
      fieldIo.readFields(sys.wFields(), wIn);
      wIn.close();

      TEST_ASSERT(sys.iterator().solve() == 0);

      // The parameters written after the solve are those read
      std::ostringstream out;
      sys.iterator().writeParam(out);
      Log::file() << out.str();
      std::string labels[2] = {"lambdaPlus", "lambdaMinus"};
      for (int i = 0; i < 2; ++i) {
         std::string const & text = out.str();
         std::size_t pos = text.find(labels[i]);
         TEST_ASSERT(pos != std::string::npos);
         std::istringstream value(text.substr(pos + labels[i].size()));
         double lambda = 0.0;
         value >> lambda;
         TEST_ASSERT(eq(lambda, 6.0));
      }
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
      int status = sys.readCommands(in);
      in.close();
      TEST_ASSERT(status == 0);
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

      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

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
      int status = sys.readCommands(in);
      in.close();
      TEST_ASSERT(status == 0);
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

      int error = sys.iterator().solve();
      TEST_ASSERT(error == 0);

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
      int status = sys.readCommands(in);
      in.close();
      TEST_ASSERT(status == 0);
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
      int status = sys.readCommands(in);
      in.close();
      TEST_ASSERT(status == 0);
   }

   /*
   * Command scripts that end without a newline, or with an argument
   * that cannot be read.
   */
   void testReadCommandsInput()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadCommandsInput.log");

      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      // A last command without a newline must be executed
      std::istringstream in1("READ_W in/planar.w\nITERATE");
      sys.readCommands(in1);
      TEST_ASSERT(sys.fHelmholtz() != 0.0);

      // A missing argument at the end of the script must be reported
      std::istringstream in2("WRITE_W out/SystemTestReadCommandsInput.w\n"
                             "REMESH_W 51");
      bool thrown = false;
      try {
         sys.readCommands(in2);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      // So must an argument that cannot be read
      std::istringstream in3("REMESH_W x out/SystemTestReadCommandsInput.w"
                             "\nFINISH\n");
      thrown = false;
      try {
         sys.readCommands(in3);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

   /*
   * Read field files with a missing or invalid last value.
   */
   void testReadFieldsIncomplete()
   {
      printMethod(TEST_FUNC);

      Domain domain;
      domain.setPlanarParameters(0.0, 1.0, 3);
      FileMaster fileMaster;
      FieldIo fieldIo;
      fieldIo.associate(domain, fileMaster);
      DArray< DArray<double> > fields;
      fields.allocate(2);
      fields[0].allocate(3);
      fields[1].allocate(3);

      // Complete files, with and without a final newline, and files
      // with a missing or invalid last value
      std::string begin = "nx 3\nnm 2\n0 1.0 2.0\n1 3.0 4.0\n2 5.0 ";
      std::string ends[4] = {"6.0\n", "6.0", "\n", "x\n"};
      for (int k = 0; k < 4; ++k) {
         std::istringstream in(begin + ends[k]);
         bool thrown = false;
         try {
            fieldIo.readFields(fields, in);
         } catch (Exception& e) {
            thrown = true;
         }
         if (k < 2) {
            TEST_ASSERT(!thrown);
            TEST_ASSERT(eq(fields[0][2], 5.0));
            TEST_ASSERT(eq(fields[1][2], 6.0));
         } else {
            TEST_ASSERT(thrown);
         }
      }
   }

   /*
   * Free energy and pressure after the COMPUTE command.
   */
   void testComputeThermo()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestComputeThermo.log");

      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      std::istringstream in1("READ_W in/planar.w\nITERATE\nFINISH\n");
      sys.readCommands(in1);
      double f0 = sys.fHelmholtz();

      // Perturb the converged w fields, and solve without iteration
      int nx = sys.domain().nx();
      for (int i = 0; i < nx; ++i) {
         sys.wField(0)[i] += 0.5*cos(Constants::Pi*double(i)/double(nx-1));
      }
      std::istringstream in2("COMPUTE\nFINISH\n");
      sys.readCommands(in2);
      double f1 = sys.fHelmholtz();
      double p1 = sys.pressure();

      // The values must be those of the perturbed fields
      sys.computeFreeEnergy();
      TEST_ASSERT(std::abs(f1 - f0) > 1.0E-4);
      TEST_ASSERT(std::abs(f1 - sys.fHelmholtz()) < 1.0E-12);
      TEST_ASSERT(std::abs(p1 - sys.pressure()) < 1.0E-12);
   }

   /*
   * Status returned by readCommands: 0 if no command failed, 1 after
   * an unknown command or an ITERATE command that did not converge.
   */
   void testReadCommandsStatus()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestReadCommandsStatus.log");

      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      std::istringstream in1("READ_W in/planar.w\nCOMPUTE\nFINISH\n");
      TEST_ASSERT(sys.readCommands(in1) == 0);
      std::istringstream in2("COMPUTE\n");
      TEST_ASSERT(sys.readCommands(in2) == 0);
      std::istringstream in3("COMPUTEE\nFINISH\n");
      TEST_ASSERT(sys.readCommands(in3) == 1);

      // An NrIterator with maxItr 1 cannot converge
      System sys2;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0 nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks 0 0.5 1 0.5 phi 1.0 }"
         " ds 0.01 }"
         " Interaction{ chi( 0 1 30.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 0.8 nx 101 }"
         " NrIterator{ epsilon 1.0e-12 maxItr 1 } }\n");
      sys2.readParam(param);
      sys2.fileMaster().setInputPrefix(filePrefix());
      sys2.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in4("READ_W in/planar.w\nITERATE\nFINISH\n");
      TEST_ASSERT(sys2.readCommands(in4) == 1);
   }

   /*
   * Free energy and pressure of a homogeneous blend of two homopolymers
   * of length 1 with volume fractions 1 - 5.0E-9 and 5.0E-9, in a field
   * w = 0, which include the terms of the dilute species.
   */
   void testThermoSmallPhi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestThermoSmallPhi.log");

      System sys;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0 nPolymer 2"
         " Polymer{ type linear nBlock 1 blocks 0 1.0 phi 0.999999995 }"
         " Polymer{ type linear nBlock 1 blocks 1 1.0 phi 5.0e-9 }"
         " ds 0.05 } Interaction{ chi( 0 1 10.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 1.0 nx 5 }"
         " NrIterator{ epsilon 0.00001 } }\n");
      sys.readParam(param);
      for (int m = 0; m < 5; ++m) {
         sys.wFields()[0][m] = 0.0;
         sys.wFields()[1][m] = 0.0;
      }
      sys.compute();
      sys.computeFreeEnergy();
      double phi = 5.0E-9;
      double fExact = (1.0 - phi)*(std::log(1.0 - phi) - 1.0)
                    + phi*(std::log(phi) - 1.0) + 10.0*phi*(1.0 - phi);
      double pExact = 1.0 - 10.0*phi*(1.0 - phi);
      TEST_ASSERT(std::abs(sys.fHelmholtz() - fExact) < 5.0E-11);
      TEST_ASSERT(std::abs(sys.pressure() - pExact) < 5.0E-11);
   }

   /*
   * Commands that write c fields after READ_W, without a solution.
   */
   void testWriteCAfterReadW()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestWriteCAfterReadW.log");

      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr2.prm", in);
      sys.readParam(in);
      in.close();
      sys.fileMaster().setInputPrefix(filePrefix());
      sys.fileMaster().setOutputPrefix(filePrefix());

      std::istringstream in1("READ_W in/planar.w\nCOMPUTE\n"
                             "WRITE_C out/SystemTestWriteCAfterReadW.c\n"
                             "FINISH\n");
      sys.readCommands(in1);

      // After READ_W, WRITE_C and WRITE_THERMO must fail
      std::string commands[2] = {
         "WRITE_C out/SystemTestWriteCAfterReadW.c\n",
         "WRITE_THERMO out/SystemTestWriteCAfterReadW.thermo\n"};
      for (int k = 0; k < 2; ++k) {
         std::istringstream in2("READ_W in/planar.w\n" + commands[k]
                                + "FINISH\n");
         bool thrown = false;
         try {
            sys.readCommands(in2);
         } catch (Exception& e) {
            thrown = true;
         }
         TEST_ASSERT(thrown);
      }

      // After COMPUTE, they succeed
      std::istringstream in3("COMPUTE\n" + commands[0] + commands[1]
                             + "FINISH\n");
      sys.readCommands(in3);
   }

   void testSweepParameterIdRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestSweepParameterIdRange.log");

      // Diblock melt. Parameters with an index out of range must be
      // rejected by getInitial, with a message that says so.
      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr1.prm", in);
      sys.readParam(in);
      in.close();
      std::string bad[] = {"phi_polymer 1 0.1", "block 0 2 0.1"};
      for (int i = 0; i < 2; ++i) {
         SweepParameter parameter(sys);
         std::istringstream in(bad[i]);
         in >> parameter;
         std::string message;
         try {
            parameter.getInitial();
         } catch (Exception& e) {
            message = e.message();
         }
         TEST_ASSERT(message.find("out of range") != std::string::npos);
      }
   }

   /*
   * SweepParameter::serialize took the archive by value, so that a
   * SweepParameter could not be written to a MemoryOArchive (whose copy
   * constructor is private).
   */
   void testSweepParameterSerialize()
   {
      printMethod(TEST_FUNC);

      SweepParameter p;
      MemoryOArchive oArchive;
      oArchive.allocate(64);
      oArchive << p;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin()
                  + sizeof(unsigned int) + sizeof(int) + 2*sizeof(double));
   }

   /*
   * Loading a SweepParameter into a new object wrote the identifiers to
   * an array that was not allocated.
   */
   void testSweepParameterLoad()
   {
      printMethod(TEST_FUNC);

      // Record of a block length parameter (type 0) with ids 0 and 1
      unsigned int type = 0;
      int nId = 2;
      int id0 = 0;
      int id1 = 1;
      double initial = 0.5;
      double change = 0.25;
      MemoryOArchive oArchive;
      oArchive.allocate(64);
      oArchive << type;
      oArchive << nId;
      oArchive << id0;
      oArchive << id1;
      oArchive << initial;
      oArchive << change;

      SweepParameter p;
      MemoryIArchive iArchive;
      iArchive = oArchive;
      iArchive >> p;
      iArchive.release();
      TEST_ASSERT(p.type() == "block");
      TEST_ASSERT(p.id(0) == 0);
      TEST_ASSERT(p.id(1) == 1);
      TEST_ASSERT(p.initial() == 0.5);
      TEST_ASSERT(p.change() == 0.25);
   }

   void testLinearSweepPhiMixed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestLinearSweepPhiMixed.log");

      // Diblock (closed) in an open solvent: a sweep of the volume
      // fraction of the diblock alone must be accepted
      System sys;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers 1.0 1.0"
         " nPolymer 1 nSolvent 1"
         " Polymer{ type linear nBlock 2 blocks 0 2.0 1 3.0 phi 0.8 }"
         " Solvent{ monomerId 1 size 0.05 mu 1.0 }"
         " vMonomer 0.05 ds 0.01 }"
         " Interaction{ chi( 0 1 4.0 ) }"
         " Domain{ mode Planar xMin 0.0 xMax 4.0 nx 33 }"
         " NrIterator{ epsilon 0.00001 } }\n");
      sys.readParam(param);
      LinearSweep sweep(sys);
      std::istringstream in("ns 2 nParameter 1"
                            " parameters phi_polymer 0 -0.1\n");
      bool thrown = false;
      try {
         sweep.readParameters(in);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
   }

   void testLinearSweepPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/SystemTestLinearSweepPhiSum.log");

      // Closed diblock melt: a change of phi of 5.0E-7 must be rejected
      System sys;
      std::ifstream in;
      openInputFile("in/planar_nr1.prm", in);
      sys.readParam(in);
      in.close();
      LinearSweep sweep(sys);
      std::istringstream sweepIn("ns 4 nParameter 1"
                                 " parameters phi_polymer 0 5.0e-7\n");
      std::string message;
      try {
         sweep.readParameters(sweepIn);
      } catch (Exception& e) {
         message = e.message();
      }
      TEST_ASSERT(message.find("add up to zero") != std::string::npos);
   }

};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor)
TEST_ADD(SystemTest, testRemeshExtendInvalid)
TEST_ADD(SystemTest, testReadParameters)
TEST_ADD(SystemTest, testSolveMdePlanar)
TEST_ADD(SystemTest, testSolveMdeSpherical)
TEST_ADD(SystemTest, testIteratorPlanarBrDiagonalChi)
TEST_ADD(SystemTest, testIteratorPlanarNr1)
TEST_ADD(SystemTest, testIteratorPlanarNrNan)
TEST_ADD(SystemTest, testIteratorPlanarNr2)
TEST_ADD(SystemTest, testIteratorPlanarBrLambda)
TEST_ADD(SystemTest, testIteratorPlanarNr3)
TEST_ADD(SystemTest, testIteratorSphericalNr1)
TEST_ADD(SystemTest, testIteratorSphericalNr3)
TEST_ADD(SystemTest, testIteratorPlanarAm1)
TEST_ADD(SystemTest, testIteratorPlanarAm3)
TEST_ADD(SystemTest, testSweepSpherical)
TEST_ADD(SystemTest, testThermoSmallPhi)
TEST_ADD(SystemTest, testWriteCAfterReadW)
TEST_ADD(SystemTest, testReadCommandsInput)
TEST_ADD(SystemTest, testReadFieldsIncomplete)
TEST_ADD(SystemTest, testComputeThermo)
TEST_ADD(SystemTest, testReadCommandsStatus)
TEST_ADD(SystemTest, testSweepParameterIdRange)
TEST_ADD(SystemTest, testSweepParameterSerialize)
TEST_ADD(SystemTest, testSweepParameterLoad)
TEST_ADD(SystemTest, testLinearSweepPhiMixed)
TEST_ADD(SystemTest, testLinearSweepPhiSum)
TEST_END(SystemTest)

#endif
