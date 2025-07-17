#ifndef RPC_FILM_ENVIRONMENT_TEST_H
#define RPC_FILM_ENVIRONMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>
#include <rpc/scft/ScftThermo.h>
#include <rpc/environment/FilmEnvironment.h>

#include <prdc/crystal/BFieldComparison.h>

#include <pscf/environment/FieldGenerator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

class FilmEnvironmentTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  
      setVerbose(0); 
   }

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      System<1> system;
      FilmEnvironment<1> ext(system);
   }

   void testReadParameters() // test FilmEnvironment::readParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/FilmEnvTestReadParameters.log");

      // Set up film environment from file
      System<1> system;
      createSystem(system, "in/system1DEnv");
      FilmEnvironment<1> env(system);

      std::ifstream in;
      openInputFile("in/environment1", in);
      env.readParam(in);
      in.close();

      // Check that the parameters that are publicly accessible were 
      // read correctly
      TEST_ASSERT(env.fieldGenerator1().type() == FieldGenerator::Mask);
      TEST_ASSERT(env.fieldGenerator2().type() == FieldGenerator::External);
      
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 0;
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_bottom",ids), 5.0));
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_top",ids), 2.0));
      ids[0] = 1;
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_bottom",ids), 0.0));
      TEST_ASSERT(eq(env.fieldGenerator2().getParameter("chi_top",ids), 10.0));
   }

   void testSolve1D() // solve a 1D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve1D.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef1D.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      system.domain().fieldIo().writeFieldsBasis("out/w1D.bf", 
                                                 system.w().basis(), 
                                                 system.domain().unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-2);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.1117881919) < 1e-4);

      // output mask field for reference
      system.domain().fieldIo().writeFieldBasis("out/mask_1D.bf", 
                                                system.mask().basis(),
                                                system.domain().unitCell());
      
      // output external field for reference
      system.domain().fieldIo().writeFieldsBasis("out/h_1D.bf", 
                                                 system.h().basis(),
                                                 system.domain().unitCell());
   }

   void testSolve2D() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2D.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/system2DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      system.w().writeBasis("out/w2D.bf");
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536665618;
      TEST_ASSERT(std::abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef2D.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
     
      double epsilon = 1.0E-2; 
      double diff = bComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.91037539493) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.8397347928) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.91037539493) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.8397347928) < 1e-4);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweep.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("out/w1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRefSweep.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      double diff = bComparison.maxDiff();

      double epsilon = 1.0E-2; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.87318676998) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.0498211637) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulk.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DEnvFBulk");

      // Read initial guess
      system.w().readBasis("in/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      TEST_ASSERT(std::abs(paramErr) < 1e-5);
      TEST_ASSERT(std::abs(system.mask().phiTot() - 0.8059299672) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef1D_2.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      system.domain().fieldIo().writeFieldsBasis("out/w1D_fBulk.bf", 
                                                 system.w().basis(), 
                                                 system.domain().unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.80033554388) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.9408830685) << "\n";
      }
      TEST_ASSERT(std::abs(system.scft().fHelmholtz() - 3.80033554388) < 1e-5);
      TEST_ASSERT(std::abs(system.scft().pressure() + 12.9408830685) < 1e-4);
   }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

};

TEST_BEGIN(FilmEnvironmentTest)
TEST_ADD(FilmEnvironmentTest, testConstructor)
TEST_ADD(FilmEnvironmentTest, testReadParameters)
TEST_ADD(FilmEnvironmentTest, testSolve1D)
TEST_ADD(FilmEnvironmentTest, testSolve2D)
TEST_ADD(FilmEnvironmentTest, testSweep)
TEST_ADD(FilmEnvironmentTest, testSolveWithFBulk)
TEST_END(FilmEnvironmentTest)

#endif
