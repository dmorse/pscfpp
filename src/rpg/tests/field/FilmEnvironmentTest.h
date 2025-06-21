#ifndef RPG_FILM_ENVIRONMENT_TEST_H
#define RPG_FILM_ENVIRONMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/field/MixAndMatchEnvs.h>
#include <rpg/System.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldComparison.h>

#include <pscf/environment/FieldGenerator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

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

      // Set up film environment from file
      System<1> system;
      createSystem(system, "in/film/system1DEnv");
      FilmEnvironment<1> env(system);

      std::ifstream in;
      openInputFile("in/film/environment1", in);
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
      createSystem(system, "in/film/system1DEnv");

      // Read initial guess
      system.readWBasis("in/film/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      Log::file() << system.mask().phiTot() << std::endl;
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/film/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free energy error = " 
                   << (system.fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.1117881919) < 1e-4);
   }

   void testSolve2D() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2D.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/film/system2DEnv");

      // Read initial guess
      system.readWBasis("in/film/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free energy error = " 
                   << (system.fHelmholtz() - 3.91037539514) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.8397354494) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.91037539514) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.8397354494) < 1e-4);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweep.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/film/system1DEnv");

      // Read initial guess
      system.readWBasis("in/film/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-4; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.87318676998) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.0498211637) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulk.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/film/system1DEnvFBulk");

      // Read initial guess
      system.readWBasis("in/film/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness error = " << paramErr << "\n";
      }
      TEST_ASSERT(abs(paramErr) < 1e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.8059299672) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/film/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.80033554388) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.9408830685) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.80033554388) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.9408830685) < 1e-4);
   }

   void testSolve1DGrid() // solve a 1D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve1DGrid.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/film/system1DEnvGrid");

      // Read initial guess
      system.readWBasis("in/film/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/film/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.1117881919) < 1e-4);
   }

   void testSolve2DGrid() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2DGrid.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/film/system2DEnvGrid");

      // Read initial guess
      system.readWBasis("in/film/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.91037539514) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.8397354494) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.91037539514) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.8397354494) < 1e-4);
   }

   void testSweepGrid() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweepGrid.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/film/system1DEnvGrid");

      // Read initial guess
      system.readWBasis("in/film/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-4; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.87318676998) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.0498211637) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulkGrid() // solve a 1D system w flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulkGrid.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/film/system1DEnvFBulkGrid");

      // Read initial guess
      system.readWBasis("in/film/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness error = " << paramErr << "\n";
      }
      TEST_ASSERT(abs(paramErr) < 1e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.8059299672) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/film/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1> > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/film/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "Max field error = " << diff << "\n";
      }
      TEST_ASSERT(diff < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "Free Energy error = " 
                   << (system.fHelmholtz() - 3.80033554388) << "\n";
         std::cout << "Pressure error = " 
                   << (system.pressure() + 12.9408830685) << "\n";
      }
      TEST_ASSERT(abs(system.fHelmholtz() - 3.80033554388) < 1e-5);
      TEST_ASSERT(abs(system.pressure() + 12.9408830685) < 1e-4);
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
TEST_ADD(FilmEnvironmentTest, testSolve1DGrid)
TEST_ADD(FilmEnvironmentTest, testSolve2DGrid)
TEST_ADD(FilmEnvironmentTest, testSweepGrid)
TEST_ADD(FilmEnvironmentTest, testSolveWithFBulkGrid)
TEST_END(FilmEnvironmentTest)

#endif
