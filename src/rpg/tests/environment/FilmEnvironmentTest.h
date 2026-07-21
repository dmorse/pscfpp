#ifndef RPG_FILM_ENVIRONMENT_TEST_H
#define RPG_FILM_ENVIRONMENT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/environment/FilmEnvironment.h>
#include <rpg/system/System.h>
#include <rpg/scft/ScftThermo.h>
#include <rpg/field/Domain.h>
#include <rpg/field/FieldIo.h>
#include <rpg/field/CFields.h>
#include <rpg/field/WFields.h>
#include <rpg/field/Mask.h>

#include <prdc/field/cuda/RField.h>
#include <prdc/field/cuda/RFieldComparison.h>
#include <prdc/environment/FieldGenerator.h>

#include <util/misc/FileMaster.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;

class FilmEnvironmentTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {  setVerbose(0); }

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
      Rp::System<1, CudaTp<1> > system;
      Rp::FilmEnvironment<1, CudaTp<1> > ext(system);
   }

   void testReadParameters() // test FilmEnvironment::readParameters()
   {
      printMethod(TEST_FUNC);

      // Set up film environment from file
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnv");
      Rp::FilmEnvironment<1, CudaTp<1> > env(system);

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
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      Log::file() << system.mask().phiTot() << std::endl;
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff();
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.87784944222);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.1117881919);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.1117881919) < 1e-4);
   }

   void testSolve2D() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2D.log");
      
      // Set up system with some data
      Rp::System<2, CudaTp<2> > system;
      createSystem(system, "in/system2DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2, CudaTp<2> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.scft().fHelmholtz() - 3.91037539514);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.8397354494);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.91037539514) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.8397354494) < 1e-4);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweep.log");
      
      // Set up system
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnv");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-4; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.87318676998);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.0498211637);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulk.log");
      
      // Set up system with some data
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnvFBulk");

      // Read initial guess
      system.w().readBasis("in/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double param = system.domain().unitCell().parameter(0);
      //double paramRef = 2.061207269;  // old value
      double paramRef = 2.06121822708647739475e+00; // changed v1.3.4
      double paramErr = std::abs(param - paramRef);
      if (verbose() > 0) {
         std::cout << "\nFilm thickness       = " << Dbl(param, 20, 12);
         std::cout << "\nFilm thickness error = " << Dbl(paramErr, 20, 12) ;
      }
      TEST_ASSERT(paramErr < 1.0e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.8059299672) < 1.0e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.80033554388);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.9408830685);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.80033554388) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.9408830685) < 1e-4);
   }

   void testSolve1DGrid() // solve a 1D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve1DGrid.log");
      
      // Set up system with some data
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnvGrid");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Iterate to a solution
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 8.0951532073e-01));

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << rComparison.maxDiff();
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.87784944222);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.1117881919);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.1117881919) < 1e-4);
   }

   void testSolve2DGrid() // solve a 2D system with an FilmEnvironment
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolve2DGrid.log");
      
      // Set up system with some data
      Rp::System<2, CudaTp<2> > system;
      createSystem(system, "in/system2DEnvGrid");

      // Read initial guess
      system.w().readBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<2, CudaTp<2> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef2D.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<2> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
     
      double epsilon = 1.0E-4; 
      double diff = rComparison.maxDiff();

      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.91037539514);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.8397354494);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.91037539514) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.8397354494) < 1e-4);
   }

   void testSweepGrid() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSweepGrid.log");
      
      // Set up system
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnvGrid");

      // Read initial guess
      system.w().readBasis("in/wIn1D.bf");

      // Run the sweep function
      system.sweep();

      #if 0
      // Check converged field is correct by comparing to reference
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRefSweep.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();

      double epsilon = 1.0E-4; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < epsilon);
      #endif

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.87318676998);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.0498211637);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.87318676998) < 1e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.0498211637) < 1e-4);
   }

   void testSolveWithFBulkGrid() // solve a 1D system w flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/FilmEnvTestSolveWithFBulkGrid.log");
      
      // Set up system with some data
      Rp::System<1, CudaTp<1> > system;
      createSystem(system, "in/system1DEnvFBulkGrid");

      // Read initial guess
      system.w().readBasis("in/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double param = system.domain().unitCell().parameter(0);
      //double paramRef = 2.061207269; // old reference value
      double paramRef = 2.06121822708647739475e+00; // changed v1.3.4
      //double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      double paramErr = param - paramRef;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness       = " << Dbl(param, 20, 12);
         std::cout << "\nFilm thickness error = " << Dbl(paramErr, 20, 12);
      }
      TEST_ASSERT(abs(paramErr) < 1e-5);
      TEST_ASSERT(abs(system.mask().phiTot() - 0.8059299672) < 1e-5);

      // Check converged field is correct by comparing to ref files in in/
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< RField<1, CudaTp<1> > > cFieldsCheck; // reference fields
      system.domain().fieldIo().readFieldsRGrid("in/cRef1DFBulk.rf", 
                                                cFieldsCheck, unitCell);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(system.c().rgrid(), cFieldsCheck);
      double diff = rComparison.maxDiff();
      if (verbose() > 0) {
         std::cout << "\nMax field error = " << diff;
      }
      TEST_ASSERT(diff < 1.0E-4);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree Energy error = " 
                   << (system.scft().fHelmholtz() - 3.80033554388);
         std::cout << "\nPressure error = " 
                   << (system.scft().pressure() + 12.9408830685);
      }
      TEST_ASSERT(abs(system.scft().fHelmholtz() - 3.80033554388) < 1.0e-5);
      TEST_ASSERT(abs(system.scft().pressure() + 12.9408830685) < 1.0e-4);
   }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(Rp::System<D, CudaTp<D> >& system, std::string fname)
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
