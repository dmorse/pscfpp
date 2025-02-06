#ifndef RPG_IMPOSED_FIELDS_GENERATOR_TEST_H
#define RPG_IMPOSED_FIELDS_GENERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/scft/iterator/ImposedFieldsGenerator.h>
#include <rpg/System.h>

#include <prdc/crystal/BFieldComparison.h>

#include <pscf/iterator/FieldGenerator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class ImposedFieldsGeneratorTest : public UnitTest
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
      ImposedFieldsGenerator<1> ext(system);
   }

   void testReadParameters() // test ExtGenFilmBase::readParameters()
   {
      printMethod(TEST_FUNC);

      // Set up external field generator from file
      System<1> system;
      createSystem(system, "in/system1DGen");
      ImposedFieldsGenerator<1> gen(system);

      std::ifstream in;
      openInputFile("in/generator1", in);
      gen.readParam(in);
      in.close();

      // Check that the parameters that are publicly accessible were 
      // read correctly
      TEST_ASSERT(gen.type() == "film");
      TEST_ASSERT(gen.fieldGenerator1().type() == FieldGenerator::Mask);
      TEST_ASSERT(gen.fieldGenerator2().type() == FieldGenerator::External);
      
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 0;
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_bottom",ids), 5.0));
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_top",ids), 2.0));
      ids[0] = 1;
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_bottom",ids), 0.0));
      TEST_ASSERT(eq(gen.fieldGenerator2().getParameter("chi_top",ids), 10.0));
   }

   void testSolve1D() // solve a 1D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve1D.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGen");

      // Read initial guess
      system.readWBasis("in/wIn1D.bf");

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
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.fHelmholtz() - 3.87784944222) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.pressure() + 12.1117881919) << "\n";
      }
      TEST_ASSERT(std::abs(system.fHelmholtz() - 3.87784944222) < 1e-5);
      TEST_ASSERT(std::abs(system.pressure() + 12.1117881919) < 1e-4);

      // output mask field for reference
      system.domain().fieldIo().writeFieldBasis("out/mask_1D.bf", 
                                                system.mask().basis(),
                                                system.domain().unitCell());
      
      // output external field for reference
      system.domain().fieldIo().writeFieldsBasis("out/h_1D.bf", 
                                                 system.h().basis(),
                                                 system.domain().unitCell());
   }

   void testSolve2D() // solve a 2D system with an ImposedFieldsGenerator
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolve2D.log");
      
      // Set up system with some data
      System<2> system;
      createSystem(system, "in/system2DGen");

      // Read initial guess
      system.readWBasis("in/wIn2D.bf");

      // Solve
      system.iterate();
      system.writeWBasis("out/w2D.bf");
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));
      
      // Check that lattice parameters are correct
      double aErr = system.domain().unitCell().parameter(0) - 1.63536608507;
      TEST_ASSERT(std::abs(aErr) < 1e-5);
      TEST_ASSERT(eq(system.domain().unitCell().parameter(1), 2.0));

      // Check converged field is correct by comparing to reference
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/wRef2D.bf", 
                                                wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
     
      double epsilon = 1.0E-4; 
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
                   << (system.fHelmholtz() - 3.91021919092) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.pressure() + 12.4986763317) << "\n";
      }
      TEST_ASSERT(std::abs(system.fHelmholtz() - 3.91021919092) < 1e-5);
      TEST_ASSERT(std::abs(system.pressure() + 12.4986763317) < 1e-4);
   }

   void testSweep() // test sweep along chiBottom and lattice parameter
   {
      // NOTE: this also tests that the ParameterModifier methods work
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSweep.log");
      
      // Set up system
      System<1> system;
      createSystem(system, "in/system1DGen");

      // Read initial guess
      system.readWBasis("out/w1D.bf");

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

      double epsilon = 1.0E-5; 
      if (verbose() > 0 || diff > epsilon) {
         std::cout << "\n";
         std::cout << "diff    = " << diff << "\n";
         std::cout << "epsilon = " << epsilon << "\n";
      }
      TEST_ASSERT(diff < epsilon);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.fHelmholtz() - 3.84231924345) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.pressure() + 11.5127288016) << "\n";
      }
      TEST_ASSERT(std::abs(system.fHelmholtz() - 3.84231924345) < 1e-5);
      TEST_ASSERT(std::abs(system.pressure() + 11.5127288016) < 1e-4);
   }

   void testSolveWithFBulk() // solve a 1D system w/ flexible film thickness
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/GeneratorTestSolveWithFBulk.log");
      
      // Set up system with some data
      System<1> system;
      createSystem(system, "in/system1DGenFBulk");

      // Read initial guess
      system.readWBasis("in/wIn1D_3.bf");

      // Iterate to a solution
      system.iterate();
      
      // Check that the right film thickness was found
      double paramErr = system.domain().unitCell().parameter(0) - 2.061207269;
      if (verbose() > 0) {
         std::cout << "\nFilm thickness error = " << paramErr << "\n";
      }
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
         std::cout << "\nMax field error = " << bComparison.maxDiff() << "\n";
      }
      system.domain().fieldIo().writeFieldsBasis("out/w1D_fBulk.bf", 
                                                 system.w().basis(), 
                                                 system.domain().unitCell());
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);

      // Check thermo parameters
      if (verbose() > 0) {
         std::cout << "\nFree energy error = " 
                   << (system.fHelmholtz() - 3.80033554388) << "\n";
         std::cout << "\nPressure error = " 
                   << (system.pressure() + 12.9408830685) << "\n";
      }
      TEST_ASSERT(std::abs(system.fHelmholtz() - 3.80033554388) < 1e-5);
      TEST_ASSERT(std::abs(system.pressure() + 12.9408830685) < 1e-4);
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

TEST_BEGIN(ImposedFieldsGeneratorTest)
TEST_ADD(ImposedFieldsGeneratorTest, testConstructor)
TEST_ADD(ImposedFieldsGeneratorTest, testReadParameters)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve1D)
TEST_ADD(ImposedFieldsGeneratorTest, testSolve2D)
TEST_ADD(ImposedFieldsGeneratorTest, testSweep)
TEST_ADD(ImposedFieldsGeneratorTest, testSolveWithFBulk)
TEST_END(ExtGenFilmTest)

#endif
