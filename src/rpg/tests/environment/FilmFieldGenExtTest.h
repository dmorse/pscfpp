#ifndef RPG_FILM_FIELD_GEN_EXT_TEST_H
#define RPG_FILM_FIELD_GEN_EXT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/environment/FilmFieldGenExt.h>
#include <rpg/system/System.h>

#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/BFieldComparison.h>
#include <prdc/environment/FieldGenerator.h>

#include <util/misc/Exception.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class FilmFieldGenExtTest : public UnitTest
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
      FilmFieldGenExt<1> ext(system);
   }

   void testReadParameters() // test FilmFieldGenExtBase::readParameters()
   {
      printMethod(TEST_FUNC);

      // Set up external field generator from file
      System<1> system;
      createSystem(system, "in/system1D");
      FilmFieldGenExt<1> ext(system);

      std::ifstream in;
      openInputFile("in/filmExt1Asym", in);
      ext.readParameters(in);
      in.close();

      // Check that everything was read in correctly
      TEST_ASSERT(ext.normalVecId() == 0);
      TEST_ASSERT(eq(ext.interfaceThickness(), 0.2));
      TEST_ASSERT(eq(ext.excludedThickness(), 0.4));
      TEST_ASSERT(ext.chiBottom().capacity() == 2);
      TEST_ASSERT(ext.chiTop().capacity() == 2);
      TEST_ASSERT(eq(ext.chiBottom(0), 5.0));
      TEST_ASSERT(eq(ext.chiBottom(1), 0.0));
      TEST_ASSERT(eq(ext.chiTop(0), 2.0));
      TEST_ASSERT(eq(ext.chiTop(1), 10.0));
      TEST_ASSERT(ext.type() == FieldGenerator::External);
   }

   void testCheckCompatibility() // test FieldGenerator::checkCompatibility()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/extTestCheckCompatibility.log");

      // Set up 1D external field with symmetric walls
      System<1> system1;
      createSystem(system1, "in/system1D");
      
      FilmFieldGenExt<1> ext1(system1);
      createFilmFieldGenExt(ext1, "in/filmExt1Sym");
      
      Log::file() << "Testing system 1:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(ext1,false));
      TEST_ASSERT(ext1.hasSymmetricWalls());
      TEST_ASSERT(!ext1.isAthermal());
      TEST_ASSERT(ext1.needsUpdate());

      // Set up 2D external field with asymmetric walls and an 
      // incompatible space group
      System<2> system2;
      createSystem(system2, "in/system2D_1");

      FilmFieldGenExt<2> ext2(system2);
      createFilmFieldGenExt(ext2, "in/filmExt2Asym");

      Log::file() << "Testing system 2:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(ext2,true));


      // Set up 3D external field with asymmetric walls and a compatible
      // space group
      System<3> system3;
      createSystem(system3, "in/system3D_3");

      FilmFieldGenExt<3> ext3(system3);
      createFilmFieldGenExt(ext3, "in/filmExt3Asym");

      Log::file() << "Testing system 3:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(ext3,false));
      TEST_ASSERT(!ext3.hasSymmetricWalls());
      TEST_ASSERT(!ext3.isAthermal());
      TEST_ASSERT(ext3.needsUpdate());
      TEST_ASSERT(ext3.normalVecId() == 2);
   }

   void testGenerate() // test FieldGenerator::generate()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/extTestGenerate.log");

      // Set up 3D external field with a compatible system
      System<3> system1;
      createSystem(system1, "in/system3D_3");
      system1.mask().allocateBasis(1920);
      system1.mask().allocateRGrid(system1.domain().mesh().dimensions());
      system1.mask().setFieldIo(system1.domain().fieldIo());
      system1.mask().readBasis("in/maskRef3.bf");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(4.0);
      //system1.setUnitCell(UnitCell<3>::Hexagonal, parameters);
      system1.setUnitCell(parameters);

      // Set up external field
      FilmFieldGenExt<3> ext1(system1);
      createFilmFieldGenExt(ext1, "in/filmExt3Asym");
      ext1.generate();
      TEST_ASSERT(!ext1.needsUpdate());
      TEST_ASSERT(system1.h().isAllocatedRGrid());
      TEST_ASSERT(system1.h().isAllocatedBasis());
      TEST_ASSERT(system1.h().hasData());
      TEST_ASSERT(system1.h().isSymmetric());

      // Check that the generated external fields are correct
      UnitCell<3> unitCell; // UnitCell object to pass to FieldIo functions
      DArray<DArray<double> > hFieldsCheck; // Copy of reference field
      system1.domain().fieldIo().readFieldsBasis("in/hRef1.bf", 
                                       hFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system1.h().basis(), hFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);


      // Set up 3D system with athermal walls
      System<3> system2;
      createSystem(system2, "in/system3D_3");
      system2.mask().allocateBasis(1920);
      system2.mask().allocateRGrid(system2.domain().mesh().dimensions());
      //system2.setUnitCell(UnitCell<3>::Hexagonal, parameters);
      system2.setUnitCell(parameters);

      // Set up external field
      FilmFieldGenExt<3> ext2(system2);
      createFilmFieldGenExt(ext2, "in/filmExt3Athermal");
      ext2.generate();

      TEST_ASSERT(!ext2.needsUpdate());
      TEST_ASSERT(!system2.h().isAllocatedRGrid());
      TEST_ASSERT(!system2.h().isAllocatedBasis());
      TEST_ASSERT(!system2.h().hasData());
   }

   void testRegenerate() // test generate() with parameter changes
   {
      // Note: this test also tests the methods getParameter and 
      // setParameter, which are used by Sweep to modify chiBottom
      // and chiTop.

      printMethod(TEST_FUNC);
      openLogFile("out/extTestRegenerate.log");

      // Set up 1D external field with a compatible system
      System<1> system;
      createSystem(system, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.9);
      //system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system.setUnitCell(parameters);
      system.mask().allocateBasis(37);
      system.mask().allocateRGrid(system.domain().mesh().dimensions());
      system.mask().setFieldIo(system.domain().fieldIo());
      system.mask().readBasis("in/maskRef1.bf");

      // Set up external field generator
      FilmFieldGenExt<1> ext(system);
      createFilmFieldGenExt(ext, "in/filmExt1Sym");
      ext.generate();

      // Change lattice parameter and update
      parameters[0] = 3.0;
      //system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system.setUnitCell(parameters);
      TEST_ASSERT(ext.needsUpdate());
      ext.generate();
      TEST_ASSERT(!ext.needsUpdate());

      // Check that updated external fields are correct
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray<DArray<double> > hFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/hRef2.bf", hFieldsCheck, 
		                                unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.h().basis(), hFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);

      // Change some chi values and update
      DArray<int> ids;
      ids.allocate(1);
      ids[0] = 1;
      TEST_ASSERT(eq(ext.getParameter("chi_top",ids), 0.0));
      TEST_ASSERT(eq(ext.getParameter("chi_bottom",ids), 0.0));
      
      ext.setParameter("chi_top", ids, -2.0);
      ext.setParameter("chi_bottom", ids, -2.0);
      TEST_ASSERT(eq(ext.getParameter("chi_top",ids), -2.0));
      TEST_ASSERT(eq(ext.getParameter("chi_bottom",ids), -2.0));

      TEST_ASSERT(ext.needsUpdate());
      ext.generate();
      TEST_ASSERT(!ext.needsUpdate());

      // Check that updated external fields are is correct
      UnitCell<1> unitCell2; // UnitCell object to pass to FieldIo functions
      DArray<DArray<double> > hFieldsCheck2; // Copy of reference field
      system.domain().fieldIo().readFieldsBasis("in/hRef3.bf", 
		                                hFieldsCheck2, 
		                                unitCell2);
      BFieldComparison bComparison2(0); // object to compare fields
      bComparison2.compare(system.h().basis(), hFieldsCheck2);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison2.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison2.maxDiff() < 1.0E-7);
   }

   void testStress() // test FilmFieldGenExt::stress
   {
      printMethod(TEST_FUNC);
      openLogFile("out/extTestStress.log");

      // Set up 2D field generator with a compatible system
      System<2> system;
      createSystem(system, "in/system2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(1.6353661975);
      parameters.append(2.0);
      //system.setUnitCell(UnitCell<2>::Rectangular, parameters);
      system.setUnitCell(parameters);

      // Set up mask
      system.mask().allocateBasis(561);
      system.mask().allocateRGrid(system.domain().mesh().dimensions());
      system.mask().setFieldIo(system.domain().fieldIo());
      system.mask().readBasis("in/maskRef2.bf");

      // Set up external field generator
      FilmFieldGenExt<2> ext(system);
      createFilmFieldGenExt(ext, "in/filmExt2Sym");
      ext.generate();

      // Read w field and solve MDEs, so system can calculate fHelmholtz
      system.w().readBasis("in/wIn2D.bf");
      system.compute();

      // Call stress and check that the result is correct
      TEST_ASSERT(eq(ext.stress(0),0.0));
      TEST_ASSERT(std::abs(ext.stress(1) + 0.163552584) < 1e-5);
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

   // Read parameter file section to create a FilmFieldGenExt object
   template <int D>
   void createFilmFieldGenExt(FilmFieldGenExt<D>& ext, std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      ext.readParameters(in);
      in.close();
   }

   // Determine if we get expected result when running checkCompatibility.
   // Function accepts a boolean indicating whether we expect it to throw 
   // an error or not, and returns a boolean indicating whether the 
   // function demonstrated the expected behavior.
   template <int D>
   bool checkCheckCompatibility(FilmFieldGenExt<D>& ext, bool expectError)
   {
      bool pass = true;
      if (expectError) {
         try {
            ext.checkCompatibility();
            // This is expected to fail. If it succeeds, the test fails.
            pass = false;
         } catch (Exception& e) {
            Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                        << std::endl;
         }
      } else {
         try {
            ext.checkCompatibility();
            // This should succeed. If not, the test fails. 
         } catch (Exception& e) {
            pass = false;
         }
      }
      return pass;
   }

};

TEST_BEGIN(FilmFieldGenExtTest)
TEST_ADD(FilmFieldGenExtTest, testConstructor)
TEST_ADD(FilmFieldGenExtTest, testReadParameters)
TEST_ADD(FilmFieldGenExtTest, testCheckCompatibility)
TEST_ADD(FilmFieldGenExtTest, testGenerate)
TEST_ADD(FilmFieldGenExtTest, testRegenerate)
TEST_ADD(FilmFieldGenExtTest, testStress)
TEST_END(FilmFieldGenExtTest)

#endif
