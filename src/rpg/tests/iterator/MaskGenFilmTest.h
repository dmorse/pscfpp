#ifndef RPG_MASK_GEN_FILM_TEST_H
#define RPG_MASK_GEN_FILM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/scft/iterator/MaskGenFilm.h>
#include <rpg/scft/iterator/Iterator.h>
#include <rpg/System.h>

#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/BFieldComparison.h>

#include <pscf/iterator/FieldGenerator.h>

#include <util/misc/Exception.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class MaskGenFilmTest : public UnitTest
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
      MaskGenFilm<1> mask(system);
   }

   void testReadParameters() // test MaskGenFilmBase::readParameters()
   {
      printMethod(TEST_FUNC);

      // Set up mask from file
      System<2> system;
      MaskGenFilm<2> mask(system);

      std::ifstream in;
      openInputFile("in/filmMask1", in);
      mask.readParameters(in);
      in.close();

      // Check that everything was read in correctly
      TEST_ASSERT(eq(mask.normalVecId(),0));
      TEST_ASSERT(eq(mask.interfaceThickness(),0.2));
      TEST_ASSERT(eq(mask.excludedThickness(),0.4));
      TEST_ASSERT(eq(mask.fBulk(), 2.98468253989));
      TEST_ASSERT(mask.hasFBulk());
      TEST_ASSERT(mask.type() == FieldGenerator::Mask);
   }

   void testCheckCompatibility() // test FieldGenerator::checkCompatibility()
   {
      // This test indirectly tests various methods that are used by 
      // checkCompatibility. In the Iterator class, isFlexible, 
      // nFlexibleParams, and flexibleParams are tested. In MaskGenFilm and
      // its parent classes, the method checkCompatibility calls the methods
      // checkSpaceGroup, setFlexibleParams, and checkLatticeVectors. These
      // methods are tested by testing the compatibility of systems that are
      // incompatible for various reasons, as well as several systems that
      // are compatible.

      printMethod(TEST_FUNC);
      openLogFile("out/maskTestCheckCompatibilityLog");

      // Set up 1D mask with a compatible system and check compatibility
      System<1> system1;
      createSystem(system1, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.9);
      system1.setUnitCell(UnitCell<1>::Lamellar, parameters);
      
      MaskGenFilm<1> mask1(system1);
      createMaskGenFilm(mask1, "in/filmMask1");
      
      Log::file() << "Testing system 1:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask1,false));
      TEST_ASSERT(!mask1.isGenerated());
      TEST_ASSERT(system1.iterator().isFlexible() == true);
      TEST_ASSERT(system1.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system1.iterator().flexibleParams()[0] == true);


      // Set up 2D mask with a compatible system and check compatibility
      System<2> system2;
      createSystem(system2, "in/system2D_1");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(3.0);
      system2.setUnitCell(UnitCell<2>::Rectangular, parameters);

      MaskGenFilm<2> mask2(system2);
      createMaskGenFilm(mask2, "in/filmMask2");

      Log::file() << "Testing system 2:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask2,false));
      TEST_ASSERT(system2.iterator().isFlexible() == true);
      TEST_ASSERT(system2.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system2.iterator().flexibleParams()[0] == true);
      TEST_ASSERT(system2.iterator().flexibleParams()[1] == false);

      // Set up 2D mask with an incompatible system and check compatibility
      System<2> system3;
      createSystem(system3, "in/system2D_2");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(2.0);
      system3.setUnitCell(UnitCell<2>::Square, parameters);

      MaskGenFilm<2> mask3(system3);
      createMaskGenFilm(mask3, "in/filmMask2");

      Log::file() << "Testing system 3:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask3,true));


      // Set up 3D mask with a compatible system and check compatibility
      System<3> system4;
      createSystem(system4, "in/system3D_1");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(4.2);
      system4.setUnitCell(UnitCell<3>::Tetragonal, parameters);

      MaskGenFilm<3> mask4(system4);
      createMaskGenFilm(mask4, "in/filmMask3");

      Log::file() << "Testing system 4:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask4,false));
      TEST_ASSERT(system4.iterator().isFlexible() == true);
      TEST_ASSERT(system4.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system4.iterator().flexibleParams()[0] == true);
      TEST_ASSERT(system4.iterator().flexibleParams()[1] == false);
   

      // Set up another 3D mask with a compatible system
      System<3> system5;
      createSystem(system5, "in/system3D_2");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(3.0);
      parameters.append(4.0);
      parameters.append(1.57079632679489662); // pi/2
      system5.setUnitCell(UnitCell<3>::Monoclinic, parameters);

      MaskGenFilm<3> mask5(system5);
      createMaskGenFilm(mask5, "in/filmMask3");

      // Test compatibility
      Log::file() << "Testing system 5:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask5,false));
      TEST_ASSERT(system5.iterator().isFlexible() == true);
      TEST_ASSERT(system5.iterator().nFlexibleParams() == 2);
      TEST_ASSERT(system5.iterator().flexibleParams()[0] == true);
      TEST_ASSERT(system5.iterator().flexibleParams()[1] == true);
      TEST_ASSERT(system5.iterator().flexibleParams()[2] == false);
      TEST_ASSERT(system5.iterator().flexibleParams()[3] == false);

      // Change beta, making the lattice vectors incompatible
      parameters[3] = 1.0;
      system5.setUnitCell(UnitCell<3>::Monoclinic, parameters);
      Log::file() << "Testing system 6:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask5,true));

      // Try a different mask with normalVecId == 1, which is incompatible
      // due to the space group
      MaskGenFilm<3> mask6(system5);
      createMaskGenFilm(mask6, "in/filmMask2");
      Log::file() << "Testing system 7:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask6,true));

   }

   void testSetup() // test FieldGenerator::setup()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestSetupLog");

      // Set up 2D mask with a compatible system
      System<2> system;
      createSystem(system, "in/system2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(1.6353661975);
      parameters.append(2.0);
      system.setUnitCell(UnitCell<2>::Rectangular, parameters);

      // Set up mask
      MaskGenFilm<2> mask(system);
      createMaskGenFilm(mask, "in/filmMask2");
      mask.setup();
      
      TEST_ASSERT(mask.isGenerated());
      TEST_ASSERT(system.mask().isAllocatedBasis());
      TEST_ASSERT(system.mask().isAllocatedRGrid());
      TEST_ASSERT(system.mask().hasData());
      TEST_ASSERT(system.mask().isSymmetric());
      TEST_ASSERT(eq(system.mask().phiTot(), 7.99990525324e-01));

      // Check that the generated mask is correct
      UnitCell<2> unitCell; // UnitCell object to pass to FieldIo functions
      DArray<double> maskCheck; // Copy of reference field
      system.domain().fieldIo().readFieldBasis("in/maskRef2.bf", 
                                               maskCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.mask().basis(), maskCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);
   }

   void testUpdate() // test FieldGenerator::update()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestUpdateLog");

      // Set up 1D mask with a compatible system
      System<1> system;
      createSystem(system, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.9);
      system.setUnitCell(UnitCell<1>::Lamellar, parameters);

      // Set up mask
      MaskGenFilm<1> mask(system);
      createMaskGenFilm(mask, "in/filmMask1");
      mask.setup();
      TEST_ASSERT(!mask.updateNeeded());

      // Change lattice parameter and update
      parameters[0] = 3.0;
      system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      TEST_ASSERT(mask.updateNeeded());
      mask.update();

      // Check that updated mask is correct
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray<double> wFieldsCheck; // Copy of reference field
      system.domain().fieldIo().readFieldBasis("in/maskRef1.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.mask().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\nMax error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);
   }

   void testStressTerm() // test MaskGenFilm::stressTerm
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestStressTermLog");

      // Set up 2D mask with a compatible system
      System<2> system;
      createSystem(system, "in/system2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(1.6353661975);
      parameters.append(2.0);
      system.setUnitCell(UnitCell<2>::Rectangular, parameters);

      // Set up mask
      MaskGenFilm<2> mask(system);
      createMaskGenFilm(mask, "in/filmMask2");
      mask.setup();

      // Read w field and solve MDEs, so system can calculate fHelmholtz
      system.readWBasis("in/wIn2D.bf");
      system.compute();

      // Call stressTerm and check that the result is correct
      double stress0 = mask.stressTerm(0);
      double stress1 = mask.stressTerm(1);
      if (verbose() > 0) {
         std::cout << "\nParameter 0:"
                   << "\nstress calculated: " << stress0
                   << "\nstress expected:   " << 0.0
                   << "\nParameter 1:"
                   << "\nstress calculated: " << stress1
                   << "\nstress expected:   " << 0.4211629642 << std::endl;
      }
      TEST_ASSERT(eq(stress0,0.0));
      TEST_ASSERT(std::abs(stress1 - 0.4211629642) < 1e-5);
   }

   void testModifyStress() // test MaskGenFilm::modifyStress
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestModifyStressLog");

      // Set up 1D mask with a compatible system
      System<1> system;
      createSystem(system, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      system.setUnitCell(UnitCell<1>::Lamellar, parameters);

      // Set up mask
      MaskGenFilm<1> mask(system);
      createMaskGenFilm(mask, "in/filmMask1");
      mask.setup();

      // Read w field and solve MDEs, so system can calculate fHelmholtz
      system.readWBasis("in/wIn1D_2.bf");
      system.compute();

      // Call modifyStress with an arbitrary input stress value 
      // and test that the result is correct
      TEST_ASSERT(eq(mask.modifyStress(0, 5.0),8.94005073083));
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

   // Read parameter file section to create a MaskGenFilm object
   template <int D>
   void createMaskGenFilm(MaskGenFilm<D>& mask, std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      mask.readParameters(in);
      in.close();
   }

   // Determine if we get expected result when running checkCompatibility.
   // Function accepts a boolean indicating whether we expect it to throw 
   // an error or not, and returns a boolean indicating whether the 
   // function demonstrated the expected behavior.
   template <int D>
   bool checkCheckCompatibility(MaskGenFilm<D>& mask, bool expectError)
   {
      bool pass = true;
      if (expectError) {
         try {
            mask.checkCompatibility();
            // This is expected to fail. If it succeeds, the test fails.
            pass = false;
         } catch (Exception& e) {
            Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                        << std::endl;
         }
      } else {
         try {
            mask.checkCompatibility();
            // This should succeed. If not, the test fails. 
         } catch (Exception& e) {
            pass = false;
         }
      }
      return pass;
   }

};

TEST_BEGIN(MaskGenFilmTest)
TEST_ADD(MaskGenFilmTest, testConstructor)
TEST_ADD(MaskGenFilmTest, testReadParameters)
TEST_ADD(MaskGenFilmTest, testCheckCompatibility)
TEST_ADD(MaskGenFilmTest, testSetup)
TEST_ADD(MaskGenFilmTest, testUpdate)
TEST_ADD(MaskGenFilmTest, testStressTerm)
TEST_ADD(MaskGenFilmTest, testModifyStress)
TEST_END(MaskGenFilmTest)

#endif