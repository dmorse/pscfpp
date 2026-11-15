#ifndef RPC_FILM_FIELD_GEN_MASK_TEST_H
#define RPC_FILM_FIELD_GEN_MASK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/environment/FilmFieldGenMask.h>
#include <rp/system/System.h>
#include <rp/scft/iterator/Iterator.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/WFields.h>
#include <rp/field/Mask.h>
#include <rp/environment/FieldGenerator.h>
#include <rp/environment/Environment.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/BFieldComparison.h>

#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <util/math/Constants.h>
#include <util/containers/FSArray.h>

#include <fstream>
#include <string>
#include <sstream>
#include <initializer_list>

using namespace Util;
using namespace Pscf;

class FilmFieldGenMaskTest : public UnitTest
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
      openLogFile("out/maskTestConstructor.log");
      Rp::System<1,CPT> system;
      Rp::FilmFieldGenMask<1,CPT> mask(system);
   }

   void testReadParameters() // test FilmFieldGenMaskBase::readParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestReadParameters.log");

      // Set up mask from file
      Rp::System<2,CPT> system;
      Rp::FilmFieldGenMask<2,CPT> mask(system);

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
      TEST_ASSERT(mask.type() == Rp::FieldGenerator::Mask);
   }

   void testCheckCompatibility() // test FieldGenerator::checkCompatibility()
   {
      // This test indirectly tests various methods that are used by 
      // checkCompatibility. In the Iterator class, isFlexible, 
      // nFlexibleParams, and flexibleParams are tested. In FilmFieldGenMask and
      // its parent classes, the method checkCompatibility calls the methods
      // checkSpaceGroup, setFlexibleParams, and checkLatticeVectors. These
      // methods are tested by testing the compatibility of systems that are
      // incompatible for various reasons, as well as several systems that
      // are compatible.

      printMethod(TEST_FUNC);
      openLogFile("out/maskTestCheckCompatibility.log");

      // Set up 1D mask with a compatible system and check compatibility
      Rp::System<1,CPT> system1;
      createSystem(system1, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.9);
      //system1.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system1.setUnitCell(parameters);
      
      Rp::FilmFieldGenMask<1,CPT> mask1(system1);
      createFilmFieldGenMask(mask1, "in/filmMask1");
      
      Log::file() << "Testing system 1:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask1,false));
      TEST_ASSERT(mask1.needsUpdate());
      TEST_ASSERT(system1.iterator().isFlexible() == true);
      TEST_ASSERT(system1.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system1.iterator().flexibleParams()[0] == true);


      // Set up 2D mask with a compatible system and check compatibility
      Rp::System<2,CPT> system2;
      createSystem(system2, "in/system2D_1");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(3.0);
      //system2.setUnitCell(UnitCell<2>::Rectangular, parameters);
      system2.setUnitCell(parameters);

      Rp::FilmFieldGenMask<2,CPT> mask2(system2);
      createFilmFieldGenMask(mask2, "in/filmMask2");

      Log::file() << "Testing system 2:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask2,false));
      TEST_ASSERT(system2.iterator().isFlexible() == true);
      TEST_ASSERT(system2.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system2.iterator().flexibleParams()[0] == true);
      TEST_ASSERT(system2.iterator().flexibleParams()[1] == false);

      // Set up 2D mask with an incompatible system and check compatibility
      Rp::System<2,CPT> system3;
      createSystem(system3, "in/system2D_2");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      //system3.setUnitCell(UnitCell<2>::Square, parameters);
      system3.setUnitCell(parameters);

      Rp::FilmFieldGenMask<2,CPT> mask3(system3);
      createFilmFieldGenMask(mask3, "in/filmMask2");

      Log::file() << "Testing system 3:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask3,true));


      // Set up 3D mask with a compatible system and check compatibility
      Rp::System<3,CPT> system4;
      createSystem(system4, "in/system3D_1");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(4.2);
      //system4.setUnitCell(UnitCell<3>::Tetragonal, parameters);
      system4.setUnitCell(parameters);

      Rp::FilmFieldGenMask<3,CPT> mask4(system4);
      createFilmFieldGenMask(mask4, "in/filmMask3");

      Log::file() << "Testing system 4:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask4,false));
      TEST_ASSERT(system4.iterator().isFlexible() == true);
      TEST_ASSERT(system4.iterator().nFlexibleParams() == 1);
      TEST_ASSERT(system4.iterator().flexibleParams()[0] == true);
      TEST_ASSERT(system4.iterator().flexibleParams()[1] == false);
   

      // Set up another 3D mask with a compatible system
      Rp::System<3,CPT> system5;
      createSystem(system5, "in/system3D_2");

      // Set unit cell parameters
      parameters.clear();
      parameters.append(2.0);
      parameters.append(3.0);
      parameters.append(4.0);
      parameters.append(1.57079632679489662); // pi/2
      //system5.setUnitCell(UnitCell<3>::Monoclinic, parameters);
      system5.setUnitCell(parameters);

      Rp::FilmFieldGenMask<3,CPT> mask5(system5);
      createFilmFieldGenMask(mask5, "in/filmMask3");

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
      //system5.setUnitCell(UnitCell<3>::Monoclinic, parameters);
      system5.setUnitCell(parameters);
      Log::file() << "Testing system 6:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask5,true));

      // Try a different mask with normalVecId == 1, which is incompatible
      // due to the space group
      Rp::FilmFieldGenMask<3,CPT> mask6(system5);
      createFilmFieldGenMask(mask6, "in/filmMask2");
      Log::file() << "Testing system 7:" << std::endl;
      TEST_ASSERT(checkCheckCompatibility(mask6,true));

   }

   void testGenerate() // test FieldGenerator::generate()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestGenerate.log");

      // Set up 2D mask with a compatible system
      Rp::System<2,CPT> system;
      createSystem(system, "in/system2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(1.6353661975);
      parameters.append(2.0);
      //system.setUnitCell(UnitCell<2>::Rectangular, parameters);
      system.setUnitCell(parameters);

      // Set up mask
      Rp::FilmFieldGenMask<2,CPT> mask(system);
      createFilmFieldGenMask(mask, "in/filmMask2");
      Log::file() << "1" << std::endl;
      mask.generate();
      Log::file() << "2" << std::endl;
      TEST_ASSERT(!mask.needsUpdate());
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

   void testRegenerate() // test generate() after parameter changes
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestRegenerate.log");

      // Set up 1D mask with a compatible system
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.9);
      //system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system.setUnitCell(parameters);

      // Set up mask
      Rp::FilmFieldGenMask<1,CPT> mask(system);
      createFilmFieldGenMask(mask, "in/filmMask1");
      mask.generate();
      TEST_ASSERT(!mask.needsUpdate());

      // Change lattice parameter and update
      parameters[0] = 3.0;
      //system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system.setUnitCell(parameters);
      TEST_ASSERT(mask.needsUpdate());
      mask.generate();

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

   void testStress() // test FilmFieldGenMask::stress
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestStress.log");

      // Set up 2D mask with a compatible system
      Rp::System<2,CPT> system;
      createSystem(system, "in/system2D_1");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(1.6353661975);
      parameters.append(2.0);
      //system.setUnitCell(UnitCell<2>::Rectangular, parameters);
      system.setUnitCell(parameters);

      // Set up mask
      Rp::FilmFieldGenMask<2,CPT> mask(system);
      createFilmFieldGenMask(mask, "in/filmMask2");
      mask.generate();

      // Read w field and solve MDEs, so system can calculate fHelmholtz
      system.w().readBasis("in/wIn2D.bf");
      system.compute();

      // Call stress and check that the result is correct
      TEST_ASSERT(eq(mask.stress(0),0.0));
      TEST_ASSERT(std::abs(mask.stress(1) - 0.4211949121) < 1e-5);
   }

   void testModifyStress() // test FilmFieldGenMask::modifyStress
   {
      printMethod(TEST_FUNC);
      openLogFile("out/maskTestModifyStress.log");

      // Set up 1D mask with a compatible system
      Rp::System<1,CPT> system;
      createSystem(system, "in/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      //system.setUnitCell(UnitCell<1>::Lamellar, parameters);
      system.setUnitCell(parameters);

      // Set up mask
      Rp::FilmFieldGenMask<1,CPT> mask(system);
      createFilmFieldGenMask(mask, "in/filmMask1");
      mask.generate();

      // Read w field and solve MDEs, so system can calculate fHelmholtz
      system.w().readBasis("in/wIn1D_2.bf");
      system.compute();

      // Call modifyStress with an arbitrary input stress value 
      // and test that the result is correct
      TEST_ASSERT(eq(mask.modifyStress(0, 5.0),8.94005073083));
   }

   // Read parameter file to create a System object
   template <int D>
   void createSystem(Rp::System<D,CPT>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   // Read parameter file section to create a FilmFieldGenMask object
   template <int D>
   void createFilmFieldGenMask(Rp::FilmFieldGenMask<D,CPT>& mask, std::string fname)
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
   bool checkCheckCompatibility(Rp::FilmFieldGenMask<D,CPT>& mask, bool expectError)
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


   /*
   * Generate the fields of a film (group p_1 or P_1) with walls normal
   * to lattice vector nv and the given unit cell parameters. Returns
   * true iff an Exception is thrown. If flexible, all cell parameters
   * are flexible on input, and flexString() then gives the flexible
   * parameters chosen for the film.
   */
   template <int D>
   bool filmThrows(std::string const & lattice, int nv,
                   std::initializer_list<double> cell,
                   bool flexible = false)
   {
      std::ostringstream param;
      param << "System{\n"
            << "  Mixture{\n"
            << "    nMonomer  2\n"
            << "    monomers[  1.0  1.0  ]\n"
            << "    nPolymer  1\n"
            << "    Polymer{\n"
            << "      type    linear\n"
            << "      nBlock  2\n"
            << "      blocks[\n  0  0.5\n  1  0.5\n  ]\n"
            << "      phi     1.0\n"
            << "    }\n"
            << "    ds   0.1\n"
            << "  }\n"
            << "  Interaction{\n"
            << "    chi(\n  1   0   10.0\n    )\n"
            << "  }\n"
            << "  Domain{\n"
            << "    mesh       " << ((D == 2) ? "4 4" : "4 4 4") << "\n"
            << "    lattice    " << lattice << "\n"
            << "    groupName  " << ((D == 2) ? "p_1" : "P_1") << "\n"
            << "  }\n"
            << "  FilmEnvironment{\n"
            << "    normalVecId        " << nv << "\n"
            << "    interfaceThickness 0.2\n"
            << "    excludedThickness  0.4\n"
            << "    chiBottom[  0.0  0.0  ]\n"
            << "    chiTop[  0.0  0.0  ]\n"
            << "  }\n"
            << "  AmIteratorBasis{\n"
            << "    epsilon      1.0e-8\n"
            << "    maxItr       1\n"
            << "    maxHist      1\n"
            << "    isFlexible   " << (flexible ? 1 : 0) << "\n"
            << "  }\n"
            << "}\n";

      Rp::System<D,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(param.str());
      system.readParam(in);
      FSArray<double, 6> parameters;
      for (double p : cell) parameters.append(p);
      system.setUnitCell(parameters);

      bool thrown = false;
      try {
         system.environment().generate();
      } catch (Exception& e) {
         thrown = true;
      }
      flexible_ = system.iterator().flexibleParams();
      return thrown;
   }

   // Flexible parameters after the last call of filmThrows, e.g. "011010"
   std::string flexString() const
   {
      std::string s;
      for (int i = 0; i < flexible_.size(); ++i) {
         s += (flexible_[i] ? '1' : '0');
      }
      return s;
   }

   /*
   * Angles on either side of 90 degrees are rejected at any length
   * scale, and right angles rounded to the 11 significant digits of
   * field file headers are accepted at any length scale (issue #292).
   */
   void testCheckOrthogonality()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmFieldGenMaskTestCheckOrthogonality.log");
      double pi = Constants::Pi;
      double up = 1.5707963268;    // cos = -5.1e-12
      double down = 1.5707963267;  // cos = +9.5e-11

      // Obtuse angles between the wall normal and a parallel vector
      TEST_ASSERT(filmThrows<3>("triclinic", 1, {2, 2.5, 2, 0, 0.25, 2.0}));
      TEST_ASSERT(filmThrows<3>("triclinic", 1,
                                {2, 2.5, 2, -0.5*pi, 0.3, up}));
      TEST_ASSERT(filmThrows<3>("monoclinic", 2, {2, 2, 2.5, 1.75}));

      // An angle 1e-6 from 90 degrees in a small cell
      TEST_ASSERT(filmThrows<2>("oblique", 0,
                                {0.05, 0.0625, 0.5*pi + 1.0E-6}));

      // Rounded right angles in large cells
      TEST_ASSERT(!filmThrows<2>("oblique", 0, {100, 125, up}));
      TEST_ASSERT(!filmThrows<3>("triclinic", 0,
                                 {50, 62.5, 50, up, 0.25, up}));
      TEST_ASSERT(!filmThrows<3>("triclinic", 1,
                                 {50, 62.5, 50, 0, 0.25, down}));
   }

   /*
   * The flexible angle of a triclinic film (a, b, c, phi, theta, gamma)
   * is theta for walls normal to a (phi = pi/2) or b (phi = 0), and
   * gamma for walls normal to c. The length normal to the walls is
   * rigid.
   */
   void testTriclinicFlexibleAngle()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmFieldGenMaskTestTriclinicFlexibleAngle.log");
      double pi = Constants::Pi;
      double right = 0.5*pi;

      TEST_ASSERT(!filmThrows<3>("triclinic", 0,
                                 {2, 2.2, 2.4, right, 0.2, right}, true));
      TEST_ASSERT(flexString() == "011010");
      TEST_ASSERT(!filmThrows<3>("triclinic", 1,
                                 {2, 2.2, 2.4, 0, 0.2, right}, true));
      TEST_ASSERT(flexString() == "101010");
      TEST_ASSERT(!filmThrows<3>("triclinic", 2,
                                 {2, 2.2, 2.4, 0.3, 0, 1.9}, true));
      TEST_ASSERT(flexString() == "110001");

      // At theta = 0 any phi passes the orthogonality check, but theta
      // is the in-plane angle only for the values of phi above.
      TEST_ASSERT(filmThrows<3>("triclinic", 0,
                                {2, 2.2, 2.4, 0, 0, right}, true));
      TEST_ASSERT(filmThrows<3>("triclinic", 1,
                                {2, 2.2, 2.4, right, 0, right}, true));
      TEST_ASSERT(!filmThrows<3>("triclinic", 1,
                                 {2, 2.2, 2.4, right, 0, right}));
   }

   // Flexible parameters after the last call of filmThrows
   FSArray<bool, 6> flexible_;

};

TEST_BEGIN(FilmFieldGenMaskTest)
TEST_ADD(FilmFieldGenMaskTest, testConstructor)
TEST_ADD(FilmFieldGenMaskTest, testReadParameters)
TEST_ADD(FilmFieldGenMaskTest, testCheckCompatibility)
TEST_ADD(FilmFieldGenMaskTest, testGenerate)
TEST_ADD(FilmFieldGenMaskTest, testRegenerate)
TEST_ADD(FilmFieldGenMaskTest, testStress)
TEST_ADD(FilmFieldGenMaskTest, testModifyStress)
TEST_ADD(FilmFieldGenMaskTest, testCheckOrthogonality)
TEST_ADD(FilmFieldGenMaskTest, testTriclinicFlexibleAngle)
TEST_END(FilmFieldGenMaskTest)

#endif

