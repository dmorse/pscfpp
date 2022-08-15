#ifndef PSPC_FILM_ITERATOR_TEST_H
#define PSPC_FILM_ITERATOR_TEST_H

#include "test/UnitTest.h"
#include "test/UnitTestRunner.h"
#include "util/misc/Exception.h"
#include "pscf/crystal/UnitCell.h"
#include "pspc/iterator/FilmIterator.h"
#include "pspc/iterator/AmIterator.h"
#include "pspc/field/BFieldComparison.h"
#include "pspc/field/RFieldComparison.h"
#include "pspc/field/FieldIo.h"
#include "pspc/System.h"
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pspc;

class FilmIteratorTest : public UnitTest
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
      System<3> system;
      FilmIterator<3, AmIterator<3> > iterator(system);
   }

   void testReadParameters() // test FilmIterator::readParameters()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestReadParameters.log");

      // Set up system with some data
      System<1> system;
      FilmIteratorTest::SetUpSystem(system, "in/film/system1D");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      FilmIteratorTest::SetUpFilmIterator(iterator, "in/film/film1D");

      // Check that everything was read in correctly
      TEST_ASSERT(eq(iterator.normalVecId(),0));
      TEST_ASSERT(eq(iterator.normalVecId(),0));
      TEST_ASSERT(eq(iterator.interfaceThickness(),0.08));
      TEST_ASSERT(eq(iterator.wallThickness(),0.2));
      TEST_ASSERT(eq(iterator.chi(0,0),3.0));
      TEST_ASSERT(eq(iterator.chi(1,0),0.0));
      TEST_ASSERT(eq(iterator.chi(0,1),0.0));
      TEST_ASSERT(eq(iterator.chi(1,1),4.0));
      TEST_ASSERT(iterator.isFlexible());
   }

   void testSetup() // test FilmIterator::setup()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestSetup.log");
      
      // Set up system with some data
      System<1> system;
      FilmIteratorTest::SetUpSystem(system, "in/film/system1D");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      FilmIteratorTest::SetUpFilmIterator(iterator, "in/film/film1D");

      // Run the setup function
      iterator.setup();
      
      // Check that the homogeneous components of the wall
      // and the blocks were adjusted correctly
      TEST_ASSERT(eq(iterator.phi(),9.0710679136e-02));
      TEST_ASSERT(eq(system.mixture().polymer(0).phi(),9.092893209e-01));

      // Check that the wall field files were output correctly by 
      // comparing them to the reference files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass into FieldIo functions
      DArray< DArray<double> > cFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/wall_ref.bf", 
                                       cFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(iterator.wallCField(), cFieldsCheck[0]);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);

      RField<1> cRGridCheck; // Array to store reference field
      system.fieldIo().readFieldRGrid("in/film/wall_ref.rf", 
                                      cRGridCheck, unitCell);
      DArray<RField<1> > cRGridFromIterator;
      cRGridFromIterator.allocate(1);
      cRGridFromIterator[0].allocate(system.domain().mesh().dimensions());
      DArray< DArray<double> > cFieldFromIterator;   // Put iterator cField inside a
      cFieldFromIterator.allocate(1);                // DArray so it can be passed 
      cFieldFromIterator[0] = iterator.wallCField(); // into convertBasisToRGrid
      system.fieldIo().convertBasisToRGrid(cFieldFromIterator,
                                          cRGridFromIterator);
      RFieldComparison<1> rComparison; // object to compare fields
      rComparison.compare(cRGridFromIterator[0], cRGridCheck);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << rComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(rComparison.maxDiff() < 1.0E-7);
   }

   void testCheckSpaceGroup() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup.log");

      // Set up 1D system with a correct space group and check it
      System<1> system1;
      FilmIteratorTest::SetUpSystem(system1, "in/film/system1D");
      FilmIterator<1, AmIterator<1> > iterator1(system1);
      FilmIteratorTest::SetUpFilmIterator(iterator1, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator1,false));

      // Set up 1D system with an incorrect space group and check it
      System<1> system2;
      FilmIteratorTest::SetUpSystem(system2, "in/film/system_bad_1D");
      FilmIterator<1, AmIterator<1> > iterator2(system2);
      FilmIteratorTest::SetUpFilmIterator(iterator2, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator2,true));

      // Set up 2D system with an incorrect space group and check it
      System<2> system3;
      FilmIteratorTest::SetUpSystem(system3, "in/film/system_bad_2D_2");
      FilmIterator<2, AmIterator<2> > iterator3(system3);
      FilmIteratorTest::SetUpFilmIterator(iterator3, "in/film/film2D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator3,true));

      // Set up 3D system with a correct space group and check it
      System<3> system4;
      FilmIteratorTest::SetUpSystem(system4, "in/film/system3D");
      FilmIterator<3, AmIterator<3> > iterator4(system4);
      FilmIteratorTest::SetUpFilmIterator(iterator4, "in/film/film3D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator4,false));
      TEST_ASSERT(iterator4.isFlexible()); // check that isFlexible works

      // Set up 3D system with an incorrect space group and check it
      System<3> system5;
      FilmIteratorTest::SetUpSystem(system5, "in/film/system_bad_3D_1");
      FilmIterator<3, AmIterator<3> > iterator5(system5);
      FilmIteratorTest::SetUpFilmIterator(iterator5, "in/film/film3D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator5,true));

      // Set up another 3D system with an incorrect space group and check it
      System<3> system6;
      FilmIteratorTest::SetUpSystem(system6, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator6(system6);
      FilmIteratorTest::SetUpFilmIterator(iterator6, "in/film/film3D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator6,true));
   }

   void testCheckLatticeVectors() // test FilmIterator::checkLatticeVectors()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckLatticeVectors.log");
      
      // Set up 2D system with incorrect lattice vectors and check it
      System<2> system1;
      FilmIteratorTest::SetUpSystem(system1, "in/film/system_bad_2D_1");
      FilmIterator<2, AmIterator<2> > iterator1(system1);
      FilmIteratorTest::SetUpFilmIterator(iterator1, "in/film/film2D");
      try {
         iterator1.checkLatticeVectors();
         // If above does not throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }

      // Set up 3D system with correct lattice vectors and check it
      System<3> system2;
      FilmIteratorTest::SetUpSystem(system2, "in/film/system_bad_3D_1");
      FilmIterator<3, AmIterator<3> > iterator2(system2);
      FilmIteratorTest::SetUpFilmIterator(iterator2, "in/film/film3D");
      iterator2.checkLatticeVectors(); // this should not throw an error

      // Set up 3D system with incorrect lattice vectors and check it
      System<3> system3;
      FilmIteratorTest::SetUpSystem(system3, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator3(system3);
      FilmIteratorTest::SetUpFilmIterator(iterator3, "in/film/film3D");
      try {
         iterator3.checkLatticeVectors();
         // If above doesn't throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }
   }
   
   void testFlexibleParams() // test FilmIterator::flexibleParams
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestFlexibleParams.log");
      
      // Set up 1D system and make sure flexibleParams is empty
      System<1> system1;
      FilmIteratorTest::SetUpSystem(system1, "in/film/system1D");
      FilmIterator<1, AmIterator<1> > iterator1(system1);
      FilmIteratorTest::SetUpFilmIterator(iterator1, "in/film/film1D");
      TEST_ASSERT(iterator1.flexibleParams().size() == 0);

      // Set up 2D system and make sure flexibleParams is correct
      System<2> system2;
      FilmIteratorTest::SetUpSystem(system2, "in/film/system2D");
      FilmIterator<2, AmIterator<2> > iterator2(system2);
      FilmIteratorTest::SetUpFilmIterator(iterator2, "in/film/film2D");
      TEST_ASSERT(iterator2.flexibleParams().size() == 1);
      TEST_ASSERT(iterator2.flexibleParams()[0] == 0);

      // Set up 3D tetragonal system, validate flexibleParams 
      System<3> system3;
      FilmIteratorTest::SetUpSystem(system3, "in/film/system3D");
      FilmIterator<3, AmIterator<3> > iterator3(system3);
      FilmIteratorTest::SetUpFilmIterator(iterator3, "in/film/film3D");
      TEST_ASSERT(iterator3.flexibleParams().size() == 1);
      TEST_ASSERT(iterator3.flexibleParams()[0] == 0);

      // Set up 3D monoclinic system (monoclinic), validate flexibleParams 
      System<3> system4;
      FilmIteratorTest::SetUpSystem(system4, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator4(system4);
      FilmIteratorTest::SetUpFilmIterator(iterator4, "in/film/film2D");
      // Using film2D here because it has normalVecId=1 which 
      // we want for this example
      TEST_ASSERT(iterator4.flexibleParams().size() == 3);
      TEST_ASSERT(iterator4.flexibleParams()[0] == 0);
      TEST_ASSERT(iterator4.flexibleParams()[1] == 2);
      TEST_ASSERT(iterator4.flexibleParams()[2] == 3);
   }

   void testSolve() // test FilmIterator::solve
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/filmTestSolve.log");
      
      // Set up system with some data
      System<1> system;
      FilmIteratorTest::SetUpSystem(system, "in/film/system1D");
      system.readWBasis("in/film/w_in.bf");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      FilmIteratorTest::SetUpFilmIterator(iterator, "in/film/film1D");

      // Run the setup function (has already been tested above)
      iterator.setup();

      // Run the solve function
      iterator.solve();

      // Check converged field is correct by comparing to files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_ref.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.wFieldsBasis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << bComparison.maxDiff() << "\n";
      }
      system.writeWBasis("out/w.bf");
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);
   }

   // Read parameter file to create a System object
   template <int D>
   void SetUpSystem(System<D>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }

   // Read parameter file section to create a FilmIterator object
   template <int D>
   void SetUpFilmIterator(FilmIterator<D, AmIterator<D>>& iterator, 
                          std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      iterator.readParam(in);
      in.close();
   }

   // Determine if we get expected result when running checkSpaceGroup.
   // Function accepts a boolean indicating whether we expect it to throw 
   // an error or not, and returns a boolean indicating whether the 
   // function demonstrated the expected behavior.
   template <int D>
   bool checkCheckSpaceGroup(FilmIterator<D, AmIterator<D>>& iterator, 
                             bool expectError)
   {
      bool pass = true;
      if (expectError) {
         try {
            iterator.checkSpaceGroup();
            // This is expected to fail. If it succeeds, the test fails.
            pass = false;
         } catch (Exception e) {
            Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                        << std::endl;
         }
      } else {
         try {
            iterator.checkSpaceGroup();
            // This should succeed. If not, the test fails. 
         } catch (Exception e) {
            pass = false;
         }
      }
      return pass;
   }
};

TEST_BEGIN(FilmIteratorTest)
TEST_ADD(FilmIteratorTest, testConstructor)
TEST_ADD(FilmIteratorTest, testReadParameters)
TEST_ADD(FilmIteratorTest, testSetup)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup)
TEST_ADD(FilmIteratorTest, testCheckLatticeVectors)
TEST_ADD(FilmIteratorTest, testFlexibleParams)
TEST_ADD(FilmIteratorTest, testSolve)
TEST_END(FilmIteratorTest)

#endif
