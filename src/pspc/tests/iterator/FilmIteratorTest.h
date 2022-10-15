#ifndef PSPC_FILM_ITERATOR_TEST_H
#define PSPC_FILM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/iterator/FilmIterator.h>
#include <pspc/iterator/AmIterator.h>
#include <pspc/field/RFieldComparison.h>
#include <pspc/field/FieldIo.h>
#include <pspc/System.h>

#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/UnitCell.h>

#include <util/misc/Exception.h>

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
      System<2> system;
      FilmIteratorTest::setUpSystem(system, "in/film/system2D");

      // Set up iterator from file
      FilmIterator<2, AmIterator<2> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film2D");

      // Check that everything was read in correctly
      TEST_ASSERT(eq(iterator.normalVecId(),1));
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
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film1D");

      system.readWBasis("in/film/w_ref.bf");

      // Run the setup function
      iterator.setup();
      
      // Check that the homogeneous components of the wall
      // and the blocks were adjusted correctly
      TEST_ASSERT(eq(iterator.maskPhi(),1.0538978363e-01));
      TEST_ASSERT(eq(system.mixture().polymer(0).phi(),8.946102164e-01));

      // Check that the wall field files were output correctly by 
      // comparing them to the reference files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass into FieldIo functions
      DArray< DArray<double> > cFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/wall_ref.bf", 
                                       cFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(iterator.maskBasis(), cFieldsCheck[0]);
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

      // Put iterator cField inside a DArray so it can be passed into 
      // convertBasisToRGrid
      DArray< DArray<double> > cFieldFromIterator;   
      cFieldFromIterator.allocate(1); 
      cFieldFromIterator[0] = iterator.maskBasis(); 

      
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

   void testCheckSpaceGroup1DA() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup1DA.log");

      // Set up 1D system with a correct space group and check it
      System<1> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system1D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      double parameter = 2.9;
      parameters.append(parameter);
      system1.setUnitCell(UnitCell<1>::Lamellar, parameters);
      //system1.readWBasis("in/film/w_ref.bf");

      FilmIterator<1, AmIterator<1> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator1,false));
   }

   void testCheckSpaceGroup1DB() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup1DB.log");

      // Set up 1D system with an incorrect space group and check it
      System<1> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system_bad_1D");
      FilmIterator<1, AmIterator<1> > iterator2(system2);

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      //double parameter = 2.2;
      parameters.append(2.2);
      system2.setUnitCell(UnitCell<1>::Lamellar, parameters);

      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film1D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator2,true));
   }

   void testCheckSpaceGroup2D() // test FilmIterator::checkSpaceGroup
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup2D.log");

      // Set up 2D system with an incorrect space group and check it
      System<2> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system_bad_2D_2");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      system3.setUnitCell(UnitCell<2>::Rectangular, parameters);

      FilmIterator<2, AmIterator<2> > iterator3(system3);

      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film2D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator3,true));
   }

   void testCheckSpaceGroup3DA() 
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup3DA.log");

      // Set up 3D system with a correct space group and check it
      System<3> system4;
      FilmIteratorTest::setUpSystem(system4, "in/film/system3D");
      FilmIterator<3, AmIterator<3> > iterator4(system4);
      FilmIteratorTest::setUpFilmIterator(iterator4, "in/film/film3D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(4.2);
      system4.setUnitCell(UnitCell<3>::Tetragonal, parameters);

      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator4,false));
      TEST_ASSERT(iterator4.isFlexible()); // check that isFlexible works
   }

   void testCheckSpaceGroup3DB() 
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckSpaceGroup3DB.log");

      // Set up 3D system with an incorrect space group and check it
      System<3> system5;
      FilmIteratorTest::setUpSystem(system5, "in/film/system_bad_3D_1");
      FilmIterator<3, AmIterator<3> > iterator5(system5);
      FilmIteratorTest::setUpFilmIterator(iterator5, "in/film/film3D");

      // Set unit cell parameter
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(4.2);
      system5.setUnitCell(UnitCell<3>::Tetragonal, parameters);

      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator5,true));

      // Set up another 3D system with an incorrect space group and check it
      System<3> system6;
      FilmIteratorTest::setUpSystem(system6, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator6(system6);
      FilmIteratorTest::setUpFilmIterator(iterator6, "in/film/film3D");
      TEST_ASSERT(FilmIteratorTest::checkCheckSpaceGroup(iterator6,true));
   }

   #if 0
   void testCheckLatticeVectors() // test FilmIterator::checkLatticeVectors()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestCheckLatticeVectors.log");
      
      // Set up 2D system with incorrect lattice vectors and check it
      System<2> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system_bad_2D_1");
      FilmIterator<2, AmIterator<2> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film2D");
      try {
         iterator1.checkLatticeVectors();
         // If above does not throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception& e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }

      // Set up 3D system with correct lattice vectors and check it
      System<3> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system_bad_3D_1");
      FilmIterator<3, AmIterator<3> > iterator2(system2);
      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film3D");
      iterator2.checkLatticeVectors(); // this should not throw an error

      // Set up 3D system with incorrect lattice vectors and check it
      System<3> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator3(system3);
      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film3D");
      try {
         iterator3.checkLatticeVectors();
         // If above doesn't throw an error, then it failed this test
         TEST_ASSERT(1 == 2);
      } catch (Exception& e) {
         Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                     << std::endl;
      }
   }
   #endif
   
   void testFlexibleParams() // test FilmIterator::flexibleParams
   {
      printMethod(TEST_FUNC);
      openLogFile("out/filmTestFlexibleParams.log");
      
      // Set up 1D system and make sure flexibleParams is empty
      System<1> system1;
      FilmIteratorTest::setUpSystem(system1, "in/film/system1D");
      FilmIterator<1, AmIterator<1> > iterator1(system1);
      FilmIteratorTest::setUpFilmIterator(iterator1, "in/film/film1D");
      iterator1.setFlexibleParams();
      TEST_ASSERT(iterator1.flexibleParams().size() == 0);

      // Set up 2D system and make sure flexibleParams is correct
      System<2> system2;
      FilmIteratorTest::setUpSystem(system2, "in/film/system2D");
      FilmIterator<2, AmIterator<2> > iterator2(system2);
      FilmIteratorTest::setUpFilmIterator(iterator2, "in/film/film2D");
      iterator2.setFlexibleParams();
      TEST_ASSERT(iterator2.flexibleParams().size() == 1);
      TEST_ASSERT(iterator2.flexibleParams()[0] == 0);

      // Set up 3D tetragonal system, validate flexibleParams 
      System<3> system3;
      FilmIteratorTest::setUpSystem(system3, "in/film/system3D");
      FilmIterator<3, AmIterator<3> > iterator3(system3);
      FilmIteratorTest::setUpFilmIterator(iterator3, "in/film/film3D");
      iterator3.setFlexibleParams();
      TEST_ASSERT(iterator3.flexibleParams().size() == 1);
      TEST_ASSERT(iterator3.flexibleParams()[0] == 0);

      // Set up 3D monoclinic system (monoclinic), validate flexibleParams 
      System<3> system4;
      FilmIteratorTest::setUpSystem(system4, "in/film/system_bad_3D_2");
      FilmIterator<3, AmIterator<3> > iterator4(system4);
      FilmIteratorTest::setUpFilmIterator(iterator4, "in/film/film2D");
      // Using film2D here because it has normalVecId=1 which 
      // we want for this example
      iterator4.setFlexibleParams();
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
      FilmIteratorTest::setUpSystem(system, "in/film/system1D");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      FilmIteratorTest::setUpFilmIterator(iterator, "in/film/film1D");

      // Read initial guess
      system.readWBasis("in/film/w_in.bf");

      // Run the solve function
      iterator.solve();

      // Check converged field is correct by comparing to files in in/film
      UnitCell<1> unitCell; // UnitCell object to pass to FieldIo functions
      DArray< DArray<double> > wFieldsCheck; // Copy of reference field
      system.fieldIo().readFieldsBasis("in/film/w_ref.bf", 
                                       wFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(system.w().basis(), wFieldsCheck);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << bComparison.maxDiff() << "\n";
      }
      system.writeWBasis("out/w.bf");
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-5);
   }

   // Read parameter file to create a System object
   template <int D>
   void setUpSystem(System<D>& system, std::string fname)
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
   void setUpFilmIterator(FilmIterator<D, AmIterator<D>>& iterator, 
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
         } catch (Exception& e) {
            Log::file() << "EXCEPTION CAUGHT, expected behavior occurred" 
                        << std::endl;
         }
      } else {
         try {
            iterator.checkSpaceGroup();
            // This should succeed. If not, the test fails. 
         } catch (Exception& e) {
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
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup1DA)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup1DB)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup2D)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup3DA)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup3DB)
//TEST_ADD(FilmIteratorTest, testCheckLatticeVectors)
//TEST_ADD(FilmIteratorTest, testFlexibleParams)
//TEST_ADD(FilmIteratorTest, testSolve)
TEST_END(FilmIteratorTest)

#endif
