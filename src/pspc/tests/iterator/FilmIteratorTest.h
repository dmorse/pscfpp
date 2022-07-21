#ifndef PSPC_FILM_ITERATOR_TEST_H
#define PSPC_FILM_ITERATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <pscf/crystal/UnitCell.h>
#include <pspc/iterator/FilmIterator.h>
#include <pspc/iterator/AmIterator.h>
#include <pspc/field/BFieldComparison.h>
#include <pspc/field/RFieldComparison.h>
#include <pspc/field/FieldIo.h>
#include <pspc/System.h>
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

      // Set up system with some data
      System<1> system;
      FilmIteratorTest::SetUpSystem(system, "in/film/system");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      std::ifstream in;
      openInputFile("in/film/film1D", in);
      iterator.readParam(in);
      in.close();

      // Check that everything was read in correctly
      TEST_ASSERT(eq(iterator.normalVec(),0));
      TEST_ASSERT(eq(iterator.normalVec(),0));
      TEST_ASSERT(eq(iterator.interfaceThickness(),0.08));
      TEST_ASSERT(eq(iterator.wallThickness(),0.2));
      TEST_ASSERT(eq(iterator.chi(0,0),3.0));
      TEST_ASSERT(eq(iterator.chi(1,0),0.0));
      TEST_ASSERT(eq(iterator.chi(0,1),0.0));
      TEST_ASSERT(eq(iterator.chi(1,1),4.0));
      TEST_ASSERT(!iterator.isFlexible());
   }

   void testSetup() // test FilmIterator::setup()
   {
      printMethod(TEST_FUNC);
      
      // Set up system with some data
      System<1> system;
      FilmIteratorTest::SetUpSystem(system, "in/film/system");

      // Set up iterator from file
      FilmIterator<1, AmIterator<1> > iterator(system);
      std::ifstream in;
      openInputFile("in/film/film1D", in);
      iterator.readParam(in);
      in.close();

      // Run the setup function
      iterator.setup();
      
      // Check that the homogeneous components of the wall
      // and the blocks were adjusted correctly
      TEST_ASSERT(eq(iterator.phi(),9.0710679136e-02));
      TEST_ASSERT(eq(system.mixture().polymer(0).phi(),9.092893209e-01));

      // Check that the wall field files were output correctly
      UnitCell<1> unitCell; // UnitCell object to pass into FieldIo functions
      DArray< DArray<double> > cFieldsCheck; // Array to store reference field
      system.fieldIo().readFieldsBasis("in/film/wall_ref.bf", cFieldsCheck, unitCell);
      BFieldComparison bComparison(0); // object to compare fields
      bComparison.compare(iterator.wallCField(), cFieldsCheck[0]);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Max error = " << bComparison.maxDiff() << "\n";
      }
      TEST_ASSERT(bComparison.maxDiff() < 1.0E-7);

      RField<1> cRGridCheck; // Array to store reference field
      system.fieldIo().readFieldRGrid("in/film/wall_ref.rf", cRGridCheck, unitCell);
      DArray<RField<1> > cRGridFromIterator;
      cRGridFromIterator.allocate(1);
      cRGridFromIterator[0].allocate(system.domain().mesh().dimensions());
      DArray< DArray<double> > cFieldFromIterator;   // Put iterator cField inside a
      cFieldFromIterator.allocate(1);                // DArray so it can be passed 
      cFieldFromIterator[0] = iterator.wallCField(); // into convertBasisToRGrid
      system.fieldIo().convertBasisToRGrid(cFieldFromIterator,cRGridFromIterator);
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
      TEST_ASSERT(1==1);
   }

   void testFlexibleParams() // test FilmIterator::flexibleParams
   {
      TEST_ASSERT(1==1);
   }

   void testSolve() // test FilmIterator::solve
   {
      TEST_ASSERT(1==1);
   }

   void SetUpSystem(System<1>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }
};

TEST_BEGIN(FilmIteratorTest)
TEST_ADD(FilmIteratorTest, testConstructor)
TEST_ADD(FilmIteratorTest, testReadParameters)
TEST_ADD(FilmIteratorTest, testSetup)
TEST_ADD(FilmIteratorTest, testCheckSpaceGroup)
TEST_ADD(FilmIteratorTest, testFlexibleParams)
TEST_ADD(FilmIteratorTest, testSolve)
TEST_END(FilmIteratorTest)

#endif
