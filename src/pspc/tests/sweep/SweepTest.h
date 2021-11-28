#ifndef PSPC_SWEEP_TEST_H
#define PSPC_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/SweepFactory.h>
#include <pspc/sweep/LinearSweep.h>

#include <fstream>


using namespace Util;
using namespace Pscf;
using namespace Pspc;

class SweepTest : public UnitTest
{

public:

   std::ofstream logFile_;

   void setUp()
   {}

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

   void testConstructors()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      LinearSweep<3> ls(system);
      SweepFactory<3> sf(system);

   }

   void testFactory() 
   {
      printMethod(TEST_FUNC);

      System<3> system;
      Sweep<3>* sweepPtr;
      SweepFactory<3> sf(system);
      
      sweepPtr = sf.factory("LinearSweep");
   }

   void testParameters() 
   {
      printMethod(TEST_FUNC);

      // set up system with some data
      System<3> system;
      SweepTest::SetUpSystem(system, "in/diblock/bcc/param.flex");
      // set up LinearSweepParameter objects 
      DArray< LinearSweepParameter<3> > ps;
      ps.allocate(4);
      for (int i = 0; i < 4; ++i) {
         ps[i].setSystem(system);
      }
      // open test input file
      std::ifstream in;
      // read in data
      openInputFile("in/param.test", in);
      for (int i = 0; i < 4; ++i) {
         in >> ps[i];
         ps[i].getInitial();
         ps[i].update(0.5);
      }
   }

   void testLinearSweepRead()
   {
      printMethod(TEST_FUNC);

      // set up system with Linear Sweep Object
      System<3> system;
      SweepTest::SetUpSystem(system, "in/param.sweep");
   }

   void testLinearSweepIterate()
   {
      printMethod(TEST_FUNC);

      // set up system with Linear Sweep Object
      System<3> system;
      SweepTest::SetUpSystem(system, "in/param.sweep");
      
      // initial guess and iterate
      system.readWBasis("in/diblock/bcc/omega.ref");
      DArray< DArray<double> > wFields_check;
      wFields_check = system.wFields();
      system.readWBasis("in/diblock/bcc/omega.in");
      system.iterate();

      //sweep
      system.sweep();

   }

   void SetUpSystem(System<3>& system, std::string fname)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile(fname, in);
      system.readParam(in);
      in.close();
   }
};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructors)
TEST_ADD(SweepTest, testFactory)
TEST_ADD(SweepTest, testParameters)
TEST_ADD(SweepTest, testLinearSweepRead)
TEST_ADD(SweepTest, testLinearSweepIterate)
TEST_END(SweepTest)

#endif
