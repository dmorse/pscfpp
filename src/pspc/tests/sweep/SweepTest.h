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
      LinearSweep<3> ls2(system);
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
};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructors)
TEST_ADD(SweepTest, testFactory)
TEST_END(BasisFieldStateTest)

#endif
