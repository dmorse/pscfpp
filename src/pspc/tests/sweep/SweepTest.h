#ifndef PSPC_SWEEP_TEST_H
#define PSPC_SWEEP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/Sweep.h>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pspc;

class SweepTest : public UnitTest
{

public:

   std::ofstream logFile_; // Log file for logging errors? Probably?

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

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testConstructor.log");

   }
};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructor)
TEST_END(BasisFieldStateTest)

#endif
