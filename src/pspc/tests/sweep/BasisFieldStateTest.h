#ifndef PSPC_BASISFIELDSTATE_TEST_H
#define PSPC_BASISFIELDSTATE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/BasisFieldState.h>
#include <pspc/field/BFieldComparison.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class BasisFieldStateTest : public UnitTest
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

      System<3> system;
      BasisFieldState<3> bfs(system);

   }

   void testAllocate()
   {
      // system.fileMaster().setInputPrefix(filePrefix());
      // system.fileMaster().setOutputPrefix(filePrefix());

      // std::ifstream in;
      // openInputFile("in/diblock/bcc/param.flex", in);
      // system.readParam(in);
      // in.close();

      // system.readWBasis("in/diblock/bcc/omega.ref");
   }

   void testRead()
   {

   }

   void testWrite()
   {

   }

   void testGetSystemState()
   {

   }

   void testSetSystemState()
   {

   }  

};

TEST_BEGIN(BasisFieldStateTest)

TEST_ADD(BasisFieldStateTest, testConstructor)
// TEST_ADD(BasisFieldStateTest, testAllocate)
// TEST_ADD(BasisFieldStateTest, testRead)
// TEST_ADD(BasisFieldStateTest, testWrite)
// TEST_ADD(BasisFieldStateTest, testGetSystemState)
// TEST_ADD(BasisFieldStateTest, testSetSystemState)

TEST_END(BasisFieldStateTest)

#endif
