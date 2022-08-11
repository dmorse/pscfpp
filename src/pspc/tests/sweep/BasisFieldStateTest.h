#ifndef PSPC_BASIS_FIELD_STATE_TEST_H
#define PSPC_BASIS_FIELD_STATE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/sweep/BasisFieldState.h>
#include <pscf/crystal/BFieldComparison.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class BasisFieldStateTest : public LogFileUnitTest
{

public:

   void setUp()
   {}

   void testConstructor()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs1(system);
      BasisFieldState<3> bfs2;
   }

   void testRead()
   {
      printMethod(TEST_FUNC);
      
      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;
   
      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in file one way
      system.readWBasis("in/bcc/omega.ref");
      // Read in file another way
      bfs.read("in/bcc/omega.ref");
      // Compare
      comparison.compare(bfs.fields(), system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);

   }

   void testWrite()
   {
      // Write tested with a read/write/read/comparison procedure
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs1(system), bfs2(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // read, write, read
      bfs1.read("in/bcc/omega.ref");
      bfs1.write("out/testBasisFieldStateWrite.ref");
      bfs2.read("out/testBasisFieldStateWrite.ref");

      // compare
      comparison.compare(bfs1.fields(),bfs2.fields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testGetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using system
      system.readWBasis("in/bcc/omega.ref");
      // get it using bfs
      bfs.getSystemState();
      // compare
      comparison.compare(bfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testSetSystemState()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs(system);
      BFieldComparison comparison;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);

      // Read in state using bfs
      bfs.read("in/bcc/omega.ref");
      // set system state
      bfs.setSystemState(true);
      // compare
      comparison.compare(bfs.fields(),system.wFields());
      // Assert small difference
      TEST_ASSERT(comparison.maxDiff() < 5.0e-7);
   }

   void testSetSystem()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      BasisFieldState<3> bfs;

      // Setup system
      BasisFieldStateTest::SetUpSystem(system);
      // Invoke setSystem
      bfs.setSystem(system);
   }

   void SetUpSystem(System<3>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/bcc/param.flex", in);
      system.readParam(in);
      in.close();
   }

};



TEST_BEGIN(BasisFieldStateTest)
TEST_ADD(BasisFieldStateTest, testConstructor)
TEST_ADD(BasisFieldStateTest, testRead)
TEST_ADD(BasisFieldStateTest, testWrite)
TEST_ADD(BasisFieldStateTest, testGetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystemState)
TEST_ADD(BasisFieldStateTest, testSetSystem)
TEST_END(BasisFieldStateTest)

#endif
