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
      std::ifstream in;
      sweepPtr->readParameters(in);
   }

   void testParameters() 
   {
      printMethod(TEST_FUNC);

      // set up system
      System<3> system;
      SweepTest::SetUpSystem(system);
      
      // set up LinearSweep<3>::Parameter object 
      LinearSweep<3> ls(system);
      LinearSweep<3>::Parameter p;

      // open test input file
      std::ifstream in;
      openInputFile("in/param.test", in);
      in >> p;
      std::cout << '\n' << typeid(in).name() << '\n';
      std::cout << typeid(p).name() << '\n';
   }

   void SetUpSystem(System<3>& system)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/diblock/bcc/param.flex", in);
      system.readParam(in);
      in.close();
   }
};



TEST_BEGIN(SweepTest)
TEST_ADD(SweepTest, testConstructors)
TEST_ADD(SweepTest, testFactory)
TEST_ADD(SweepTest, testParameters)
TEST_END(BasisFieldStateTest)

#endif
