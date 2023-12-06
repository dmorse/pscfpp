#ifndef PSPC_SIMULATOR_TEST_H
#define PSPC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/System.h>
#include <pspc/simulate/Simulator.h>

#include <prdc/cpu/RFieldComparison.h>
#include <prdc/crystal/BFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Pspc;

class SimulatorTest : public LogFileUnitTest
{

   System<3> system;

public:


   SimulatorTest()
    : system()
   {}

   void setUp()
   {  
      setVerbose(0); 
   }

   void testInitialize()
   {
      printMethod(TEST_FUNC);
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testSystem.log");
      ParamComponent::setEcho(true);

      std::ifstream in;
      openInputFile("in/param_simulator", in);
      system.readParam(in);

      Simulator<3> simulator(system);
      simulator.readParameters(in);
      simulator.analyzeChi();
      in.close();
   }

};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testInitialize)
TEST_END(SimulatorTest)

#endif
