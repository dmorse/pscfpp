#ifndef RPG_MC_SIMULATOR_TEST_H
#define RPG_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/compressor/Compressor.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class McSimulatorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }
   
   template <int D>
   void initSystem(System<D>& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

   }
   
   template <int D>
   void initSimulator(McSimulator<D>& simulator, std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }
   
   void testMcSimulate()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulate.log");
      
      System<3> system;
      initSystem(system, "in/param_system_disordered");
      
      McSimulator<3> simulator(system);
      initSimulator(simulator, "in/param_McSimulator");
      
      system.readWRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
   }

};

TEST_BEGIN(McSimulatorTest)
TEST_ADD(McSimulatorTest, testMcSimulate)
TEST_END(McSimulatorTest)

#endif
