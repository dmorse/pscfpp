#ifndef RPG_BD_SIMULATOR_TEST_H
#define RPG_BD_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/fts/compressor/Compressor.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class BdSimulatorTest : public LogFileUnitTest
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
   void initSimulator(BdSimulator<D>& simulator, std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }
   
   void testLMBdSimulateDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLMBdSimulateDiblocks.log");
      
      System<3> system;
      initSystem(system, "in/param_system_disordered");
      
      BdSimulator<3> simulator(system);
      initSimulator(simulator, "in/param_BdSimulator");
      
      system.readWRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
   }

   void testLMBdSimulateTriblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLMBdSimulateTriblocks.log");
      
      System<3> system;
      initSystem(system, "in/param_system_triblock");
      
      BdSimulator<3> simulator(system);
      initSimulator(simulator, "in/param_BdSimulator");
      
      system.readWRGrid("in/w_triblock.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
   }

};

TEST_BEGIN(BdSimulatorTest)
TEST_ADD(BdSimulatorTest, testLMBdSimulateDiblocks)
TEST_ADD(BdSimulatorTest, testLMBdSimulateTriblocks)
TEST_END(BdSimulatorTest)

#endif
