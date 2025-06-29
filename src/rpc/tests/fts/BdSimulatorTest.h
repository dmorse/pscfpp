#ifndef RPC_BD_SIMULATOR_TEST_H
#define RPC_BD_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

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
      
      system.w().readRGrid("in/w_dis.rf");
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
      
      system.w().readRGrid("in/w_triblock.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
   }

};

TEST_BEGIN(BdSimulatorTest)
TEST_ADD(BdSimulatorTest, testLMBdSimulateDiblocks)
TEST_ADD(BdSimulatorTest, testLMBdSimulateTriblocks)
TEST_END(BdSimulatorTest)

#endif
