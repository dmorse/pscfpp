#ifndef RPC_ANALYZER_TEST_H
#define RPC_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/analyzer/AnalyzerManager.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

class AnalyzerTest : public LogFileUnitTest
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

   void testAnalyzeTrajectory()
   {
      printMethod(TEST_FUNC);
      
      openLogFile("out/testAnalyzer.log");
      System<3> system;
      initSystem(system, "in/param_system_disordered");
      BdSimulator<3> simulator(system);
      initSimulator(simulator, "in/param_BdSimulator_analyzer");
      std::string filename = filePrefix() + "in/w_dis_trajectory.rf";
      simulator.analyze(0, 10, "RGridTrajectoryReader", filename);
   }
   
   void testFourthOrderParameter()
   {
      printMethod(TEST_FUNC);
      std::string filename = filePrefix() + "out/fourthOrder_analyzer";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/fourthOrder" 
                  << std::endl;

      }
      std::string line;
      
      // Skip the first line
      std::getline(file, line);
      
      double chi; double fourthorder;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> chi >> fourthorder;
      double diff = fabs(0.42383968 - fourthorder);
      TEST_ASSERT(diff < 1.0E-2);
   }

};

TEST_BEGIN(AnalyzerTest)
TEST_ADD(AnalyzerTest, testAnalyzeTrajectory)
TEST_ADD(AnalyzerTest, testFourthOrderParameter)
TEST_END(AnalyzerTest)

#endif
