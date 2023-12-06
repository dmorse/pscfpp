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
   {  setVerbose(0); }

   void testAnalyzeChi()
   {
      printMethod(TEST_FUNC);
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      openLogFile("out/testSystem.log");
      ParamComponent::setEcho(true);

      std::ifstream in;
      openInputFile("in/param_simulator", in);
      system.readParam(in);
      in.close();

      Simulator<3> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      double chi = system.interaction().chi(0,1);
      TEST_ASSERT( fabs(system.interaction().chi(0,0)) < 1.0E-8);
      TEST_ASSERT( fabs(system.interaction().chi(1,1)) < 1.0E-8);

      DArray<double> vals = simulator.chiEvals();
      TEST_ASSERT( fabs((vals[0] - simulator.chiEval(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((vals[1] - simulator.chiEval(1))/chi) < 1.0E-8);
      TEST_ASSERT(fabs((vals[0] + chi)/chi) < 1.0E-8);
      TEST_ASSERT(fabs(vals[1]/chi) < 1.0E-8);

      DMatrix<double> vecs = simulator.chiEvecs();
      TEST_ASSERT(fabs(vecs(0,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(0,1) + 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,0) - 1.0) < 1.0E-8);
      TEST_ASSERT(fabs(vecs(1,1) - 1.0) < 1.0E-8);

      DArray<double>  sc = simulator.sc();
      TEST_ASSERT( fabs((sc[0] - simulator.sc(0))/chi) < 1.0E-8);
      TEST_ASSERT( fabs((sc[1] - simulator.sc(1))/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(0)/chi) < 1.0E-8);
      TEST_ASSERT( fabs(simulator.sc(1)/chi - 0.5) < 1.0E-8);

      #if 0
      std::cout << std::endl;
      std::cout << "vals  = " << vals[0] << "  " << vals[1] << std::endl;
      std::cout << "vec0  = " << vecs(0,0) << "  " << vecs(0,1) << std::endl;
      std::cout << "vec1  = " << vecs(1,0) << "  " << vecs(1,1) << std::endl;
      #endif

   }

};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testAnalyzeChi)
TEST_END(SimulatorTest)

#endif
