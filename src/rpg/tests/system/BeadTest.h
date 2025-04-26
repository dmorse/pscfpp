#ifndef RPG_BEAD_TEST_H
#define RPG_BEAD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/System.h>
#include <rpg/solvers/Polymer.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class BeadTest : public LogFileUnitTest
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
   
   /*
   * Compare Helmoltze free energies to prior result.
   */
   void testComputeFreeEnergyBead()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testComputeFreeEnergyBead.log");
      System<1> system;
      initSystem(system, "in/bead/param_system_1D_N100");
      system.readWBasis("in/bead/omegaN100.in");
      system.compute();
      system.computeFreeEnergy();
      double diff = fabs(1.92537673380e-02 - system.fHelmholtz());
      TEST_ASSERT(diff < 1.0E-5);
   }
   
   /*
   * Compare lnq to prior result.
   */
   void testComputelnqBead()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testComputelnqBead.log");
      System<1> system;
      initSystem(system, "in/bead/param_system_1D_N100");
      system.readWBasis("in/bead/omegaN100.in");
      system.compute();
      Polymer<1> & polymer = system.mixture().polymer(0);
      double q; double lnq;
      q = polymer.q();
      lnq = log(q);
      double diff = fabs(-5.3088407702611 - lnq);
      TEST_ASSERT(diff < 1.0E-3);
   }

};

TEST_BEGIN(BeadTest)
TEST_ADD(BeadTest, testComputeFreeEnergyBead)
TEST_ADD(BeadTest, testComputelnqBead)
TEST_END(BeadTest)

#endif
