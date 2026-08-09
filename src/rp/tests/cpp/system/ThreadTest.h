#ifndef RPC_THREAD_TEST_H
#define RPC_THREAD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Polymer.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class ThreadTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CPT>& system, std::string filename)
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
   void testComputeFreeEnergyThread()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testComputeFreeEnergyThread.log");
      Rp::System<1,CPT> system;
      initSystem(system, "in/thread/param_system_1D");
      system.w().readBasis("in/thread/omega.in");
      system.compute();
      system.scft().compute();
      double diff = fabs(1.92567500293 - system.scft().fHelmholtz());
      TEST_ASSERT(diff < 1.0E-3);
   }
   
   /*
   * Compare lnq to prior result.
   */
   void testComputelnqThread()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testComputelnqThread.log");
      Rp::System<1,CPT> system;
      initSystem(system, "in/thread/param_system_1D");
      system.w().readBasis("in/thread/omega.in");
      system.compute();
      Rp::Polymer<1,CPT> const & polymer = system.mixture().polymer(0);
      double q; double lnq;
      q = polymer.q();
      lnq = log(q);
      double diff = fabs(-5.309150997608 - lnq);
      TEST_ASSERT(diff < 1.0E-3);
   }

};

TEST_BEGIN(ThreadTest)
TEST_ADD(ThreadTest, testComputeFreeEnergyThread)
TEST_ADD(ThreadTest, testComputelnqThread)
TEST_END(ThreadTest)

#endif
