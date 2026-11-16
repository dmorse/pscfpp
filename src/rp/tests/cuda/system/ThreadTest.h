#ifndef RPG_THREAD_TEST_H
#define RPG_THREAD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Polymer.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;

class ThreadTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   template <int D>
   void initSystem(Rp::System<D,CUT>& system, std::string filename)
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
      Rp::System<1,CUT> system;
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
      Rp::System<1,CUT> system;
      initSystem(system, "in/thread/param_system_1D");
      system.w().readBasis("in/thread/omega.in");
      system.compute();
      Rp::Polymer<1,CUT> const & polymer = system.mixture().polymer(0);
      double q; double lnq;
      q = polymer.q();
      lnq = log(q);
      double diff = fabs(-5.309150997608 - lnq);
      TEST_ASSERT(diff < 1.0E-3);
      Log::file() << lnq << std::endl; 
   }

   /*
   * Read a parameter file with a contour step for which the number of
   * contour steps of a block does not fit an int.
   */
   void testReadParamTinyDs()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testReadParamTinyDs.log");

      // Parameter file in/thread/param_system_1D, with ds = 1.0E-10
      // instead of 0.01 (5.6E+9 contour steps for the first block)
      std::ifstream file;
      openInputFile("in/thread/param_system_1D", file);
      std::ostringstream buffer;
      buffer << file.rdbuf();
      file.close();
      std::string text = buffer.str();
      std::size_t pos = text.find("0.01", text.find("ds"));
      TEST_ASSERT(pos != std::string::npos);
      text.replace(pos, 4, "1.0E-10");

      Rp::System<1,CUT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(text);
      std::string message;
      try {
         system.readParam(in);
      } catch (Exception& e) {
         message = e.message();
      }
      TEST_ASSERT(message.find("too small") != std::string::npos);
   }

};

TEST_BEGIN(ThreadTest)
TEST_ADD(ThreadTest, testComputeFreeEnergyThread)
TEST_ADD(ThreadTest, testComputelnqThread)
TEST_ADD(ThreadTest, testReadParamTinyDs)
TEST_END(ThreadTest)

#endif
