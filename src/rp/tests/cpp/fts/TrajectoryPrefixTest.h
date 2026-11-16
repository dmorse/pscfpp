#ifndef RPC_TRAJECTORY_PREFIX_TEST_H
#define RPC_TRAJECTORY_PREFIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <pscf/field/cpp/RField.h>

#include <util/misc/FileMaster.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

class TrajectoryPrefixTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Read a 1D diblock melt (16 grid points) with the given simulator
   * block, and set w fields with a cosine exchange component.
   */
   void initSystem1D(Rp::System<1,CPT>& system,
                     std::string const & simulatorBlock)
   {
      Rp::Analyzer<1,CPT>::initStatic();
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.05 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar } "
         + simulatorBlock + " }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(1.5);
      system.setUnitCell(parameters);
      DArray< RField<1,CPT> > w;
      w.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < 16; ++k) {
            w[i][k] = (i == 0) ? std::cos(2.0*M_PI*k/16.0) : 0.0;
         }
      }
      system.w().setRGrid(w);
   }

   /*
   * Write a trajectory of the 1D diblock melt with frames 0 to nStep,
   * compressed only to an error of 1e-2.
   */
   void writeTrajectory(std::string const & filename, int nStep)
   {
      Rp::System<1,CPT> system;
      initSystem1D(system,
         "BdSimulator{ seed 1234 LMBdStep{ mobility 0.01 }"
         " LrAmCompressor{ epsilon 1.0e-2 maxItr 100 maxHist 20 }"
         " AnalyzerManager{\n"
         " TrajectoryWriter{\n"
         " interval 1"
         " outputFileName " + filename + " }\n"
         " }\n"
         " }");
      system.simulate(nStep);
   }

   /*
   * ANALYZE with the input prefix out/ must find a trajectory file
   * named relative to the prefix, and still accept a name that already
   * includes the prefix.
   */
   void testAnalyzeInputPrefix()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeInputPrefix.log");

      writeTrajectory("out/trajectory1D.rf", 1);
      std::string names[2] = {"trajectory1D.rf",
                              filePrefix() + "out/trajectory1D.rf"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, "BdSimulator{ AnalyzerManager{\n"
            " HamiltonianAnalyzer{\n"
            " interval 1"
            " outputFileName out/analyzePrefix }\n"
            " }\n"
            " }");
         system.fileMaster().setInputPrefix(filePrefix() + "out/");
         bool success = true;
         try {
            system.simulator().analyze(0, 1, "RGridTrajectoryReader",
                                       names[i]);
         } catch (Exception& e) {
            success = false;
         }
         TEST_ASSERT(success);
      }
   }

};

TEST_BEGIN(TrajectoryPrefixTest)
TEST_ADD(TrajectoryPrefixTest, testAnalyzeInputPrefix)
TEST_END(TrajectoryPrefixTest)

#endif
