#ifndef RPC_AVERAGE_ANALYZER_TEST_H
#define RPC_AVERAGE_ANALYZER_TEST_H

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

class AverageAnalyzerTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }

   /*
   * Return the number of nonempty lines in a file.
   */
   int countLines(std::string filename)
   {
      std::ifstream file;
      openInputFile(filename, file);
      std::string line;
      int n = 0;
      while (std::getline(file, line)) {
         if (line.find_first_not_of(" \t\r") != std::string::npos) ++n;
      }
      return n;
   }

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
   * Two SIMULATE runs with averaging analyzers that write a data file
   * for every sample: an AverageListAnalyzer (HamiltonianAnalyzer) and
   * an AverageAnalyzer (MaxOrderParameter). The second run must succeed,
   * and both data files must hold the samples of that run.
   */
   void testSimulateTwice()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSimulateTwice.log");

      Rp::System<1,CPT> system;
      initSystem1D(system,
         "BdSimulator{ seed 1234 LMBdStep{ mobility 0.001 }"
         " LrAmCompressor{ epsilon 1.0e-8 maxItr 100 maxHist 20 }"
         " AnalyzerManager{\n"
         " HamiltonianAnalyzer{\n"
         " interval 1"
         " outputFileName out/twice_hamiltonian nSamplePerOutput 1 }\n"
         " MaxOrderParameter{\n"
         " interval 1"
         " outputFileName out/twice_maxOrder nSamplePerOutput 1 }\n"
         " }\n"
         " }");

      bool success = true;
      try {
         system.simulate(3);
         system.simulate(3);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);

      // Steps 0 to 3 of the second run, one line per sample
      TEST_ASSERT(countLines("out/twice_hamiltonian.dat") == 4);
      TEST_ASSERT(countLines("out/twice_maxOrder.dat") == 4);
   }

};

TEST_BEGIN(AverageAnalyzerTest)
TEST_ADD(AverageAnalyzerTest, testSimulateTwice)
TEST_END(AverageAnalyzerTest)

#endif
