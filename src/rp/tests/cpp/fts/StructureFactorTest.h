#ifndef RPC_STRUCTURE_FACTOR_TEST_H
#define RPC_STRUCTURE_FACTOR_TEST_H

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

class StructureFactorTest : public LogFileUnitTest
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
   * SIMULATE with a BinaryStructureFactor in a rectangular cell 3.2 x
   * 3.2000000064 on a 4 x 4 mesh, in which two bunches of waves differ
   * by 1.5E-8 in |k|^2, but only by 4E-9 in |k|: the 9 bunches must be
   * written in order of increasing wavenumber.
   */
   void testStructureFactorCloseBunches()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testStructureFactorCloseBunches.log");

      Rp::Analyzer<2,CPT>::initStatic();
      Rp::System<2,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.05 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 4 4 lattice rectangular }"
         " BdSimulator{ seed 1234 LMBdStep{ mobility 0.01 }"
         " LrCompressor{ epsilon 1.0e-8 maxItr 100 }"
         " AnalyzerManager{\n"
         " BinaryStructureFactor{\n"
         " interval 1"
         " outputFileName out/closeBunches }\n"
         " }\n"
         " } }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(3.2);
      parameters.append(3.2000000064);
      system.setUnitCell(parameters);
      DArray< RField<2,CPT> > w;
      w.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < 16; ++k) {
            w[i][k] = (i == 0) ? std::cos(2.0*M_PI*(k/4)/4.0) : 0.0;
         }
      }
      system.w().setRGrid(w);

      bool success = true;
      try {
         system.simulate(1);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(countLines("out/closeBunches") == 9);
      std::ifstream file;
      openInputFile("out/closeBunches", file);
      double q, s, qOld = -1.0;
      for (int i = 0; i < 9; ++i) {
         file >> q >> s;
         TEST_ASSERT(q > qOld);
         qOld = q;
      }
   }

};

TEST_BEGIN(StructureFactorTest)
TEST_ADD(StructureFactorTest, testStructureFactorCloseBunches)
TEST_END(StructureFactorTest)

#endif
