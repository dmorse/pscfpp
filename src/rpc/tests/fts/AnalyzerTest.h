#ifndef RPC_ANALYZER_TEST_H
#define RPC_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/analyzer/AnalyzerManager.h>

#include <util/format/Dbl.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class AnalyzerTest : public LogFileUnitTest
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

   template <int D>
   void initSimulator(Rp::BdSimulator<D,CPT>& simulator, std::string filename)
   {
      Rp::Analyzer<D,CPT>::initStatic();
      int baseInterval = Rp::Analyzer<D,CPT>::baseInterval;
      TEST_ASSERT(baseInterval == 1);
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   void analyzeTrajectory()
   {
      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_disordered");
      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/param_BdSimulator_analyzer");
      std::string filename = filePrefix() + "in/w_dis_trajectory.rf";
      simulator.analyze(0, 50, "RGridTrajectoryReader", filename);
   }

   bool checkValue(std::string label, 
                   double value, double ref, 
                   double epsilon = 1.0E-4)
   {
      double diff = std::abs(value - ref);
      bool success = (diff < epsilon);
      if (verbose() > 0 || !success) {
         std::cout << "\nlabel         = " << label;
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      return success;
   }

   bool checkAverage(std::string filename, 
                     double ref, 
                     double epsilon = 1.0E-4)
   {
      // Open file containing an average value
      std::ifstream file;
      openInputFile(filename, file);

      // Read the average value 
      double value;
      std::string x;
      file >> x >> x >> value;
      file.close();

      // Compare to reference value, return success value
      return checkValue(filename, value, ref, epsilon);
   }

   void testAnalyzeTrajectory()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzer.log");
      analyzeTrajectory();

      std::string filename;
      double ref;

      filename = "out/fourthOrder.ave";
      ref = 3.5756726e-01;
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/maxOrder.ave";
      ref = 2.6173130e-02;
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/perturbationDerivative.ave";
      ref = 4.3812680e+03;
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/concentrationDerivative.ave";
      ref = 1.7407118e+01;
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/chiDerivative.ave";
      ref = 8.2574275e+02;
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/cubicLengthDerivative.ave";
      //ref = 4.2481444e+03;      // old, before v1.3.4
      //ref = 4.234712800000e+03; // bugfix, without perturbationHamiltonian
      ref = 1.238404900000e+03;   // new in v1.3.4
      TEST_ASSERT(checkAverage(filename, ref));

      checkHamiltonianAnalyzer();
      checkBinaryStructureFactor();
   }

   void checkHamiltonianAnalyzer()
   {
      // Open data file
      std::string filename = "out/hamiltonian.ave";
      std::ifstream file;
      openInputFile(filename, file);

      // Read average values 
      std::string str;
      double ideal, field, total, error;
      file >> str >> ideal >> str >> error;
      file >> str >> field >> str >> error;
      file >> str >> total >> str >> error;
      file.close();
     
      double idealRef = -4.5215543e+03;
      TEST_ASSERT(checkValue("H_ideal", ideal, idealRef, 1.0E-3));

      double fieldRef = 1.1596217e+04;
      TEST_ASSERT(checkValue("H_field", field, fieldRef, 1.0E-3));

      double totalRef = 3.7887118e+03;
      TEST_ASSERT(checkValue("H_total", total, totalRef, 1.0E-3));
   }

   #if 0
   void checkBinaryStructureFactorGrid()
   {
      // Open data file
      std::string filename = "out/binaryStructureFactorGrid";
      std::ifstream file;
      openInputFile(filename, file);

      // Read and discard line with column labels
      std::string line;
      std::getline(file, line);

      // Obtain the first three lines of q and S(q)
      double q, Sq;

      file >> q >> Sq;
      double qDiff = fabs(0.0 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      double SqDiff = fabs(-1.13103186e+00 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      file >> q >> Sq;
      qDiff = fabs(1.90978277e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = fabs(-5.13521218e-01 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      file >> q >> Sq;
      qDiff = fabs(2.70084069e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = fabs(2.32910285e+00 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      file.close();
   }
   #endif

   void checkBinaryStructureFactor()
   {
      // Open data file
      std::string filename = "out/binaryStructureFactor_ave";
      std::ifstream file;
      openInputFile(filename, file);

      // Read and discard line with column labels
      // std::string line;
      // std::getline(file, line);

      // Obtain the first three lines of q and S(q)
      double q, Sq;

      file >> q >> Sq;
      double qDiff = std::abs(0.0 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      double SqDiff = std::abs(-1.13103186e+00 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      file >> q >> Sq;
      qDiff = std::abs(1.90978277e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = std::abs(-1.31782364E-01 - Sq);     // for proper unequal weighting
      //SqDiff = std::abs(-5.13521218e-01 - Sq);   // for equally weighted waves
      TEST_ASSERT(SqDiff < 1.0E-4);

      file >> q >> Sq;
      qDiff = std::abs(2.70084069e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = std::abs(2.76882867e+00 - Sq); // for proper unequal weighting
      //SqDiff = fabs(2.32910285e+00 - Sq);   // for equally weighted waves
      TEST_ASSERT(SqDiff < 1.0E-4);

      file.close();
   }

};

TEST_BEGIN(AnalyzerTest)
TEST_ADD(AnalyzerTest, testAnalyzeTrajectory)
TEST_END(AnalyzerTest)

#endif
