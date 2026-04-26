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

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;
using namespace Pscf::Rpc;

class AnalyzerTest : public LogFileUnitTest
{

public:

   void setUp()
   {  
      setVerbose(0); 
   }

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
      Analyzer<D>::initStatic();
      TEST_ASSERT(Analyzer<D>::baseInterval == 1);
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }
   
   void analyzeTrajectory()
   {
      System<3> system;
      initSystem(system, "in/param_system_disordered");
      BdSimulator<3> simulator(system);
      initSimulator(simulator, "in/param_BdSimulator_analyzer");
      std::string filename = filePrefix() + "in/w_dis_trajectory.rf";
      simulator.analyze(0, 50, "RGridTrajectoryReader", filename);
   }

   void testAnalyzeTrajectory()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzer.log");
      
      analyzeTrajectory();
   }

   void testFourthOrderParameter()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzer.log");
      
      analyzeTrajectory();
      
      std::string filename = filePrefix() + "out/fourthOrder_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/fourthOrder_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of fourthOrder parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
     
      double ref = 3.5756726e-01;
      //double diff = fabs(3.5756726e-01 - value);
      double diff = fabs(value - ref);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }
   
   void testMaxOrderParameter()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMaxOrderParameter.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/maxOrder_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/maxOrder_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of maxOrder parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
      
      double ref = 2.6173130e-02;
      double diff = fabs(ref - value);
      //double diff = fabs(2.6173130e-02 - value);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }

   void testPerturbationDerivative()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testPerturbationDerivative.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/perturbationDerivative_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/perturbationDerivative_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
      
      double ref = 4.3812680e+03;
      //double diff = fabs(4.3812680e+03 - value);
      double diff = fabs(ref - value);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }

   void testConcentrationDerivative()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testConcentrationDerivative.log");
      analyzeTrajectory();

      std::string filename = filePrefix() 
	                   + "out/concentrationDerivative_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout 
          << "Error: Could not open file out/concentrationDerivative_analyzer.ave" 
          << std::endl;

      }
      
      // Obtain the average value of parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
      
      double ref = 1.7407118e+01;
      //double diff = fabs(1.7407118e+01 - value);
      double diff = fabs(ref - value);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }

   void testChiDerivative()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testChiDerivative.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/chiDerivative_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/chiDerivative_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
      
      double ref = 8.2574275e+02;
      //double diff = fabs(8.2574275e+02 - value);
      double diff = fabs(ref - value);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }


   void testHamiltonianAnalyzer()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testHamiltonianAnalyzer.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/hamiltonian_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/hamiltonian_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of each component of Hamiltonian
      std::string line, line2, line3;
      double ideal; double field; double total;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x >> ideal;
      double idealDiff = fabs(-4.5215543e+03 - ideal);
      TEST_ASSERT(idealDiff < 1.0E-4);

      std::getline(file, line2);
      std::istringstream iss2(line2);
      iss2 >> x >> field;
      double fieldDiff = fabs(1.1596217e+04 - field);
      TEST_ASSERT(fieldDiff < 1.0E-4);

      std::getline(file, line3);
      std::istringstream iss3(line3);
      iss3 >> x >> total;
      double totalDiff = fabs(3.7887118e+03 - total);
      TEST_ASSERT(totalDiff < 1.0E-4);
      
   }

   void testBinaryStructureFactorGrid()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testBinaryStructureFactorGrid.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/binaryStructureFactorGrid_analyzer";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/binaryStructureFactorGrid_analyzer" 
                  << std::endl;

      }
      
      // Obtain the first three lines of q and S(q)
      std::string line, line2, line3;
      double q; double Sq;
      std::string x;
      std::getline(file, line);
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> q >> Sq;
      double qDiff = fabs(0.0 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      double SqDiff = fabs(-1.13103186e+00 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      std::getline(file, line2);
      std::istringstream iss2(line2);
      iss2 >> q >> Sq;
      qDiff = fabs(1.90978277e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = fabs(-5.13521218e-01 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);

      std::getline(file, line3);
      std::istringstream iss3(line3);
      iss3 >> q >> Sq;
      qDiff = fabs(2.70084069e+00 - q);
      TEST_ASSERT(qDiff < 1.0E-4);
      SqDiff = fabs(2.32910285e+00 - Sq);
      TEST_ASSERT(SqDiff < 1.0E-4);
      
   }

   void testBoxLengthDerivative()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testBoxLengthDerivative.log");
      analyzeTrajectory();

      std::string filename = filePrefix() + "out/boxLengthDerivative_analyzer.ave";
      std::ifstream file(filename);
      if (!file.is_open()) {
        std::cout << "Error: Could not open file out/boxLengthDerivative_analyzer.ave" 
                  << std::endl;

      }
      
      // Obtain the average value of parameter
      std::string line;
      double value;
      std::string x;
      std::getline(file, line);
      std::istringstream iss(line);
      iss >> x  >> x >> value;
      
      //double ref = 4.2481444e+03; // old, < v1.3.4
      double ref = 4.234712800000e+03; // new in v1.3.4
      //double diff = fabs(4.2481444e+03 - value);
      double diff = fabs(ref - value);
      if (verbose() > 0) {
         std::cout << "\nAverage value = " << Dbl(value, 20, 12);
         std::cout << "\nAverage diff  = " << Dbl(diff, 20, 12);
      }
      TEST_ASSERT(diff < 1.0E-4);
   }

};

TEST_BEGIN(AnalyzerTest)
TEST_ADD(AnalyzerTest, testAnalyzeTrajectory)
TEST_ADD(AnalyzerTest, testFourthOrderParameter)
TEST_ADD(AnalyzerTest, testMaxOrderParameter)
TEST_ADD(AnalyzerTest, testPerturbationDerivative)
TEST_ADD(AnalyzerTest, testConcentrationDerivative)
TEST_ADD(AnalyzerTest, testChiDerivative)
TEST_ADD(AnalyzerTest, testHamiltonianAnalyzer)
TEST_ADD(AnalyzerTest, testBinaryStructureFactorGrid)
TEST_ADD(AnalyzerTest, testBoxLengthDerivative)
TEST_END(AnalyzerTest)

#endif
