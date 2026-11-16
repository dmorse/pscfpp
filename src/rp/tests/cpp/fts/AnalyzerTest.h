#ifndef RPC_ANALYZER_TEST_H
#define RPC_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/analyzer/BinaryChiDerivative.h>
#include <rp/fts/analyzer/CubicLengthDerivative.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/field/FieldIo.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <util/format/Dbl.h>
#include <util/misc/FileMaster.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

/*
* BinaryChiDerivative with a public compute function.
*/
class TestChiDerivative : public Rp::BinaryChiDerivative<1,CPT>
{
public:

   TestChiDerivative(Rp::Simulator<1,CPT>& simulator,
                     Rp::System<1,CPT>& system)
    : Rp::BinaryChiDerivative<1,CPT>(simulator, system)
   {}

   using Rp::BinaryChiDerivative<1,CPT>::compute;

};

/*
* CubicLengthDerivative with a public compute function.
*/
class TestLengthDerivative : public Rp::CubicLengthDerivative<3,CPT>
{
public:

   TestLengthDerivative(Rp::Simulator<3,CPT>& simulator,
                        Rp::System<3,CPT>& system)
    : Rp::CubicLengthDerivative<3,CPT>(simulator, system)
   {}

   using Rp::CubicLengthDerivative<3,CPT>::compute;

};

class AnalyzerTest : public LogFileUnitTest
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
      //ref = 3.5756726e-01;   // old code, with computePrefactor
      ref = 3.5317562e-01;     // new code, with WaveList::implicitInverse
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
      //ref = 8.2574275e+02;  // without the perturbation term
      ref = 1.5803108e+02;    // with the Einstein crystal term
      TEST_ASSERT(checkAverage(filename, ref));

      filename = "out/cubicLengthDerivative.ave";
      //ref = 4.2481444e+03;      // old, before v1.3.4
      //ref = 4.234712800000e+03; // bugfix, without perturbationHamiltonian
      //ref = 1.238404900000e+03; // new in v1.3.4
      ref = 1.500057600000e+03;   // stress term times lambda
      TEST_ASSERT(checkAverage(filename, ref));

      checkHamiltonianAnalyzer();
      checkBinaryStructureFactor();
   }

   /*
   * Set chi of the 1D diblock melt (which has an Einstein crystal
   * perturbation with the default epsilon = chi), and return the
   * Hamiltonian of the current w fields.
   */
   double hamiltonianAtChi(Rp::System<1,CPT>& system, double chi)
   {
      Rp::Simulator<1,CPT>& simulator = system.simulator();
      system.interaction().setChi(0, 1, chi);
      simulator.analyzeChi();
      simulator.perturbation().setup();
      simulator.clearData();
      simulator.computeWc();
      simulator.computeHamiltonian();
      return simulator.hamiltonian();
   }

   /*
   * BinaryChiDerivative with an EinsteinCrystalPerturbation (lambda
   * 0.25, default epsilon): the value must be the derivative of H with
   * respect to chi at fixed w fields, computed here by central finite
   * differences, plus the term G/(2 chi) of the normalization.
   */
   void testChiDerivativePerturbation()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testChiDerivativePerturbation.log");

      // Reference field of the perturbation: half of the w fields
      {
         Rp::System<1,CPT> system;
         initSystem1D(system, "");
         DArray< RField<1,CPT> > w;
         w.allocate(2);
         for (int i = 0; i < 2; ++i) {
            w[i].allocate(system.domain().mesh().dimensions());
            for (int k = 0; k < 16; ++k) {
               w[i][k] = 0.5*system.w().rgrid(i)[k];
            }
         }
         system.domain().fieldIo().writeFieldsRGrid("out/chiDerivRef.rf",
                                  w, system.domain().unitCell(), false);
      }

      Rp::System<1,CPT> system;
      initSystem1D(system, "BdSimulator{ seed 1234"
         " EinsteinCrystalPerturbation{ lambda 0.25"
         " fieldFileName out/chiDerivRef.rf } }");
      system.compute();

      double chi = 10.0;
      double delta = 1.0E-3;
      double hPlus = hamiltonianAtChi(system, chi + delta);
      double hMinus = hamiltonianAtChi(system, chi - delta);
      hamiltonianAtChi(system, chi);
      TestChiDerivative analyzer(system.simulator(), system);
      double value = analyzer.compute();
      double expected = (hPlus - hMinus)/(2.0*delta) + 16.0/(2.0*chi);
      if (verbose() > 0) {
         std::cout << "\nanalyzer = " << value
                   << ", finite differences = " << expected;
      }
      TEST_ASSERT(std::abs(value - expected) < 1.0E-6*std::abs(expected));
   }

   /*
   * Set the cubic cell length of a 3D diblock melt, keeping the values
   * of the w fields on the mesh, and return the Hamiltonian.
   */
   double hamiltonianAtLength(Rp::System<3,CPT>& system, double length)
   {
      FSArray<double, 6> parameters;
      parameters.append(length);
      system.setUnitCell(parameters);
      system.compute();
      Rp::Simulator<3,CPT>& simulator = system.simulator();
      simulator.clearData();
      simulator.computeWc();
      simulator.computeHamiltonian();
      return simulator.hamiltonian();
   }

   /*
   * CubicLengthDerivative of a diblock melt on a 6^3 mesh with an
   * EinsteinCrystalPerturbation (lambda 0.25, reference field half of
   * the w fields): the value must be the central difference of H with
   * respect to the cell length L at fixed w fields on the mesh, plus
   * the term -3 G/(2 L) of the normalization.
   */
   void testLengthDerivativePerturbation()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLengthDerivativePerturbation.log");

      Rp::System<3,CPT> system;
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::string mixture =
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.05 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 6 6 6 lattice cubic } ";
      std::istringstream param(mixture
         + "BdSimulator{ seed 1234"
         " EinsteinCrystalPerturbation{ lambda 0.25"
         " fieldFileName out/lengthDerivRef.rf } } }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      system.setUnitCell(parameters);

      // w fields, and the reference field of the perturbation (half)
      DArray< RField<3,CPT> > w, ref;
      w.allocate(2);
      ref.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         ref[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < 216; ++k) {
            double x = 2.0*M_PI*(k/36)/6.0;
            double y = 2.0*M_PI*((k/6)%6)/6.0;
            w[i][k] = (i == 0) ? std::cos(x) + 0.5*std::sin(y) : 0.0;
            ref[i][k] = 0.5*w[i][k];
         }
      }
      system.domain().fieldIo().writeFieldsRGrid("out/lengthDerivRef.rf",
                               ref, system.domain().unitCell(), false);
      system.w().setRGrid(w);
      system.simulator().analyzeChi();
      system.simulator().perturbation().setup();

      double length = 2.0;
      double delta = 1.0E-4;
      double hPlus = hamiltonianAtLength(system, length + delta);
      double hMinus = hamiltonianAtLength(system, length - delta);
      hamiltonianAtLength(system, length);
      TestLengthDerivative analyzer(system.simulator(), system);
      double value = analyzer.compute();
      double expected = (hPlus - hMinus)/(2.0*delta)
                      - 3.0*216.0/(2.0*length);
      if (verbose() > 0) {
         std::cout << "\nanalyzer = " << value
                   << ", finite differences = " << expected;
      }
      TEST_ASSERT(std::abs(value - expected) < 1.0E-5*std::abs(expected));
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
TEST_ADD(AnalyzerTest, testChiDerivativePerturbation)
TEST_ADD(AnalyzerTest, testLengthDerivativePerturbation)
TEST_END(AnalyzerTest)

#endif
