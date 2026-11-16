#ifndef RPC_ANALYZER_TEST_H
#define RPC_ANALYZER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/analyzer/AnalyzerManager.h>
#include <rp/fts/trajectory/RGridTrajectoryReader.h>
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
         " AnalyzerManager{"
         " HamiltonianAnalyzer{ interval 1"
         " outputFileName out/twice_hamiltonian nSamplePerOutput 1 }"
         " MaxOrderParameter{ interval 1"
         " outputFileName out/twice_maxOrder nSamplePerOutput 1 }"
         " } }");

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
         " AnalyzerManager{ TrajectoryWriter{ interval 1"
         " outputFileName " + filename + " } } }");
      system.simulate(nStep);
   }

   /*
   * Copy a trajectory file, ending the copy after the first nLine
   * lines of frame iFrame (a damaged last frame).
   */
   void writeTruncatedTrajectory(std::string const & inFileName,
                                 std::string const & outFileName,
                                 int iFrame, int nLine)
   {
      std::ifstream in;
      openInputFile(inFileName, in);
      std::ofstream out;
      openOutputFile(outFileName, out);
      std::string label = "i = " + std::to_string(iFrame);
      std::string line;
      int n = -1;
      while (std::getline(in, line)) {
         if (line == label) n = 0;
         if (n >= 0) {
            if (n == nLine) break;
            ++n;
         }
         out << line << "\n";
      }
   }

   /*
   * Copy a trajectory file, replacing the line after each line whose
   * first word is label (e.g. "mesh") by text.
   */
   void writeModifiedTrajectory(std::string const & inFileName,
                                std::string const & outFileName,
                                std::string const & label,
                                std::string const & text)
   {
      std::ifstream in;
      openInputFile(inFileName, in);
      std::ofstream out;
      openOutputFile(outFileName, out);
      std::string line, first;
      bool replace = false;
      while (std::getline(in, line)) {
         out << (replace ? text : line) << "\n";
         std::istringstream words(line);
         first = "";
         words >> first;
         replace = (first == label);
      }
   }

   /*
   * Read frame iFrame of a trajectory file into the w fields of system.
   */
   void readFrame(Rp::System<1,CPT>& system, std::string const & filename,
                  int iFrame)
   {
      Rp::RGridTrajectoryReader<1,CPT> reader(system);
      reader.open(filePrefix() + filename);
      reader.readHeader();
      for (int i = 0; i <= iFrame; ++i) {
         TEST_ASSERT(reader.readFrame());
      }
      reader.close();
   }

   /*
   * Run ANALYZE with an RGridTrajectoryReader. Return false if it
   * throws an Exception.
   */
   bool analyze(Rp::System<1,CPT>& system, int min, int max,
                std::string const & filename)
   {
      try {
         system.simulator().analyze(min, max, "RGridTrajectoryReader",
                                    filePrefix() + filename);
      } catch (Exception& e) {
         return false;
      }
      return true;
   }

   /*
   * ANALYZE 0 1 of a trajectory with frames 0 to 3, by BD and MC
   * simulators that have only an analyzer: the system must be left
   * with frame 1, also if frame 2 is damaged.
   */
   void testAnalyzeLastFrame()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeLastFrame.log");

      writeTrajectory("out/trajectory1D.rf", 3);
      writeTruncatedTrajectory("out/trajectory1D.rf",
                               "out/trajectory1D_damaged.rf", 2, 5);
      std::string simulators[2] = {"BdSimulator", "McSimulator"};
      std::string files[2] = {"out/trajectory1D.rf",
                              "out/trajectory1D_damaged.rf"};
      for (int i = 0; i < 2; ++i) {
         for (int j = 0; j < 2; ++j) {
            Rp::System<1,CPT> system;
            initSystem1D(system, simulators[i]
               + "{ AnalyzerManager{ HamiltonianAnalyzer{ interval 1"
               " outputFileName out/analyzeLastFrame } } }");
            TEST_ASSERT(analyze(system, 0, 1, files[j]));

            // Compare with frame 1
            Rp::System<1,CPT> reference;
            initSystem1D(reference, "");
            readFrame(reference, files[0], 1);
            RFieldComparison<1,CPT> comparison;
            comparison.compare(reference.w().rgrid(), system.w().rgrid());
            TEST_ASSERT(comparison.maxDiff() == 0.0);
         }
      }
   }

   /*
   * ANALYZE 0 0 of a trajectory whose frames are compressed only to an
   * error of 1e-2, by BD and MC simulators with a compressor of
   * epsilon 1e-10 (and a BD step or an MC move): the analyzed frame
   * must be left unchanged.
   */
   void testAnalyzeFirstFrame()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeFirstFrame.log");

      writeTrajectory("out/trajectory1D.rf", 2);
      std::string compressor
            = " LrAmCompressor{ epsilon 1.0e-10 maxItr 200 maxHist 20 }";
      std::string analyzers
            = " AnalyzerManager{ HamiltonianAnalyzer{ interval 1"
              " outputFileName out/analyzeFirstFrame } } }";
      std::string simulators[2];
      simulators[0] = "BdSimulator{ LMBdStep{ mobility 0.01 }"
                      + compressor + analyzers;
      simulators[1] = "McSimulator{ McMoveManager{"
                      " RealMove{ probability 1.0 sigma 0.1 } }"
                      + compressor + analyzers;
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, simulators[i]);
         TEST_ASSERT(analyze(system, 0, 0, "out/trajectory1D.rf"));

         // Compare with frame 0
         Rp::System<1,CPT> reference;
         initSystem1D(reference, "");
         readFrame(reference, "out/trajectory1D.rf", 0);
         RFieldComparison<1,CPT> comparison;
         comparison.compare(reference.w().rgrid(), system.w().rgrid());
         TEST_ASSERT(comparison.maxDiff() == 0.0);
      }
   }

   /*
   * ANALYZE 0 3 of a trajectory with frames 0 to 3, with baseInterval
   * 2 and an analyzer of interval 2: frames 0 and 2 must be sampled.
   */
   void testAnalyzeBaseInterval()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeBaseInterval.log");

      writeTrajectory("out/trajectory1D.rf", 3);
      std::string simulators[2] = {"BdSimulator", "McSimulator"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, simulators[i]
            + "{ AnalyzerManager{ baseInterval 2"
            " HamiltonianAnalyzer{ interval 2"
            " outputFileName out/analyzeBase } } }");
         TEST_ASSERT(analyze(system, 0, 3, "out/trajectory1D.rf"));
         TEST_ASSERT(countLines("out/analyzeBase.dat") == 2);
      }
      Rp::Analyzer<1,CPT>::initStatic();
   }

   /*
   * ANALYZE of a trajectory with frames 0 to 3 must fail if no frame
   * is sampled: with no frame in [min, max] (ANALYZE 10 20), or with
   * an analyzer whose interval (4) divides no index in [min, max]
   * (ANALYZE 1 3).
   */
   void testAnalyzeNoSample()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testAnalyzeNoSample.log");

      writeTrajectory("out/trajectory1D.rf", 3);
      std::string simulators[2] = {"BdSimulator", "McSimulator"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, simulators[i]
            + "{ AnalyzerManager{ HamiltonianAnalyzer{ interval 1"
            " outputFileName out/analyzeNoFrame } } }");
         TEST_ASSERT(!analyze(system, 10, 20, "out/trajectory1D.rf"));

         Rp::System<1,CPT> system2;
         initSystem1D(system2, simulators[i]
            + "{ AnalyzerManager{ HamiltonianAnalyzer{ interval 4"
            " outputFileName out/analyzeNoSample } } }");
         TEST_ASSERT(!analyze(system2, 1, 3, "out/trajectory1D.rf"));
      }
   }

   /*
   * Timers of BD and MC simulators that have only an analyzer, and
   * no compressor (WRITE_TIMERS and CLEAR_TIMERS after ANALYZE).
   */
   void testTimersWithoutCompressor()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testTimersWithoutCompressor.log");

      writeTrajectory("out/trajectory1D.rf", 1);
      std::string simulators[2] = {"BdSimulator", "McSimulator"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, simulators[i]
            + "{ AnalyzerManager{ HamiltonianAnalyzer{ interval 1"
            " outputFileName out/timers } } }");
         TEST_ASSERT(!system.simulator().hasCompressor());
         TEST_ASSERT(analyze(system, 0, 1, "out/trajectory1D.rf"));
         bool success = true;
         try {
            std::ostringstream out;
            system.writeTimers(out);
            system.clearTimers();
         } catch (Exception& e) {
            success = false;
         }
         TEST_ASSERT(success);
      }
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
         initSystem1D(system, "BdSimulator{ AnalyzerManager{"
            " HamiltonianAnalyzer{ interval 1"
            " outputFileName out/analyzePrefix } } }");
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

   /*
   * SIMULATE 4 BD steps of the 1D diblock melt with a ramp of the cell
   * parameter from 1.5 to 1.9, with the given analyzer blocks.
   */
   void simulateCellRamp(Rp::System<1,CPT>& system,
                         std::string const & analyzers)
   {
      initSystem1D(system,
         "BdSimulator{ seed 1234 LMBdStep{ mobility 0.01 }"
         " LrAmCompressor{ epsilon 1.0e-8 maxItr 100 maxHist 20 }"
         " LinearRamp{ nParameter 1 parameters[ cell_param 0 0.4 ] }"
         " AnalyzerManager{ " + analyzers + " } }");
      system.simulate(4);
   }

   /*
   * SIMULATE with a ramp of the cell parameter (1.5 to 1.9 in 4 BD
   * steps), with a TrajectoryWriter and a ConcentrationWriter. Frame k
   * of both files, read by ANALYZE or by an RGridTrajectoryReader, must
   * set the unit cell of step k.
   */
   void testTrajectoryCellRamp()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testTrajectoryCellRamp.log");

      {
         Rp::System<1,CPT> system;
         simulateCellRamp(system,
            "TrajectoryWriter{ interval 1 outputFileName out/cellRamp_w.rf }"
            " ConcentrationWriter{ interval 1"
            " outputFileName out/cellRamp_c.rf }");
      }
      for (int k = 0; k <= 4; ++k) {
         double cell = 1.5 + 0.1*k;

         Rp::System<1,CPT> system;
         initSystem1D(system, "BdSimulator{ AnalyzerManager{"
            " HamiltonianAnalyzer{ interval 1"
            " outputFileName out/analyzeCellRamp } } }");
         TEST_ASSERT(analyze(system, 0, k, "out/cellRamp_w.rf"));
         TEST_ASSERT(std::abs(system.domain().unitCell().parameter(0)
                              - cell) < 1.0E-10);

         Rp::System<1,CPT> system2;
         initSystem1D(system2, "");
         readFrame(system2, "out/cellRamp_c.rf", k);
         TEST_ASSERT(std::abs(system2.domain().unitCell().parameter(0)
                              - cell) < 1.0E-10);
      }
   }

   /*
   * ANALYZE of a trajectory whose header gives another number of
   * monomer types, or whose frames give other mesh dimensions, than
   * the system must fail.
   */
   void testTrajectoryMismatch()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testTrajectoryMismatch.log");

      writeTrajectory("out/trajectory1D.rf", 1);
      writeModifiedTrajectory("out/trajectory1D.rf",
                              "out/trajectory1D_nMonomer.rf",
                              "N_monomer", "          3");
      writeModifiedTrajectory("out/trajectory1D.rf",
                              "out/trajectory1D_mesh.rf",
                              "mesh", "          32");
      std::string files[2] = {"out/trajectory1D_nMonomer.rf",
                              "out/trajectory1D_mesh.rf"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, "BdSimulator{ AnalyzerManager{"
            " HamiltonianAnalyzer{ interval 1"
            " outputFileName out/analyzeMismatch } } }");
         TEST_ASSERT(!analyze(system, 0, 1, files[i]));
      }
   }

   /*
   * BinaryStructureFactor with a ramp of the cell parameter: SIMULATE
   * must stop with an error. ANALYZE of a trajectory written
   * under such a ramp must use the unit cell of the first sampled frame
   * (frame 2, for ANALYZE 1 2 with interval 2), and must stop if a later
   * sampled frame has another unit cell (ANALYZE 0 4).
   */
   void testStructureFactorCellRamp()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testStructureFactorCellRamp.log");

      {
         Rp::System<1,CPT> system;
         bool success = true;
         try {
            simulateCellRamp(system, "BinaryStructureFactor{ interval 1"
                                     " outputFileName out/cellRampS }");
         } catch (Exception& e) {
            success = false;
         }
         TEST_ASSERT(!success);
      }

      {
         Rp::System<1,CPT> system;
         simulateCellRamp(system, "TrajectoryWriter{ interval 1"
                                  " outputFileName out/cellRamp_w.rf }");
      }
      int min[2] = {1, 0};
      int max[2] = {2, 4};
      std::string interval[2] = {"2", "1"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system, "BdSimulator{ AnalyzerManager{"
            " baseInterval " + interval[i] +
            " BinaryStructureFactor{ interval " + interval[i] +
            " outputFileName out/cellRampS } } }");
         bool success = analyze(system, min[i], max[i],
                                "out/cellRamp_w.rf");
         TEST_ASSERT(success == (i == 0));
      }
      Rp::Analyzer<1,CPT>::initStatic();

      // Wavenumber of the first nonzero bunch, for the cell of frame 2
      std::ifstream file;
      openInputFile("out/cellRampS", file);
      double q, s;
      file >> q >> s >> q >> s;
      TEST_ASSERT(std::abs(q - 2.0*M_PI/1.7) < 1.0E-6);
   }

};

TEST_BEGIN(AnalyzerTest)
TEST_ADD(AnalyzerTest, testAnalyzeTrajectory)
TEST_ADD(AnalyzerTest, testSimulateTwice)
TEST_ADD(AnalyzerTest, testAnalyzeLastFrame)
TEST_ADD(AnalyzerTest, testAnalyzeFirstFrame)
TEST_ADD(AnalyzerTest, testAnalyzeBaseInterval)
TEST_ADD(AnalyzerTest, testAnalyzeNoSample)
TEST_ADD(AnalyzerTest, testTimersWithoutCompressor)
TEST_ADD(AnalyzerTest, testAnalyzeInputPrefix)
TEST_ADD(AnalyzerTest, testTrajectoryCellRamp)
TEST_ADD(AnalyzerTest, testTrajectoryMismatch)
TEST_ADD(AnalyzerTest, testStructureFactorCellRamp)
TEST_END(AnalyzerTest)

#endif
