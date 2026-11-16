#ifndef RPC_MC_SIMULATOR_TEST_H
#define RPC_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <rp/fts/compressor/Compressor.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/backend/cpp/VecOp.h>
#include <pscf/backend/cpp/Reduce.h>

#include <util/misc/FileMaster.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

class McSimulatorTest : public LogFileUnitTest
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
   void initSimulator(Rp::McSimulator<D,CPT>& simulator, std::string filename)
   {
      Rp::Analyzer<D,CPT>::initStatic();
      int baseInterval = Rp::Analyzer<D,CPT>::baseInterval;
      TEST_ASSERT(baseInterval == 1);

      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * Allocate an array of rgrid fields.
   */
   template <int D>
   void allocateRGridFields(Rp::System<D,CPT> const & system,
                            DArray< RField<D,CPT> >& fields)
   {
      // Check and allocate outer DArray
      int nMonomer = system.mixture().nMonomer();
      UTIL_CHECK(nMonomer > 0);
      if (!fields.isAllocated()) {
         fields.allocate(nMonomer);
      }
      UTIL_CHECK(fields.capacity() == nMonomer);

      // Allocate fields
      Mesh<D> const & mesh = system.domain().mesh();
      IntVec<D> const & meshDimensions = mesh.dimensions();
      int meshSize = mesh.size();
      UTIL_CHECK(meshSize > 0);
      for (int i = 0; i < nMonomer; ++i) {
         if (!fields[i].isAllocated()) {
            fields[i].allocate(meshDimensions);
         }
         UTIL_CHECK(fields[i].capacity() == meshSize);
      }
   }

   /*
   * Read r-grid fields into an array.
   */
   template <int D>
   void readRGridFields(Rp::System<D,CPT> const & system,
                        std::string filename,
                        DArray< RField<D,CPT> >& fields,
                        UnitCell<D>& unitCell)
   {
      allocateRGridFields(system, fields);
      Rp::FieldIo<D,CPT> const & fieldIo = system.domain().fieldIo();
      fieldIo.readFieldsRGrid(filename, fields, unitCell);
   }

   void testMcSimulateDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateDiblocks.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_disordered");

      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/param_McSimulator");

      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock.rf");

      // Read reference field
      DArray< RField<3,CPT> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

   void testMcSimulateBdMoveDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateBdMoveDiblocks.log");
      
      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_disordered");
      
      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/param_McSimulator_BdMove");
      
      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock_bdMove.rf");

      // Read reference field
      DArray< RField<3,CPT> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_bdMove_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

   /*
   * Noise of a BdMove with a single LM step (nStep = 1).
   *
   * The step adds eta - mobility*dc to the exchange component wc(0),
   * where eta, the sum of the "old" and "new" noise, must have the
   * Euler-Maruyama variance 2*mobility*meshSize/vSystem at each grid
   * point. Without the larger variance of the first old noise, eta
   * would have half of this variance.
   */
   void testBdMoveFirstStepNoise()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testBdMoveFirstStepNoise.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_disordered");
      Rp::McSimulator<3,CPT> simulator(system);
      std::istringstream in(
         "McSimulator{ seed 1784106389"
         " McMoveManager{\n"
         " BdMove{ probability 1.0 mobility 1.0e-4"
         " nStep 1 }\n"
         " }\n"
         " LrAmCompressor{ epsilon 1.0e-6 maxItr 200 maxHist 30 } }\n");
      simulator.readParam(in);
      const double mobility = 1.0E-4;

      // Compressed initial state, its wc(0) and force dc(0)
      system.w().readRGrid("in/w_dis.rf");
      simulator.analyzeChi();
      TEST_ASSERT(simulator.compressor().compress() == 0);
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();

      // eta = wc(0) after - wc(0) before + mobility*dc(0) before
      RField<3,CPT> eta;
      eta.allocate(system.domain().mesh().dimensions());
      VecOp::eqV(eta, simulator.wc(0));
      VecOp::mulEqS(eta, -1.0);
      VecOp::addEqVc(eta, simulator.dc(0), mobility);
      simulator.simulate(1);
      simulator.computeWc();
      VecOp::addEqV(eta, simulator.wc(0));

      // Variance of eta over the 16^3 grid points (relative standard
      // error of the estimate about 0.022)
      const int meshSize = system.domain().mesh().size();
      const double vSystem = system.domain().unitCell().volume();
      double mean = Reduce::sum(eta)/double(meshSize);
      double variance = Reduce::sumSq(eta)/double(meshSize) - mean*mean;
      double expected = 2.0*mobility*double(meshSize)/vSystem;
      Log::file() << "Variance of eta / expected value = "
                  << variance/expected << std::endl;
      TEST_ASSERT(std::abs(variance/expected - 1.0) < 0.1);
   }

   void testMcSimulateShiftDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateShiftDiblocks.log");
      
      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_disordered");
      
      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/param_McSimulator_ShiftMove");
      
      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock_shift.rf");

      // Read reference field
      DArray< RField<3,CPT> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_shift_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

   void testMcSimulateTriblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateTriblocks.log");

      Rp::System<3,CPT> system;
      initSystem(system, "in/param_system_triblock");

      Rp::McSimulator<3,CPT> simulator(system);
      initSimulator(simulator, "in/param_triblock_McSimulator");

      system.w().readRGrid("in/w_triblock.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_triblock.rf");

      // Read reference field
      DArray< RField<3,CPT> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_triblock_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);
   }


   /*
   * Are the c fields of the system those of its current w fields?
   */
   template <int D>
   bool hasCurrentCFields(Rp::System<D,CPT>& system)
   {
      if (!system.c().hasData()) return false;
      DArray< RField<D,CPT> > c;
      allocateRGridFields(system, c);
      for (int i = 0; i < c.capacity(); ++i) {
         c[i] = system.c().rgrid()[i];
      }
      system.compute();
      RFieldComparison<D,CPT> comparison;
      comparison.compare(c, system.c().rgrid());
      return (comparison.maxDiff() < 1.0E-10);
   }

   /*
   * A BdMove whose compressor fails in every BD step (maxItr = 1), so
   * also in the last one (issue #274). The McSimulator is that of the
   * System, as in pscf_rpc. The move must complete and leave the c
   * fields of the final w fields.
   */
   void testBdMoveFailedLastStep()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testBdMoveFailedLastStep.log");

      Rp::System<1,CPT> system;
      initSystem(system, "in/param_system_1D_BdMove_fail");
      TEST_ASSERT(system.hasSimulator());
      FSArray<double, 6> parameters;
      parameters.append(2.0);
      system.setUnitCell(parameters);
      DArray< RField<1,CPT> > w;
      allocateRGridFields(system, w);
      for (int i = 0; i < w.capacity(); ++i) {
         for (int j = 0; j < w[i].capacity(); ++j) {
            w[i][j] = 0.0;
         }
      }
      system.w().setRGrid(w);

      bool success = true;
      try {
         system.simulator().simulate(1);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      TEST_ASSERT(hasCurrentCFields(system));
   }

};

TEST_BEGIN(McSimulatorTest)
TEST_ADD(McSimulatorTest, testMcSimulateDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateTriblocks)
TEST_ADD(McSimulatorTest, testMcSimulateBdMoveDiblocks)
TEST_ADD(McSimulatorTest, testBdMoveFirstStepNoise)
TEST_ADD(McSimulatorTest, testMcSimulateShiftDiblocks)
TEST_ADD(McSimulatorTest, testBdMoveFailedLastStep)
TEST_END(McSimulatorTest)

#endif
