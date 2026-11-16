#ifndef RPC_MC_SIMULATOR_TEST_H
#define RPC_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/montecarlo/McMove.h>
#include <rp/fts/montecarlo/McMoveManager.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <rp/fts/compressor/Compressor.h>
#include <pscf/field/cpp/RFieldComparison.h>

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

   /*
   * Read a 1D diblock melt (16 grid points) with the given McSimulator
   * block, and set w fields with a cosine exchange component of the
   * given amplitude.
   */
   void initSystem1D(Rp::System<1,CPT>& system,
                     std::string const & simulatorBlock,
                     double amplitude = 1.0)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.02 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar } "
         + simulatorBlock + " }\n");
      system.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(1.5);
      system.setUnitCell(parameters);
      DArray< RField<1,CPT> > w;
      allocateRGridFields(system, w);
      for (int i = 0; i < 2; ++i) {
         for (int k = 0; k < 16; ++k) {
            w[i][k] = (i == 0) ? amplitude*std::cos(2.0*M_PI*k/16.0) : 0.0;
         }
      }
      system.w().setRGrid(w);
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

   /*
   * A run in which every MC move is rejected (RealMove with sigma 1),
   * so also the last one: the system must be left with the c fields of
   * its final w fields.
   */
   void testRejectedLastMove()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testRejectedLastMove.log");

      Rp::System<1,CPT> system;
      initSystem1D(system,
          "McSimulator{ seed 1234 McMoveManager{\n"
          " RealMove{ probability 1.0 sigma 1.0 }\n"
          " }\n"
          " LrAmCompressor{ epsilon 1.0e-8 maxItr 100 maxHist 20 } }");
      Rp::McSimulator<1,CPT>& simulator
          = dynamic_cast< Rp::McSimulator<1,CPT>& >(system.simulator());
      system.simulate(3);
      Rp::McMove<1,CPT>& move = simulator.mcMoveManager()[0];
      TEST_ASSERT(move.nAttempt() == 3);
      TEST_ASSERT(move.nAccept() == 0);
      TEST_ASSERT(hasCurrentCFields(system));
   }

};

TEST_BEGIN(McSimulatorTest)
TEST_ADD(McSimulatorTest, testMcSimulateDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateTriblocks)
TEST_ADD(McSimulatorTest, testMcSimulateBdMoveDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateShiftDiblocks)
TEST_ADD(McSimulatorTest, testBdMoveFailedLastStep)
TEST_ADD(McSimulatorTest, testRejectedLastMove)
TEST_END(McSimulatorTest)

#endif
