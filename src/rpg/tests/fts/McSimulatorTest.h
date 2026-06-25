#ifndef RPG_MC_SIMULATOR_TEST_H
#define RPG_MC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/field/FieldIo.h>

#include <rp/fts/compressor/Compressor.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class McSimulatorTest : public LogFileUnitTest
{

public:

   void setUp()
   {  setVerbose(0); }
   
   template <int D>
   void initSystem(Rp::System<D, Rpg::Types<D> >& system, std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();

   }
   
   template <int D>
   void initSimulator(Rp::McSimulator<D, Rpg::Types<D> >& simulator, std::string filename)
   {
      std::ifstream in;
      openInputFile(filename, in);
      simulator.readParam(in);
      in.close();
   }

   /*
   * Allocate an array of rgrid fields.
   */
   template <int D>
   void allocateRGridFields(Rp::System<D, Rpg::Types<D> > const & system,
                            DArray< RField<D> >& fields)
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
   void readRGridFields(Rp::System<D, Rpg::Types<D> > const & system,
                        std::string filename,
                        DArray< RField<D> >& fields,
                        UnitCell<D>& unitCell)
   {
      allocateRGridFields(system, fields);
      Rp::FieldIo<D, Types<D> > const & fieldIo = system.domain().fieldIo();
      fieldIo.readFieldsRGrid(filename, fields, unitCell);
   }
   
   void testMcSimulateDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateDiblocks.log");
      
      Rp::System<3, Rpg::Types<3> > system;
      initSystem(system, "in/param_system_disordered");
      
      Rp::McSimulator<3, Rpg::Types<3> > simulator(system);
      initSimulator(simulator, "in/param_McSimulator");
      
      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock.rf");

      // Read reference field
      DArray< RField<3> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

   void testMcSimulateBdMoveDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateBdMoveDiblocks.log");
      
      Rp::System<3, Rpg::Types<3> > system;
      initSystem(system, "in/param_system_disordered");
      
      Rp::McSimulator<3, Rpg::Types<3> > simulator(system);
      initSimulator(simulator, "in/param_McSimulator_BdMove");
      
      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock_bdMove.rf");

      // Read reference field
      DArray< RField<3> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_bdMove_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

   void testMcSimulateShiftDiblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateShiftDiblocks.log");
      
      Rp::System<3, Rpg::Types<3> > system;
      initSystem(system, "in/param_system_disordered");
      
      Rp::McSimulator<3, Rpg::Types<3> > simulator(system);
      initSimulator(simulator, "in/param_McSimulator_ShiftMove");
      
      system.w().readRGrid("in/w_dis.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_diblock_shift.rf");

      // Read reference field
      DArray< RField<3> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_diblock_shift_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }
   
   void testMcSimulateTriblocks()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testMcSimulateTriblocks.log");
      
      Rp::System<3, Rpg::Types<3> > system;
      initSystem(system, "in/param_system_triblock");
      
      Rp::McSimulator<3, Rpg::Types<3> > simulator(system);
      initSimulator(simulator, "in/param_triblock_McSimulator");
      
      system.w().readRGrid("in/w_triblock.rf");
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid("out/w_mc_triblock.rf");
      
      // Read reference field
      DArray< RField<3> > rf_0;
      UnitCell<3> unitCell;
      readRGridFields(system,"in/w_mc_triblock_ref.rf", rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);

   }

};

TEST_BEGIN(McSimulatorTest)
TEST_ADD(McSimulatorTest, testMcSimulateDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateBdMoveDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateShiftDiblocks)
TEST_ADD(McSimulatorTest, testMcSimulateTriblocks)
TEST_END(McSimulatorTest)

#endif
