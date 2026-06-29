#ifndef RPG_EINSTEIN_CRYSTAL_PERTURBATION_TEST_H
#define RPG_EINSTEIN_CRYSTAL_PERTURBATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/fts/brownian/BdSimulator.h>

#include <prdc/field/cuda/RField.h>
#include <prdc/field/cuda/RFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;
using namespace Pscf::Rpg;

class PerturbationTest : public LogFileUnitTest
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
   void initSimulator(Rp::BdSimulator<D, Rpg::Types<D> >& simulator, std::string filename)
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

   /*
   * Generic BdSimulator test function template.
   */ 
   void testPerturbation(Rp::System<3, Rpg::Types<3> >& system, 
                        std::string systemfilename,
                        std::string simulatorfilename,
                        std::string infieldsfilename,
                        std::string outfieldsfilename,
                        std::string reffieldsfilename, 
                        char const * outfilename)
   {
      openLogFile(outfilename);
      initSystem(system, systemfilename);

      Rp::BdSimulator<3, Rpg::Types<3> > simulator(system);
      initSimulator(simulator, simulatorfilename);

      system.w().readRGrid(infieldsfilename);
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid(outfieldsfilename);

      // Read reference field
      DArray< RField<3> > rf_0;
      UnitCell<3> unitCell;

      readRGridFields(system,reffieldsfilename, rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);
   }
   
   void testEinsteinCrystalPerturbation()
   {
      printMethod(TEST_FUNC);
      Rp::System<3, Rpg::Types<3> > system;
      testPerturbation(system, "in/param_system_disordered",
                     "in/param_LMBdStep_EinsteinCrystalPerturbation",
                     "in/w_dis.rf",
                     "out/w_LMBdStep_EinsteinCrystalPerturbation.rf",
                     "in/w_LMBdStep_EinsteinCrystalPerturbation_ref.rf", 
                     "out/testEinsteinCrystalPerturbation.log");                                                               
   }

};

TEST_BEGIN(PerturbationTest)
TEST_ADD(PerturbationTest, testEinsteinCrystalPerturbation)
TEST_END(PerturbationTest)

#endif
