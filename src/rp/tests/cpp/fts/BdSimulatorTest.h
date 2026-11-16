#ifndef RPC_BD_SIMULATOR_TEST_H
#define RPC_BD_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/ExplicitBdStep.h>
#include <rp/fts/brownian/PredCorrBdStep.h>
#include <rp/fts/brownian/LMBdStep.h>
#include <rp/system/System.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <rp/fts/compressor/Compressor.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/interaction/Interaction.h>

#include <util/misc/FileMaster.h>
#include <util/containers/FSArray.h>
#include <util/tests/LogFileUnitTest.h>

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace Util;
using namespace Pscf;

class BdSimulatorTest : public LogFileUnitTest
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
   * Generic BdSimulator test function template.
   */ 
   void testBdSimulator(Rp::System<3,CPT>& system, 
                        std::string systemfilename,
                        std::string simulatorfilename,
                        std::string infieldsfilename,
                        std::string outfieldsfilename,
                        std::string reffieldsfilename, 
                        char const * outfilename)
   {
      openLogFile(outfilename);
      initSystem(system, systemfilename);

      Rp::BdSimulator<3,CPT> simulator(system);
      initSimulator(simulator, simulatorfilename);

      system.w().readRGrid(infieldsfilename);
      simulator.compressor().compress();
      simulator.simulate(50);
      system.w().writeRGrid(outfieldsfilename);

      // Read reference field
      DArray< RField<3,CPT> > rf_0;
      UnitCell<3> unitCell;

      readRGridFields(system,reffieldsfilename, rf_0, unitCell);

      // Compare with reference fields
      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, system.w().rgrid());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-7);
   }

   
   void testLMBdStepDiblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_disordered",
                     "in/param_LMBdStep",
                     "in/w_dis.rf",
                     "out/w_LMBd_diblock.rf",
                     "in/w_LMBd_diblock_ref.rf", 
                     "out/testLMBdStepDiblocks.log");                                                               
   }
   
   void testLMBdStepTriblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_triblock", 
                     "in/param_LMBdStep",
                     "in/w_triblock.rf", 
                     "out/w_LMBd_triblock.rf",
                     "in/w_LMBd_triblock_ref.rf", 
                     "out/testLMBdStepTriblocks.log");  
   }

   void testExplicitBdStepDiblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_disordered",
                     "in/param_explicitBdStep",
                     "in/w_dis.rf",
                     "out/w_explicitBd_diblock.rf",
                     "in/w_explicitBd_diblock_ref.rf", 
                     "out/testExplicitBdStepDiblocks.log");                                                               
   }
   
   void testExplicitBdStepTriblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_triblock", 
                     "in/param_explicitBdStep",
                     "in/w_triblock.rf", 
                     "out/w_explicitBd_triblock.rf",
                     "in/w_explicitBd_triblock_ref.rf", 
                     "out/testExplicitBdStepTriblocks.log");  
   }

   void testPredCorrBdStepDiblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_disordered",
                     "in/param_predCorrBdStep",
                     "in/w_dis.rf",
                     "out/w_predCorrBdStep_diblock.rf",
                     "in/w_predCorrBdStep_diblock_ref.rf", 
                     "out/testPredCorrBdStepDiblocks.log");                                                               
   }
   
   void testPredCorrBdStepTriblocks()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testBdSimulator(system, "in/param_system_triblock", 
                     "in/param_predCorrBdStep",
                     "in/w_triblock.rf", 
                     "out/w_predCorrBdStep_triblock.rf",
                     "in/w_predCorrBdStep_triblock_ref.rf", 
                     "out/testPredCorrBdStepTriblocks.log");  
   }

   /*
   * Ramp of vMonomer, with an LrCompressor: every step must converge.
   */
   void testVMonomerRampLrCompressor()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testVMonomerRampLrCompressor.log");

      // Diblock melt, vMonomer ramp from 0.01 to 0.002
      Rp::System<1,CPT> system;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.02 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar }"
         " BdSimulator{ seed 1234 LMBdStep{ mobility 0.001 }"
         " LrCompressor{ epsilon 1.0e-8 maxItr 100 }"
         " LinearRamp{ nParameter 1 parameters v_monomer -0.008 } } }\n");
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

      system.simulate(20);
      TEST_ASSERT(system.simulator().iStep() == 20);
      TEST_ASSERT(std::abs(system.mixture().vMonomer() - 0.002) < 1.0E-12);
   }

   /*
   * Read a 1D diblock melt (16 grid points) with the given BdSimulator
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
      w.allocate(2);
      for (int i = 0; i < 2; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < 16; ++k) {
            w[i][k] = (i == 0) ? amplitude*std::cos(2.0*M_PI*k/16.0) : 0.0;
         }
      }
      system.w().setRGrid(w);
   }

   /*
   * Check the state stored by the simulator after a run: the w fields
   * must be compressed for the current parameters, and the stored
   * forces dc and Hamiltonian must be those of the current fields and
   * parameters.
   */
   void checkFinalState(Rp::System<1,CPT>& system)
   {
      Rp::Simulator<1,CPT>& simulator = system.simulator();
      TEST_ASSERT(simulator.hasDc());
      TEST_ASSERT(simulator.hasHamiltonian());
      DArray< RField<1,CPT> > dc;
      dc.allocate(1);
      dc[0].allocate(system.domain().mesh().dimensions());
      dc[0] = simulator.dc(0);
      double hamiltonian = simulator.hamiltonian();

      // Incompressibility of the c fields of the final w fields
      system.compute();
      double maxResid = 0.0;
      for (int k = 0; k < 16; ++k) {
         double resid = system.c().rgrid(0)[k] + system.c().rgrid(1)[k];
         resid = std::abs(resid - 1.0);
         if (resid > maxResid) maxResid = resid;
      }

      // Recompute dc and the Hamiltonian
      simulator.clearData();
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();
      simulator.computeHamiltonian();
      RFieldComparison<1,CPT> comparison;
      comparison.compare(dc, simulator.dc());
      double hDiff = std::abs(simulator.hamiltonian() - hamiltonian);
      if (verbose() > 0) {
         std::cout << "\nmax |sum c - 1| = " << maxResid
                   << ", dc diff = " << comparison.maxDiff()
                   << ", H diff = " << hDiff;
      }
      TEST_ASSERT(maxResid < 1.0E-6);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-8);
      TEST_ASSERT(hDiff < 1.0E-10*std::abs(hamiltonian));
   }

   /*
   * Ramps of chi (no compression needed after a ramp step) and of a
   * statistical segment length (the c fields of given w fields change):
   * after the run, the state must be that of the final parameters.
   */
   void testRampState()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testRampState.log");

      std::string ramps[2] = {"chi 0 1 2.0", "kuhn 0 0.2"};
      for (int i = 0; i < 2; ++i) {
         Rp::System<1,CPT> system;
         initSystem1D(system,
             "BdSimulator{ seed 1234 LMBdStep{ mobility 0.001 }"
             " LrAmCompressor{ epsilon 1.0e-8 maxItr 100 maxHist 20 }"
             " LinearRamp{ nParameter 1 parameters[ " + ramps[i]
             + " ] } }");
         system.simulate(4);
         TEST_ASSERT(system.simulator().iStep() == 4);
         checkFinalState(system);
      }
   }

   /*
   * A chi ramp (10 to 12) in a run in which every BD step fails (the
   * compressor has maxItr = 1, the run starts from homogeneous fields).
   * The ramp must still end at the final value after the last step.
   */
   void testRampFailedSteps()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testRampFailedSteps.log");

      Rp::System<1,CPT> system;
      initSystem1D(system,
          "BdSimulator{ seed 1234 LMBdStep{ mobility 0.05 }"
          " LrAmCompressor{ epsilon 1.0e-8 maxItr 1 maxHist 20 }"
          " LinearRamp{ nParameter 1 parameters[ chi 0 1 2.0 ] } }",
          0.0);
      system.simulate(3);
      TEST_ASSERT(system.simulator().iTotalStep() == 3);
      TEST_ASSERT(system.simulator().iStep() == 0);
      double chi = system.interaction().chi(0,1);
      TEST_ASSERT(std::abs(chi - 12.0) < 1.0E-12);
   }

};

TEST_BEGIN(BdSimulatorTest)
TEST_ADD(BdSimulatorTest, testLMBdStepDiblocks)
TEST_ADD(BdSimulatorTest, testLMBdStepTriblocks)
TEST_ADD(BdSimulatorTest, testExplicitBdStepDiblocks)
TEST_ADD(BdSimulatorTest, testExplicitBdStepTriblocks)
TEST_ADD(BdSimulatorTest, testPredCorrBdStepDiblocks)
TEST_ADD(BdSimulatorTest, testPredCorrBdStepTriblocks)
TEST_ADD(BdSimulatorTest, testVMonomerRampLrCompressor)
TEST_ADD(BdSimulatorTest, testRampState)
TEST_ADD(BdSimulatorTest, testRampFailedSteps)
TEST_END(BdSimulatorTest)

#endif
