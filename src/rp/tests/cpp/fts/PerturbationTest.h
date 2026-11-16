#ifndef RPC_EINSTEIN_CRYSTAL_PERTURBATION_TEST_H
#define RPC_EINSTEIN_CRYSTAL_PERTURBATION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/perturbation/Perturbation.h>
#include <rp/field/WFields.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldComparison.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/containers/FSArray.h>

#include <fstream>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class PerturbationTest : public LogFileUnitTest
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
   void testPerturbation(Rp::System<3,CPT>& system, 
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
   
   void testEinsteinCrystalPerturbation()
   {
      printMethod(TEST_FUNC);
      Rp::System<3,CPT> system;
      testPerturbation(system, "in/param_system_disordered",
                     "in/param_LMBdStep_EinsteinCrystalPerturbation",
                     "in/w_dis.rf",
                     "out/w_LMBdStep_EinsteinCrystalPerturbation.rf",
                     "in/w_LMBdStep_EinsteinCrystalPerturbation_ref.rf", 
                     "out/testEinsteinCrystalPerturbation.log");                                                               
   }

   /*
   * Read a 1D system with 16 grid points from a parameter string, and
   * set w fields w_i(x) = cos(2 pi x/L + i).
   */
   void initSystem1D(Rp::System<1,CPT>& system, std::string const & param)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());
      std::istringstream in(param);
      system.readParam(in);
      FSArray<double, 6> parameters;
      parameters.append(1.5);
      system.setUnitCell(parameters);
      int nMonomer = system.mixture().nMonomer();
      DArray< RField<1,CPT> > w;
      w.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(system.domain().mesh().dimensions());
         for (int k = 0; k < 16; ++k) {
            w[i][k] = std::cos(2.0*M_PI*k/16.0 + i);
         }
      }
      system.w().setRGrid(w);
   }

   /*
   * Recompute and return the Hamiltonian of the final state of a run.
   */
   double finalHamiltonian(Rp::System<1,CPT>& system)
   {
      Rp::Simulator<1,CPT>& simulator = system.simulator();
      simulator.clearData();
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();
      simulator.computeHamiltonian();
      return simulator.hamiltonian();
   }

   /*
   * Diblock melt with an EinsteinCrystalPerturbation with lambda = 0
   * (so that H is the Einstein crystal Hamiltonian) and the default
   * epsilon = chi, under a ramp of chi from 10 to 12. After the run, H
   * must be that of a perturbation set up for the final chi.
   */
   void testEinsteinCrystalChiRamp()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testEinsteinCrystalChiRamp.log");

      std::string mixture =
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.05 }"
         " Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 16 lattice lamellar } ";

      // Reference field: 3 times the initial w fields
      {
         Rp::System<1,CPT> system;
         initSystem1D(system, mixture + "}");
         DArray< RField<1,CPT> > w;
         w.allocate(2);
         for (int i = 0; i < 2; ++i) {
            w[i].allocate(system.domain().mesh().dimensions());
            for (int k = 0; k < 16; ++k) {
               w[i][k] = 3.0*system.w().rgrid(i)[k];
            }
         }
         system.w().setRGrid(w);
         system.w().writeRGrid("out/ecChiRampRef.rf");
      }

      Rp::System<1,CPT> system;
      initSystem1D(system, mixture
         + "BdSimulator{ seed 1234 LMBdStep{ mobility 0.001 }"
         " LrAmCompressor{ epsilon 1.0e-8 maxItr 100 maxHist 20 }"
         " EinsteinCrystalPerturbation{ lambda 0.0"
         " fieldFileName out/ecChiRampRef.rf }"
         " LinearRamp{ nParameter 1 parameters[ chi 0 1 2.0 ] } } }");
      system.simulate(4);

      double hamiltonian = finalHamiltonian(system);
      system.simulator().perturbation().setup();
      system.simulator().computeHamiltonian();
      double hamiltonianSetup = system.simulator().hamiltonian();
      TEST_ASSERT(std::abs(hamiltonianSetup) > 1.0);
      TEST_ASSERT(std::abs(hamiltonian - hamiltonianSetup)
                  < 1.0E-10*std::abs(hamiltonianSetup));
   }

   /*
   * ABC triblock melt with an EinsteinCrystalPerturbation with lambda
   * = 0 and given epsilon values, whose reference field is the initial
   * w field, under a ramp of chi_BC from 12 to 22, which rotates the
   * eigenvectors of the projected chi matrix. With a mobility of
   * 1.0E-14, the composition components of the fields stay at those of
   * the reference field, so H (the Einstein crystal Hamiltonian) must
   * stay close to zero.
   */
   void testEinsteinCrystalChiRampTernary()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testEinsteinCrystalChiRampTernary.log");

      Rp::System<1,CPT> system;
      initSystem1D(system,
         "System{ Mixture{ nMonomer 3 monomers[ 1.0 1.0 1.0 ]"
         " nPolymer 1"
         " Polymer{ type linear nBlock 3 blocks[ 0 0.25 1 0.5 2 0.25 ]"
         " phi 1.0 } vMonomer 0.01 ds 0.05 }"
         " Interaction{ chi( 1 0 10.0 2 0 20.0 2 1 12.0 ) }"
         " Domain{ mesh 16 lattice lamellar }"
         " BdSimulator{ seed 1234 LMBdStep{ mobility 1.0e-14 }"
         " LrAmCompressor{ epsilon 1.0e-10 maxItr 300 maxHist 20 }"
         " EinsteinCrystalPerturbation{ lambda 0.0 epsilon[ 8.0 6.0 ]"
         " fieldFileName out/ecChiRampTernaryRef.rf }"
         " LinearRamp{ nParameter 1 parameters[ chi 2 1 10.0 ] } } }");
      system.w().writeRGrid("out/ecChiRampTernaryRef.rf");
      system.simulate(5);

      double hamiltonian = finalHamiltonian(system);
      TEST_ASSERT(std::abs(hamiltonian) < 1.0E-6);
   }

};

TEST_BEGIN(PerturbationTest)
TEST_ADD(PerturbationTest, testEinsteinCrystalPerturbation)
TEST_ADD(PerturbationTest, testEinsteinCrystalChiRamp)
TEST_ADD(PerturbationTest, testEinsteinCrystalChiRampTernary)
TEST_END(PerturbationTest)

#endif
