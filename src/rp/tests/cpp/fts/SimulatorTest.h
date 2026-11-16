#ifndef RPC_SIMULATOR_TEST_H
#define RPC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#include <rp/fts/brownian/BdSimulator.h>
#include <rp/fts/ramp/LinearRamp.h>

#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/interaction/Interaction.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <fstream>

using namespace Util;
using namespace Pscf;

class SimulatorTest : public LogFileUnitTest
{

   Rp::System<3,CPT> system;

public:


   SimulatorTest()
    : system()
   {}

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {  setVerbose(0); }

   void initSystem(std::string filename)
   {
      system.fileMaster().setInputPrefix(filePrefix());
      system.fileMaster().setOutputPrefix(filePrefix());

      //openLogFile("out/testSystem.log");
      //ParamComponent::setEcho(true);

      std::ifstream in;
      openInputFile(filename, in);
      system.readParam(in);
      in.close();
   }

   void testAnalyzeChi()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param_Simulator");
      Rp::Simulator<3,CPT> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();
      double const eps = 1.0E-8;

      double chi = system.interaction().chi(0,1);
      TEST_ASSERT(std::abs(system.interaction().chi(0,0)) < eps);
      TEST_ASSERT(std::abs(system.interaction().chi(1,1)) < eps);

      DArray<double> vals = simulator.chiEvals();
      TEST_ASSERT(std::abs((vals[0] - simulator.chiEval(0))/chi) < eps);
      TEST_ASSERT(std::abs((vals[1] - simulator.chiEval(1))/chi) < eps);
      TEST_ASSERT(std::abs((vals[0] + chi)/chi) < eps);
      TEST_ASSERT(std::abs(vals[1]/chi) < eps);

      DMatrix<double> vecs = simulator.chiEvecs();
      TEST_ASSERT(std::abs(vecs(0,0) - 1.0) < eps);
      TEST_ASSERT(std::abs(vecs(0,1) + 1.0) < eps);
      TEST_ASSERT(std::abs(vecs(1,0) - 1.0) < eps);
      TEST_ASSERT(std::abs(vecs(1,1) - 1.0) < eps);

      DArray<double>  sc = simulator.sc();
      TEST_ASSERT(std::abs((sc[0] - simulator.sc(0))/chi) < eps);
      TEST_ASSERT(std::abs((sc[1] - simulator.sc(1))/chi) < eps);
      TEST_ASSERT(std::abs(simulator.sc(0)/chi) < eps);
      TEST_ASSERT(std::abs(simulator.sc(1)/chi - 0.5) < eps);

      #if 0
      std::cout << std::endl;
      std::cout << "vals  = " << vals[0] << "  " << vals[1]
                << std::endl;
      std::cout << "vec0  = " << vecs(0,0) << "  " << vecs(0,1)
                << std::endl;
      std::cout << "vec1  = " << vecs(1,0) << "  " << vecs(1,1)
                << std::endl;
      #endif

   }

   /*
   * A LinearRamp that changes all chi entries of a ternary blend of
   * homopolymers from 8 to 80 in one step. The complete new matrix is
   * valid, but the matrix after the change of only its first entry has
   * a positive eigenvalue of the projected chi matrix.
   */
   void testSimultaneousChiRamp()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSimultaneousChiRamp.log");

      Rp::System<1,CPT> system1;
      std::istringstream in(
         "System{ Mixture{ nMonomer 3 monomers[ 1 1 1 ] nPolymer 3"
         " Polymer{ type linear nBlock 1 blocks[ 0 1 ]"
         " phi 0.3333333333333333 }"
         " Polymer{ type linear nBlock 1 blocks[ 1 1 ]"
         " phi 0.3333333333333333 }"
         " Polymer{ type linear nBlock 1 blocks[ 2 1 ]"
         " phi 0.3333333333333334 }"
         " ds 0.05 } Interaction{ chi( 0 1 8 0 2 8 1 2 8 ) }"
         " Domain{ mesh 8 lattice lamellar } }\n");
      system1.readParam(in);
      Rp::BdSimulator<1,CPT> simulator(system1);
      std::istringstream simulatorIn("BdSimulator{ seed 1234 }\n");
      simulator.readParam(simulatorIn);
      simulator.analyzeChi();

      Rp::LinearRamp<1,CPT> ramp(simulator);
      std::istringstream rampIn(
         "LinearRamp{ nParameter 3 parameters["
         " chi 0 1 72 chi 0 2 72 chi 1 2 72 ] }\n");
      ramp.readParam(rampIn);
      ramp.setup(1);
      bool success = true;
      try {
         ramp.setParameters(1);
      } catch (Exception& e) {
         success = false;
      }
      TEST_ASSERT(success);
      for (int i = 0; i < 3; ++i) {
         for (int j = i + 1; j < 3; ++j) {
            TEST_ASSERT(std::abs(system1.interaction().chi(i,j) - 80.0)
                        < 1.0E-12);
         }
      }
      for (int i = 0; i < 2; ++i) {
         TEST_ASSERT(std::abs(simulator.chiEval(i) + 80.0) < 1.0E-10);
      }
   }

   void testSaddlePointField()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testSaddlePointField.log");

      initSystem("in/param_Simulator");
      Rp::Simulator<3,CPT> simulator(system);
      simulator.allocate();
      simulator.analyzeChi();

      system.w().readRGrid("in/w_gyr.rf");
      DArray< RField<3,CPT> > const & w = system.w().rgrid();

      system.compute();
      DArray< RField<3,CPT> > const & c = system.c().rgrid();

      int nMonomer = system.mixture().nMonomer();
      int meshSize = system.domain().mesh().size();

      simulator.computeWc();
      DArray< RField<3,CPT> > const & wc = simulator.wc();

      simulator.computeCc();
      DArray< RField<3,CPT> > const & cc = simulator.cc();

      simulator.computeDc();
      DArray< RField<3,CPT> > const & dc = simulator.dc();

      // Check allocation and capacities
      TEST_ASSERT(c.capacity() == nMonomer);
      TEST_ASSERT(w.capacity() == nMonomer);
      TEST_ASSERT(wc.capacity() == nMonomer);
      TEST_ASSERT(cc.capacity() == nMonomer);
      int i;
      for (i=0; i < nMonomer; ++i) {
         TEST_ASSERT(c[i].capacity() == meshSize);
         TEST_ASSERT(w[i].capacity() == meshSize);
         TEST_ASSERT(wc[i].capacity() == meshSize);
         TEST_ASSERT(cc[i].capacity() == meshSize);
      }
      TEST_ASSERT(dc.capacity() == nMonomer - 1);
      for (i=0; i < nMonomer - 1; ++i) {
         TEST_ASSERT(dc[i].capacity() == meshSize);
      }

      // Test wc field
      double eps = 1.0E-6;
      for (i = 0; i < meshSize; ++i) {
         TEST_ASSERT(fabs(w[0][i] - wc[0][i] - wc[1][i]) < eps);
         TEST_ASSERT(fabs(w[0][i] - w[1][i] - 2.0*wc[0][i]) < eps);
      }

      // Test cc field
      for (i = 0; i < meshSize; ++i) {
         TEST_ASSERT(fabs( c[0][i] - c[1][i] - cc[0][i] ) < eps);
      }

      // Test dc field
      double diff;
      for (i = 0; i < meshSize; ++i) {
         diff = fabs( dc[0][i] );
         TEST_ASSERT(diff < 1.0E-8);
      }

      double volume = system.domain().unitCell().volume();
      double vMonomer = system.mixture().vMonomer();
      double ratio = volume/vMonomer;

      // SCFT free energy for converged solution
      system.scft().compute();
      double fHelmholtz = system.scft().fHelmholtz();

      // FTS Hamiltonian at saddle-point
      simulator.computeHamiltonian();
      double hamiltonian = simulator.hamiltonian();

      // Compare FTS Hamiltonian to SCFT free energy
      diff = fabs((hamiltonian - ratio*fHelmholtz)/hamiltonian);
      TEST_ASSERT(diff < 1.0E-8);

      if (verbose() > 0) {
         Log::file() << "Hamiltonian difference (fractional) = "
                     << diff << std::endl;

         Log::file() << "fHelmholtz = "
                     << fHelmholtz << "  "
                     << ratio*fHelmholtz  << std::endl;

         Log::file() << "Hamiltonian = "
                     << hamiltonian/ratio << "  "
                     << hamiltonian << std::endl;
      }

   }

   void testComputeHamiltonian()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param_system_disordered");
      Rp::Simulator<3,CPT> simulator(system);

      simulator.allocate();
      simulator.analyzeChi();

      system.w().readRGrid("in/w_dis.rf");
      system.compute();
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeHamiltonian();

      double diff;
      double idealHamiltonian = simulator.idealHamiltonian();
      diff = fabs(-4784.86 - idealHamiltonian);
      TEST_ASSERT(diff < 1.0E-1);

      double fieldHamiltonian = simulator.fieldHamiltonian();
      diff = fabs(12081.8 - fieldHamiltonian);
      TEST_ASSERT(diff < 1.0E-1);

      double totalHamiltonian = simulator.hamiltonian();
      diff = fabs(7296.89 - totalHamiltonian);
      TEST_ASSERT(diff < 1.0E-1);

      if (verbose() > 0) {
         openLogFile("out/testComputeHamiltonian.log");
         Log::file() << "ideal Hamiltonian: "
                     << Dbl(idealHamiltonian, 20, 12) << std::endl;
         Log::file() << "field Hamiltonian: "
                     << Dbl(fieldHamiltonian, 20, 12) << std::endl;
         Log::file() << "total Hamiltonian: "
                     << Dbl(totalHamiltonian, 20, 12) << std::endl;
      }
   }

   void testDc()
   {
      printMethod(TEST_FUNC);

      initSystem("in/param_system_disordered");
      Rp::Simulator<3,CPT> simulator(system);

      simulator.allocate();
      simulator.analyzeChi();

      system.w().readRGrid("in/w_dis.rf");
      system.compute();
      simulator.computeWc();
      simulator.computeCc();
      simulator.computeDc();

      int nMonomer = system.mixture().nMonomer();
      IntVec<3> dimensions = system.domain().mesh().dimensions();
      DArray< RField<3,CPT> > dc0;
      dc0.allocate(nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         dc0[i].allocate(dimensions);
      }
      UnitCell<3> refUnitCell;
      system.domain().fieldIo().readFieldsRGrid("in/dc_dis.rf", dc0, 
                                                refUnitCell);


      #if 0
      // Gernerate reference dc
      std::ofstream out;
      openOutputFile("in/dc_dis.rf", out);
      system.domain().fieldIo().writeFieldsRGrid(out, simulator.dc(),
                                                 system.domain().unitCell(),
                                                 true, false, true);
      out.close();
      #endif

      RFieldComparison<3,CPT> comparison;
      comparison.compare(dc0, simulator.dc());
      TEST_ASSERT(comparison.maxDiff() < 1.0E-2);

   }

};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testAnalyzeChi)
TEST_ADD(SimulatorTest, testSimultaneousChiRamp)
TEST_ADD(SimulatorTest, testSaddlePointField)
TEST_ADD(SimulatorTest, testComputeHamiltonian)
TEST_ADD(SimulatorTest, testDc)
TEST_END(SimulatorTest)

#endif
