#ifndef RPC_SIMULATOR_TEST_H
#define RPC_SIMULATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/ramp/LinearRamp.h>
#include <rp/fts/ramp/RampParameter.h>
#include <rp/scft/ScftThermo.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/interaction/Interaction.h>

#include <util/tests/LogFileUnitTest.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>
#include <util/misc/Exception.h>

#include <fstream>
#include <sstream>
#include <string>

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

   /*
   * The Hamiltonian of a homogeneous blend of two homopolymers of
   * length 1 with volume fractions 1 - 5.0E-9 and 5.0E-9 must not
   * depend on a uniform shift of the fields (w = 0 or w = 20).
   */
   void testHamiltonianSmallPhi()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testHamiltonianSmallPhi.log");

      Rp::System<1,CPT> system1;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 2"
         " Polymer{ type linear nBlock 1 blocks[ 0 1.0 ] phi 0.999999995 }"
         " Polymer{ type linear nBlock 1 blocks[ 1 1.0 ] phi 5.0e-9 }"
         " vMonomer 0.01 ds 0.05 } Interaction{ chi( 1 0 10.0 ) }"
         " Domain{ mesh 8 lattice lamellar } }\n");
      system1.readParam(param);
      FSArray<double, 6> parameters;
      parameters.append(1.0);
      system1.setUnitCell(parameters);
      Rp::Simulator<1,CPT> simulator(system1);
      simulator.allocate();
      simulator.analyzeChi();
      DArray< RField<1,CPT> > w;
      w.allocate(2);
      w[0].allocate(system1.domain().mesh().dimensions());
      w[1].allocate(system1.domain().mesh().dimensions());
      double h[2];
      for (int j = 0; j < 2; ++j) {
         for (int m = 0; m < 8; ++m) {
            w[0][m] = w[1][m] = 20.0*j;
         }
         system1.w().setRGrid(w);
         system1.compute();
         simulator.clearData();
         simulator.computeWc();
         simulator.computeHamiltonian();
         h[j] = simulator.hamiltonian();
      }
      TEST_ASSERT(std::abs(h[1] - h[0]) < 1.0E-9);
   }

   void testRampParameterIdRange()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testRampParameterIdRange.log");

      // Diblock melt, 1 lattice parameter. Parameters with an index out
      // of range must be rejected by getInitial, with a message that
      // says so.
      Rp::System<1,CPT> system1;
      system1.fileMaster().setInputPrefix(filePrefix());
      std::ifstream in;
      openInputFile("in/param_system_1D_diblock_thread", in);
      system1.readParam(in);
      in.close();
      Rp::Simulator<1,CPT> simulator(system1);
      std::string bad[] = {"phi_polymer 1 0.1", "cell_param 1 0.1"};
      for (int i = 0; i < 2; ++i) {
         Rp::RampParameter<1,CPT> parameter(simulator);
         std::istringstream in(bad[i]);
         in >> parameter;
         std::string message;
         try {
            parameter.getInitial();
         } catch (Exception& e) {
            message = e.message();
         }
         TEST_ASSERT(message.find("out of range") != std::string::npos);
      }
   }

   /*
   * A ramp parameter written with operator << must be read back by
   * operator >> (as for the parameter file block of a LinearRamp
   * written by writeParam).
   */
   void testRampParameterIo()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testRampParameterIo.log");

      char const * parameters[] = {"v_monomer 0.1", "chi 0 1 0.1"};
      for (int i = 0; i < 2; ++i) {
         Rp::RampParameter<1,CPT> original, copy;
         std::istringstream in(parameters[i]);
         in >> original;
         std::ostringstream out;
         out << original;
         std::istringstream reread(out.str());
         bool success = true;
         try {
            reread >> copy;
         } catch (Exception& e) {
            success = false;
         }
         TEST_ASSERT(success);
         TEST_ASSERT(copy.type() == original.type());
         TEST_ASSERT(std::abs(copy.change() - original.change()) < 1.0E-15);
      }
   }

   void testLinearRampPhiMixed()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearRampPhiMixed.log");

      // Diblock (closed) in an open solvent: a ramp of the volume
      // fraction of the diblock alone must be accepted
      Rp::System<1,CPT> system1;
      std::istringstream param(
         "System{ Mixture{ nMonomer 2 monomers[ 1.0 1.0 ]"
         " nPolymer 1 nSolvent 1"
         " Polymer{ type linear nBlock 2 blocks[ 0 0.5 1 0.5 ]"
         " phi 0.8 }"
         " Solvent{ monomerId 1 size 0.1 mu 1.0 }"
         " vMonomer 0.01 ds 0.01 }"
         " Interaction{ chi( 1 0 12.0 ) }"
         " Domain{ mesh 40 lattice lamellar } }\n");
      system1.readParam(param);
      Rp::Simulator<1,CPT> simulator(system1);
      Rp::LinearRamp<1,CPT> ramp(simulator);
      std::istringstream in("LinearRamp{"
                            " nParameter 1 parameters phi_polymer 0 -0.1 }\n");
      bool thrown = false;
      try {
         ramp.readParam(in);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
   }

   void testLinearRampPhiSum()
   {
      printMethod(TEST_FUNC);
      openLogFile("out/testLinearRampPhiSum.log");

      // Closed diblock melt: a change of phi of 1.0E-8 must be rejected
      Rp::System<1,CPT> system1;
      std::ifstream in;
      openInputFile("in/param_system_1D_diblock_thread", in);
      system1.readParam(in);
      in.close();
      Rp::Simulator<1,CPT> simulator(system1);
      Rp::LinearRamp<1,CPT> ramp(simulator);
      std::istringstream rampIn("LinearRamp{ nParameter 1"
                                " parameters phi_polymer 0 1.0e-8 }\n");
      std::string message;
      try {
         ramp.readParam(rampIn);
      } catch (Exception& e) {
         message = e.message();
      }
      TEST_ASSERT(message.find("add up to zero") != std::string::npos);
   }

};

TEST_BEGIN(SimulatorTest)
TEST_ADD(SimulatorTest, testAnalyzeChi)
TEST_ADD(SimulatorTest, testSaddlePointField)
TEST_ADD(SimulatorTest, testComputeHamiltonian)
TEST_ADD(SimulatorTest, testDc)
TEST_ADD(SimulatorTest, testHamiltonianSmallPhi)
TEST_ADD(SimulatorTest, testRampParameterIdRange)
TEST_ADD(SimulatorTest, testRampParameterIo)
TEST_ADD(SimulatorTest, testLinearRampPhiMixed)
TEST_ADD(SimulatorTest, testLinearRampPhiSum)
TEST_END(SimulatorTest)

#endif
