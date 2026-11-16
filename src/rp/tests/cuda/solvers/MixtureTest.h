#ifndef RPG_MIXTURE_TEST_H
#define RPG_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Solvent.h>
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>

#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/Reduce.h>

#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/FFT.h>
#include <pscf/field/cuda/WaveList.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>
#include <pscf/backend/VecOp.h>

#include <util/math/Constants.h>


#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace Util;
using namespace Pscf;

class MixtureTest : public UnitTest
{

public:

   void setUp()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   /*
   * Read a 1D mixture, unit cell and mesh, set up and allocate the
   * mixture, and allocate arrays of w and c fields.
   */
   void setupMixture1D(std::istream& in, Rp::Mixture<1,CUT>& mixture,
                       UnitCell<1>& unitCell, Mesh<1>& mesh,
                       FFT<1,CUT>& fft, WaveList<1,CUT>& waveList,
                       DArray< RField<1,CUT> >& wFields,
                       DArray< RField<1,CUT> >& cFields)
   {
      mixture.readParam(in);
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      mesh.setDimensions(d);
      fft.setup(d);
      waveList.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();
      int nMonomer = mixture.nMonomer();
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }
   }

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CUT> mixture;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CUT> mixture;

      std::ifstream in;
      openInputFile("in/Mixture1d", in);
      mixture.readParam(in);
      in.close();
   }

   void testReadParameters1D_bead()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CUT> mixture;
      PolymerModel::setModel(PolymerModel::Bead);

      std::ifstream in;
      openInputFile("in/Mixture1d_bead", in);
      mixture.readParam(in);
      in.close();

      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testSolver1D()
   {
      printMethod(TEST_FUNC);

      Rp::Mixture<1,CUT> mixture;
      Mesh<1> mesh;
      FFT<1,CUT> fft;
      UnitCell<1> unitCell;
      WaveList<1,CUT> wavelist;
      IntVec<1> d;

      // Read parameter block, unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture1d", in);
      mixture.readParam(in);
      in >> unitCell;
      in >> d;
      in.close();

      // Set up objects
      mesh.setDimensions(d);
      fft.setup(d);
      wavelist.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      // Allocate w and c field arrays
      int nMonomer = mixture.nMonomer();
      DArray< RField<1,CUT> > wFields;
      DArray< RField<1,CUT> > cFields;
      DArray< HostSendArray<cudaReal,CUT> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields_h[0][i] = 0.5 + cs;
         wFields_h[1][i] = 0.5 - cs;
      }
      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 0).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(1,1), Q = " << q11 << "\n";
      }

   }

   void testSolver1D_bead()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead); 

      // Read parameter block, unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture1d_bead", in);
      Rp::Mixture<1,CUT> mixture;
      mixture.readParam(in);
      UnitCell<1> unitCell;
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      in.close();

      // Set up associated objects and allocate
      Mesh<1> mesh;
      mesh.setDimensions(d);
      FFT<1,CUT> fft;
      fft.setup(d);
      WaveList<1,CUT> wavelist;
      wavelist.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      // Check polymer block sizes
      Rp::Polymer<1,CUT>& polymer = mixture.polymer(0);
      TEST_ASSERT(polymer.block(0).nBead() == 20);
      TEST_ASSERT(polymer.block(1).nBead() == 30);
      TEST_ASSERT(polymer.nBead() == 50);

      // Check end flags for a diblock
      TEST_ASSERT(polymer.block(0).propagator(0).isHeadEnd());
      TEST_ASSERT(!polymer.block(0).propagator(0).isTailEnd());
      TEST_ASSERT(!polymer.block(0).propagator(1).isHeadEnd());
      TEST_ASSERT(polymer.block(0).propagator(1).isTailEnd());
      TEST_ASSERT(!polymer.block(1).propagator(0).isHeadEnd());
      TEST_ASSERT(polymer.block(1).propagator(0).isTailEnd());
      TEST_ASSERT(polymer.block(1).propagator(1).isHeadEnd());
      TEST_ASSERT(!polymer.block(1).propagator(1).isTailEnd());

      // Allocate w and c field arrays
      int nMonomer = mixture.nMonomer();
      DArray< RField<1,CUT> > wFields;
      DArray< RField<1,CUT> > cFields;
      DArray< HostSendArray<cudaReal,CUT> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      // Initialize w fields
      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields_h[0][i] = 0.5 + cs;
         wFields_h[1][i] = 0.5 - cs;
      }
      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      // Solve MDE
      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 0).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(1,1), Q = " << q11 << "\n";
      }


   }

   /*
   * Bead chains in fields equal to h at some grid points and 0 at the
   * others: a homopolymer of one bead of type 0 (phi 0.25) and a
   * diblock of 4 + 4 beads of types 0 and 1 (phi 0.75), with h at the
   * odd grid points for type 0, and at grid points 0 and 4 for type 1.
   * For h = 2000, the reciprocal Boltzmann weight exp(h - s) of a bead
   * at those points (s the field shift) is above the range of a double.
   * The concentrations must be finite, those of the homopolymer must be
   * phi exp(-w)/<exp(-w)>, those of each block must add up to its share
   * of phi, and they must agree with those for h = 300 (for which the
   * weight exp(-300) is negligible).
   */
   void testSolverBeadStrongField()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead);
      double h[2] = {2000.0, 300.0};
      double c[2][2][8];
      for (int k = 0; k < 2; ++k) {
         Rp::Mixture<1,CUT> mixture;
         UnitCell<1> unitCell;
         Mesh<1> mesh;
         FFT<1,CUT> fft;
         WaveList<1,CUT> waveList;
         DArray< RField<1,CUT> > wFields;
         DArray< RField<1,CUT> > cFields;
         std::istringstream in(
            "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 2 "
            "Polymer{ type linear nBlock 1 blocks[ 0 1 ] phi 0.25 } "
            "Polymer{ type linear nBlock 2 blocks[ 0 4 1 4 ] phi 0.75 } "
            "} lamellar 2.0 8");
         setupMixture1D(in, mixture, unitCell, mesh, fft, waveList,
                        wFields, cFields);
         int nx = mesh.size();
         HostSendArray<cudaReal,CUT> w_h(nx);
         HostRecvArray<cudaReal,CUT> c_h(nx);
         for (int i = 0; i < nx; ++i) {
            w_h[i] = (i % 2 == 1) ? h[k] : 0.0;
         }
         wFields[0] = w_h;
         for (int i = 0; i < nx; ++i) {
            w_h[i] = (i % 4 == 0) ? h[k] : 0.0;
         }
         wFields[1] = w_h;
         mixture.compute(wFields, cFields);

         c_h = mixture.polymer(0).block(0).cField();
         for (int i = 0; i < nx; ++i) {
            double cRef = (i % 2 == 0) ? 0.5 : 0.0;
            TEST_ASSERT(std::abs(c_h[i] - cRef) < 1.0E-12);
         }
         for (int j = 0; j < 2; ++j) {
            RField<1,CUT> const & cBlock
                                = mixture.polymer(1).block(j).cField();
            double sum = Reduce::sum(cBlock)/double(nx);
            TEST_ASSERT(std::abs(sum - 0.375) < 1.0E-10);
         }
         for (int m = 0; m < 2; ++m) {
            c_h = cFields[m];
            for (int i = 0; i < nx; ++i) {
               c[k][m][i] = c_h[i];
               TEST_ASSERT(std::isfinite(c[k][m][i]));
            }
         }
      }
      for (int m = 0; m < 2; ++m) {
         for (int i = 0; i < 8; ++i) {
            TEST_ASSERT(std::abs(c[0][m][i] - c[1][m][i]) < 1.0E-12);
         }
      }
   }

   void testSolver2D()
   {
      printMethod(TEST_FUNC);

      Rp::Mixture<2,CUT> mixture;
      Mesh<2> mesh;
      FFT<2,CUT> fft;
      UnitCell<2> unitCell;
      WaveList<2,CUT> wavelist;
      IntVec<2> d;

      // Read parameter block, unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture2d", in);
      mixture.readParam(in);
      in >> unitCell;
      in >> d;
      in.close();

      mesh.setDimensions(d);
      fft.setup(d);
      wavelist.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      int nMonomer = mixture.nMonomer();
      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!

      // Allocate w and c field arrays on device and host
      DArray< RField<2,CUT> > wFields;
      DArray< RField<2,CUT> > cFields;
      DArray< HostSendArray<cudaReal,CUT> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      // Generate oscillatory w field
      int dx = mesh.dimension(0);
      int dy = mesh.dimension(1);
      double fx = 2.0*Constants::Pi/double(dx);
      double fy = 2.0*Constants::Pi/double(dy);
      double cx, cy;
      int k = 0;
      for (int i = 0; i < dx; ++i) {
         cx = cos(fx*double(i));
         for (int j = 0; j < dy; ++j) {
            cy = cos(fy*double(j));
            wFields_h[0][k] = 0.5 + cx + cy;
            wFields_h[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }
      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      // Perform computation (solve MDE and compute concentrations)
      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 0).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(1,1), Q = " << q11 << "\n";
      }

   }

   void testSolver2D_hex()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<2,CUT> mixture;
      Mesh<2> mesh;
      FFT<2,CUT> fft;
      WaveList<2,CUT> wavelist;
      IntVec<2> d;

      // Read file: param block, unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture2d_hex", in);
      mixture.readParam(in);
      UnitCell<2> unitCell;
      in >> unitCell;
      in >> d;
      in.close();

      mesh.setDimensions(d);
      fft.setup(d);
      wavelist.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      int nMonomer = mixture.nMonomer();
      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!

      // Allocate w and c field arrays on device and host
      DArray< RField<2,CUT> > wFields;
      DArray< RField<2,CUT> > cFields;
      DArray< HostSendArray<cudaReal,CUT> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      // Generate oscillatory w field on host
      int dx = mesh.dimension(0);
      int dy = mesh.dimension(1);
      double fx = 2.0*Constants::Pi/double(dx);
      double fy = 2.0*Constants::Pi/double(dy);
      double cx, cy;
      int k = 0;
      for (int i = 0; i < dx; ++i) {
         cx = cos(fx*double(i));
         for (int j = 0; j < dy; ++j) {
            cy = cos(fy*double(j));
            wFields_h[0][k] = 0.5 + cx + cy;
            wFields_h[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }

      // Copy w fields to device from host
      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      // Perform computation
      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 0).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(1,1), Q = " << q11 << "\n";
      }

   }

   void testSolver3D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<3,CUT> mixture;

      std::ifstream in;
      openInputFile("in/Mixture3d", in);
      mixture.readParam(in);

      UnitCell<3> unitCell;
      in >> unitCell;

      IntVec<3> d;
      in >> d;
      in.close();

      Mesh<3> mesh;
      mesh.setDimensions(d);
      FFT<3,CUT> fft;
      fft.setup(d);

      // Construct wavelist
      WaveList<3,CUT> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Set up mixture
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      int nMonomer = mixture.nMonomer();
      DArray< RField<3,CUT> > wFields;
      DArray< RField<3,CUT> > cFields;
      DArray< HostSendArray<cudaReal,CUT> > wFields_h;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      wFields_h.allocate(nMonomer);
      int nx = mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
         wFields_h[i].allocate(nx);
      }

      UTIL_CHECK(nMonomer == 2); // Hard-coded in here!
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields_h[0][i] = 0.5 + cs;
         wFields_h[1][i] = 0.5 - cs;
      }

      wFields[0] = wFields_h[0];
      wFields[1] = wFields_h[1];

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 0).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(1,1), Q = " << q11 << "\n";
      }

   }

   /*
   * Solvents of size 1800 in a uniform field 0.45, whose partition
   * function q = exp(-810) is below the range of a double: an open
   * solvent of monomer type 0 with mu = 809.5, for which exp(mu) is
   * above that range, and a closed solvent of monomer type 1 with
   * phi = 0.4. Then the closed solvent in a field alternating between
   * -0.05 and 0.95, which changes its Boltzmann weight by exp(1800).
   */
   void testSolverLargeSolvent()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CUT> mixture;
      UnitCell<1> unitCell;
      Mesh<1> mesh;
      FFT<1,CUT> fft;
      WaveList<1,CUT> waveList;
      DArray< RField<1,CUT> > wFields;
      DArray< RField<1,CUT> > cFields;
      std::istringstream in(
         "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1 nSolvent 2 "
         "Polymer{ type linear nBlock 1 blocks[ 0 1.0 ] phi 0.1 } "
         "Solvent{ monomerId 0 size 1800.0 mu 809.5 } "
         "Solvent{ monomerId 1 size 1800.0 phi 0.4 } ds 0.01 } "
         "lamellar 1.0 8");
      setupMixture1D(in, mixture, unitCell, mesh, fft, waveList,
                     wFields, cFields);
      int nx = mesh.size();
      HostSendArray<cudaReal,CUT> w_h(nx);
      HostRecvArray<cudaReal,CUT> c_h(nx);

      for (int i = 0; i < nx; ++i) {
         w_h[i] = 0.45;
      }
      wFields[0] = w_h;
      wFields[1] = w_h;
      mixture.compute(wFields, cFields);
      Rp::Solvent<1,CUT> const & open = mixture.solvent(0);
      Rp::Solvent<1,CUT> const & closed = mixture.solvent(1);
      TEST_ASSERT(std::abs(open.phi() - std::exp(-0.5)) < 1.0E-10);
      TEST_ASSERT(std::abs(closed.mu() - std::log(0.4) - 810.0) < 1.0E-10);
      c_h = cFields[0];
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::abs(c_h[i] - 0.1 - std::exp(-0.5)) < 1.0E-10);
      }
      c_h = cFields[1];
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(std::abs(c_h[i] - 0.4) < 1.0E-10);
      }

      // Half of the grid points have the field -0.05, so q = exp(90)/2
      for (int i = 0; i < nx; ++i) {
         w_h[i] = (i % 2 == 0) ? -0.05 : 0.95;
      }
      wFields[1] = w_h;
      mixture.compute(wFields, cFields);
      TEST_ASSERT(std::abs(closed.mu() - std::log(0.4) + 90.0
                           - std::log(2.0)) < 1.0E-10);
      c_h = cFields[1];
      for (int i = 0; i < nx; ++i) {
         double c = (i % 2 == 0) ? 0.8 : 0.0;
         TEST_ASSERT(std::abs(c_h[i] - c) < 1.0E-10);
      }
   }

   /*
   * Logarithm of the partition function q of a linear polymer (thread
   * model, or bead model with blocks of equal kuhn lengths), evaluated
   * without field shifts: the blocks are set up without a shift, and
   * their steps are applied in order to a uniform head, with the result
   * rescaled after every step.
   */
   double lnQSteps(Rp::Polymer<1,CUT>& polymer,
                   DArray< RField<1,CUT> > const & wFields)
   {
      RField<1,CUT> q, qNew;
      q.allocate(wFields[0].meshDimensions());
      qNew.allocate(wFields[0].meshDimensions());
      VecOp::eqS(q, 1.0);
      double lnQ = 0.0;
      for (int j = 0; j < polymer.nBlock(); ++j) {
         Rp::Block<1,CUT>& block = polymer.block(j);
         block.setupSolver(wFields[block.monomerId()]);
         if (PolymerModel::isThread()) {
            for (int k = 1; k < block.ns(); ++k) {
               block.stepThread(q, qNew);
               double scale = Reduce::max(qNew);
               VecOp::mulVS(q, qNew, 1.0/scale);
               lnQ += std::log(scale);
            }
         } else {
            // One bond before each bead, except the first of the chain
            for (int k = 0; k < block.nBead(); ++k) {
               if (j == 0 && k == 0) {
                  VecOp::eqV(qNew, q);
                  block.stepFieldBead(qNew);
               } else {
                  block.stepBead(q, qNew);
               }
               double scale = Reduce::max(qNew);
               VecOp::mulVS(q, qNew, 1.0/scale);
               lnQ += std::log(scale);
            }
         }
      }
      return lnQ + std::log(Reduce::sum(q)/double(q.capacity()));
   }

   /*
   * Homopolymers in a uniform field 0.45, whose partition functions
   * are below the range of a double: thread chains of lengths 1800
   * (closed, phi = 0.1, q = exp(-810)) and 1790 (open, mu = 805.0),
   * and bead chains of 1700 beads (closed, q = exp(-765)) and 1690
   * beads (open, mu = 760.0). Each open species has phi = exp(-0.5).
   */
   void testSolverLongChain()
   {
      printMethod(TEST_FUNC);
      std::string texts[2] = {
         "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 2 nSolvent 1 "
         "Polymer{ type linear nBlock 1 blocks[ 0 1800.0 ] phi 0.1 } "
         "Polymer{ type linear nBlock 1 blocks[ 0 1790.0 ] mu 805.0 } "
         "Solvent{ monomerId 1 size 1.0 phi 0.9 } ds 1.0 } lamellar 1.0 8",
         "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 2 nSolvent 1 "
         "Polymer{ type linear nBlock 1 blocks[ 0 1700 ] phi 0.1 } "
         "Polymer{ type linear nBlock 1 blocks[ 0 1690 ] mu 760.0 } "
         "Solvent{ monomerId 1 size 1.0 phi 0.9 } } lamellar 1.0 8"};
      double lnQ0[2] = {-810.0, -765.0};
      double lnQ1[2] = {-805.5, -760.5};
      for (int k = 0; k < 2; ++k) {
         if (k == 0) {
            PolymerModel::setModel(PolymerModel::Thread);
         } else {
            PolymerModel::setModel(PolymerModel::Bead);
         }
         Rp::Mixture<1,CUT> mixture;
         UnitCell<1> unitCell;
         Mesh<1> mesh;
         FFT<1,CUT> fft;
         WaveList<1,CUT> waveList;
         DArray< RField<1,CUT> > wFields;
         DArray< RField<1,CUT> > cFields;
         std::istringstream in(texts[k]);
         setupMixture1D(in, mixture, unitCell, mesh, fft, waveList,
                        wFields, cFields);
         int nx = mesh.size();
         HostRecvArray<cudaReal,CUT> c_h(nx);

         VecOp::eqS(wFields[0], 0.45);
         VecOp::eqS(wFields[1], 0.0);
         mixture.compute(wFields, cFields);
         Rp::Polymer<1,CUT> const & closed = mixture.polymer(0);
         Rp::Polymer<1,CUT> const & open = mixture.polymer(1);
         TEST_ASSERT(std::abs(closed.lnQ() - lnQ0[k]) < 1.0E-9);
         TEST_ASSERT(std::abs(closed.mu() - std::log(0.1) + lnQ0[k])
                     < 1.0E-9);
         TEST_ASSERT(std::abs(open.lnQ() - lnQ1[k]) < 1.0E-9);
         TEST_ASSERT(std::abs(open.phi() - std::exp(-0.5)) < 1.0E-9);
         c_h = cFields[0];
         for (int i = 0; i < nx; ++i) {
            TEST_ASSERT(std::abs(c_h[i] - 0.1 - std::exp(-0.5)) < 1.0E-9);
         }
         c_h = cFields[1];
         for (int i = 0; i < nx; ++i) {
            TEST_ASSERT(std::abs(c_h[i] - 0.9) < 1.0E-12);
         }
      }
   }

   /*
   * Thread chains in fields that alternate between neighbouring grid
   * points, which the chains cannot follow, so that the shifted
   * partition functions are below the range of a double and the shifts
   * must be raised: polymer 0 is a homopolymer of length 1000 in a field
   * alternating between 0 and 2 (monomer type 2), and polymer 1 is a
   * diblock 1000 + 1000 with w_A = 0 (type 0) and w_B alternating
   * between 0 and 6 (type 1). The values of ln q must agree with an
   * evaluation without shifts, and the concentration of each block must
   * add up to its share of phi (with the shift of the B block applied to
   * the A block, the A propagators would overflow). Then the same for a
   * homopolymer of 2000 beads in a field alternating between 0 and 1,
   * whose propagator ends at a chain end.
   */
   void testSolverRaisedShift()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CUT> mixture;
      UnitCell<1> unitCell;
      Mesh<1> mesh;
      FFT<1,CUT> fft;
      WaveList<1,CUT> waveList;
      DArray< RField<1,CUT> > wFields;
      DArray< RField<1,CUT> > cFields;
      std::istringstream in(
         "Mixture{ nMonomer 3 monomers[ 1.0 1.0 1.0 ] nPolymer 2 "
         "Polymer{ type linear nBlock 1 blocks[ 2 1000.0 ] phi 0.5 } "
         "Polymer{ type linear nBlock 2 blocks[ 0 1000.0 1 1000.0 ] "
         "phi 0.5 } ds 1.0 } lamellar 2.0 8");
      setupMixture1D(in, mixture, unitCell, mesh, fft, waveList,
                     wFields, cFields);
      int nx = mesh.size();
      HostSendArray<cudaReal,CUT> w_h(nx);

      VecOp::eqS(wFields[0], 0.0);
      for (int i = 0; i < nx; ++i) {
         w_h[i] = (i % 2 == 0) ? 0.0 : 6.0;
      }
      wFields[1] = w_h;
      for (int i = 0; i < nx; ++i) {
         w_h[i] = (i % 2 == 0) ? 0.0 : 2.0;
      }
      wFields[2] = w_h;
      mixture.compute(wFields, cFields);

      double lnQ[2];
      for (int ip = 0; ip < 2; ++ip) {
         Rp::Polymer<1,CUT> const & polymer = mixture.polymer(ip);
         lnQ[ip] = polymer.lnQ();
         for (int j = 0; j < polymer.nBlock(); ++j) {
            Rp::Block<1,CUT> const & block = polymer.block(j);
            double phi = 0.5*block.length()/polymer.length();
            double sum = Reduce::sum(block.cField())/double(nx);
            TEST_ASSERT(std::abs(sum - phi) < 1.0E-10);
         }
      }
      for (int ip = 0; ip < 2; ++ip) {
         double lnQRef = lnQSteps(mixture.polymer(ip), wFields);
         TEST_ASSERT(std::abs(lnQ[ip] - lnQRef) < 1.0E-10*std::abs(lnQRef));
      }

      // Bead model
      PolymerModel::setModel(PolymerModel::Bead);
      Rp::Mixture<1,CUT> beadMixture;
      UnitCell<1> beadCell;
      Mesh<1> beadMesh;
      FFT<1,CUT> beadFft;
      WaveList<1,CUT> beadWaveList;
      DArray< RField<1,CUT> > beadW;
      DArray< RField<1,CUT> > beadC;
      std::istringstream beadIn(
         "Mixture{ nMonomer 1 monomers[ 1.0 ] nPolymer 1 "
         "Polymer{ type linear nBlock 1 blocks[ 0 2000 ] phi 1.0 } } "
         "lamellar 2.0 8");
      setupMixture1D(beadIn, beadMixture, beadCell, beadMesh, beadFft,
                     beadWaveList, beadW, beadC);
      for (int i = 0; i < nx; ++i) {
         w_h[i] = (i % 2 == 0) ? 0.0 : 1.0;
      }
      beadW[0] = w_h;
      beadMixture.compute(beadW, beadC);
      double lnQBead = beadMixture.polymer(0).lnQ();
      TEST_ASSERT(std::abs(Reduce::sum(beadC[0])/double(nx) - 1.0) < 1.0E-10);
      double lnQRef = lnQSteps(beadMixture.polymer(0), beadW);
      TEST_ASSERT(std::abs(lnQBead - lnQRef) < 1.0E-10*std::abs(lnQRef));
   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor1D)
TEST_ADD(MixtureTest, testReadParameters1D)
TEST_ADD(MixtureTest, testReadParameters1D_bead)
TEST_ADD(MixtureTest, testSolver1D)
TEST_ADD(MixtureTest, testSolver1D_bead)
TEST_ADD(MixtureTest, testSolverBeadStrongField)
TEST_ADD(MixtureTest, testSolver2D)
TEST_ADD(MixtureTest, testSolver2D_hex)
TEST_ADD(MixtureTest, testSolver3D)
TEST_ADD(MixtureTest, testSolverLargeSolvent)
TEST_ADD(MixtureTest, testSolverLongChain)
TEST_ADD(MixtureTest, testSolverRaisedShift)
TEST_END(MixtureTest)

#endif
