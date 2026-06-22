#ifndef RPG_MIXTURE_TEST_H
#define RPG_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/solvers/Mixture.h>
#include <rpg/solvers/Polymer.h>
#include <rpg/solvers/Solvent.h>
#include <rpg/solvers/Block.h>
#include <rpg/solvers/Propagator.h>

#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/cudaTypes.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/WaveList.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/math/Constants.h>


#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class MixtureTest : public UnitTest
{

public:

   void setUp()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1, Rpg::Types<1> > mixture;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1, Rpg::Types<1> > mixture;

      std::ifstream in;
      openInputFile("in/Mixture1d", in);
      mixture.readParam(in);
      in.close();
   }

   void testReadParameters1D_bead()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1, Rpg::Types<1> > mixture;
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

      Rp::Mixture<1, Rpg::Types<1> > mixture;
      Mesh<1> mesh;
      Cuda::FFT<1> fft;
      UnitCell<1> unitCell;
      Cuda::WaveList<1> wavelist;
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
      DArray< Cuda::RField<1> > wFields;
      DArray< Cuda::RField<1> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
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
      Rp::Mixture<1, Rpg::Types<1> > mixture;
      mixture.readParam(in);
      UnitCell<1> unitCell;
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      in.close();

      // Set up associated objects and allocate
      Mesh<1> mesh;
      mesh.setDimensions(d);
      Cuda::FFT<1> fft;
      fft.setup(d);
      Cuda::WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      // Check polymer block sizes
      Rp::Polymer<1, Types<1> >& polymer = mixture.polymer(0);
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
      DArray< Cuda::RField<1> > wFields;
      DArray< Cuda::RField<1> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
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

   void testSolver2D()
   {
      printMethod(TEST_FUNC);

      Rp::Mixture<2, Rpg::Types<2> > mixture;
      Mesh<2> mesh;
      Cuda::FFT<2> fft;
      UnitCell<2> unitCell;
      Cuda::WaveList<2> wavelist;
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
      DArray< Cuda::RField<2> > wFields;
      DArray< Cuda::RField<2> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
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
      Rp::Mixture<2, Rpg::Types<2> > mixture;
      Mesh<2> mesh;
      Cuda::FFT<2> fft;
      Cuda::WaveList<2> wavelist;
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
      DArray< Cuda::RField<2> > wFields;
      DArray< Cuda::RField<2> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
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
      Rp::Mixture<3, Rpg::Types<3> > mixture;

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
      Cuda::FFT<3> fft;
      fft.setup(d);

      // Construct wavelist
      Cuda::WaveList<3> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Set up mixture
      mixture.associate(mesh, fft, unitCell, wavelist);
      mixture.allocate();

      int nMonomer = mixture.nMonomer();
      DArray< Cuda::RField<3> > wFields;
      DArray< Cuda::RField<3> > cFields;
      DArray< HostDArray<cudaReal> > wFields_h;
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

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor1D)
TEST_ADD(MixtureTest, testReadParameters1D)
TEST_ADD(MixtureTest, testReadParameters1D_bead)
TEST_ADD(MixtureTest, testSolver1D)
TEST_ADD(MixtureTest, testSolver1D_bead)
TEST_ADD(MixtureTest, testSolver2D)
TEST_ADD(MixtureTest, testSolver2D_hex)
TEST_ADD(MixtureTest, testSolver3D)
TEST_END(MixtureTest)

#endif
