#ifndef RPC_MIXTURE_TEST_H
#define RPC_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/solvers/Mixture.h>
#include <rp/solvers/Polymer.h>
#include <rp/solvers/Block.h>
#include <rp/solvers/Propagator.h>

#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/WaveList.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/param/BracketPolicy.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class MixtureTest : public UnitTest
{

public:

   void setUp() 
   {  
      PolymerModel::setModel(PolymerModel::Thread); 
      BracketPolicy::set(BracketPolicy::Optional); 
   }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   template <int D>
   bool tracePath(Rp::Polymer<D,CPT> const & polymer, int is, int it)
   {
      if (is == it) return true;

      Pair<int> pair;
      int ib, id;
      bool done = false;
      while (!done) {
         // std::cout << std::endl << is;
         if (is == it) return false;
         pair = polymer.path(is, it);
         ib = pair[0];
         id = pair[1];
         if (polymer.block(ib).vertexId(id) != is) return false;
         if (id == 0) {
            is = polymer.block(ib).vertexId(1);
         } else {
            is = polymer.block(ib).vertexId(0);
         }
         if (is == it) {
           // std::cout << std::endl << is;
           // std::cout << std::endl;
           done = true;
         } else {
           if (polymer.vertex(is).size() <= 1) return false;
         }
      }
      return true;
   }

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CPT> mixture;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CPT> mixture;
    
      std::ifstream in;
      openInputFile("in/Mixture1d", in);
      mixture.readParam(in);
      in.close();

      Rp::Polymer<1,CPT>& polymer = mixture.polymer(0);
      TEST_ASSERT(tracePath(polymer, 2, 0));
   }

   void testReadParameters1DBranched()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CPT> mixture;

      std::ifstream in;
      openInputFile("in/MixtureBranched", in);
      mixture.readParam(in);
      in.close();

      /*
      * Polymer graph topology (MixtureBranched file)
      *  
      *  0       2
      *   \     /
      *    4 - 5 
      *   /     \
      *  1       3
      * 
      */

      Rp::Polymer<1,CPT>& polymer = mixture.polymer(0);
      Pair<int> pair;

      pair = polymer.path(0, 1);
      TEST_ASSERT(pair[0] == 0);
      TEST_ASSERT(pair[1] == 0);

      pair = polymer.path(5, 1);
      TEST_ASSERT(pair[0] == 4);
      TEST_ASSERT(pair[1] == 1);

      TEST_ASSERT(tracePath(polymer, 1, 3));
      TEST_ASSERT(tracePath(polymer, 2, 4));
      TEST_ASSERT(tracePath(polymer, 5, 4));

   }

   void testReadParameters1D_bead()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<1,CPT> mixture;
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
      Rp::Mixture<1,CPT> mixture;

      std::ifstream in;
      openInputFile("in/Mixture1d", in);
      mixture.readParam(in);
      UnitCell<1> unitCell;
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      in.close();

      Mesh<1> mesh;
      mesh.setDimensions(d);
      FFT<1,CPT> fft;
      fft.setup(d);

      WaveList<1,CPT> waveList;
      waveList.allocate(mesh, unitCell);

      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();

      #if 0
      std::cout << "\n";
      mixture.writeParam(std::cout);
      std::cout << "unitCell  " << unitCell << std::endl;
      std::cout << "mesh      " << mesh.dimensions() << std::endl;
      #endif

      int nMonomer = mixture.nMonomer();
      DArray< RField<1,CPT> > wFields;
      DArray< RField<1,CPT> > cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }

      double cs;
      for (int i = 0; i < nx; ++i) {
         //cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         //cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 1).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));

      #if 0
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(0,1), Q = " << q11 << "\n";
      }
      
      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      #endif
      
   }

   void testSolver1D_bead()
   {
      printMethod(TEST_FUNC);
      PolymerModel::setModel(PolymerModel::Bead); 

      // Read parameter block, unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture1d_bead", in);
      Rp::Mixture<1,CPT> mixture;
      mixture.readParam(in);
      UnitCell<1> unitCell;
      in >> unitCell;
      IntVec<1> d;
      in >> d;
      in.close();

      // Setup up associated objects and allocate
      Mesh<1> mesh;
      mesh.setDimensions(d);
      FFT<1,CPT> fft;
      fft.setup(d);
      WaveList<1,CPT> waveList;
      waveList.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();

      // Check polymer block sizes
      Rp::Polymer<1,CPT>& polymer = mixture.polymer(0);
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

      // Allocate w and field arrays
      int nMonomer = mixture.nMonomer();
      DArray< RField<1,CPT> > wFields;
      DArray< RField<1,CPT> > cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }

      // Initialize w fields
      double cs;
      for (int i = 0; i < nx; ++i) {
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }

      // Solve MDE
      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different vertices
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 1).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));

      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(0,1), Q = " << q11 << "\n";
      }
      
   }

   void testSolver2D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<2,CPT> mixture;

      // Read param file block unit cell and mesh dimensions
      std::ifstream in;
      openInputFile("in/Mixture2d", in);
      mixture.readParam(in);
      UnitCell<2> unitCell;
      in >> unitCell;
      IntVec<2> d;
      in >> d;
      in.close();

      // Create associated objects and allocate mixture
      Mesh<2> mesh;
      mesh.setDimensions(d);
      FFT<2,CPT> fft;
      fft.setup(d);
      WaveList<2,CPT> waveList;
      waveList.allocate(mesh, unitCell);
      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();

      #if 0
      std::cout << "\n";
      mixture.writeParam(std::cout);
      std::cout << "unitCell  " << unitCell << std::endl;
      std::cout << "mesh      " << mesh.dimensions() << std::endl;
      #endif

      // Allocate w and c field arrays
      int nMonomer = mixture.nMonomer();
      DArray< RField<2,CPT> > wFields;
      DArray< RField<2,CPT> > cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }

      #if 0
      double cs;
      for (int i = 0; i < nx; ++i) {
         //cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         //cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }
      #endif

      // Generate oscillatory wField
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
            wFields[0][k] = 0.5 + cx + cy;
            wFields[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }
      TEST_ASSERT(k == nx);

      // Solve MDE
      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 1).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));

      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(0,1), Q = " << q11 << "\n";
      }
      
   }

   void testSolver2D_hex()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<2,CPT> mixture;

      std::ifstream in;
      openInputFile("in/Mixture2d_hex", in);
      mixture.readParam(in);
      UnitCell<2> unitCell;
      in >> unitCell;
      IntVec<2> d;
      in >> d;
      in.close();

      Mesh<2> mesh;
      mesh.setDimensions(d);
      FFT<2,CPT> fft;
      fft.setup(d);

      WaveList<2,CPT> waveList;
      waveList.allocate(mesh, unitCell);

      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();

      #if 0
      std::cout << "\n";
      mixture.writeParam(std::cout);
      std::cout << "unitCell  " << unitCell << std::endl;
      std::cout << "mesh      " << mesh.dimensions() << std::endl;
      #endif

      int nMonomer = mixture.nMonomer();
      DArray< RField<2,CPT> > wFields;
      DArray< RField<2,CPT> > cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }

      // Generate oscillatory wField
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
            wFields[0][k] = 0.5 + cx + cy;
            wFields[1][k] = 0.5 - cx - cy;
            ++k;
         }
      }
      TEST_ASSERT(k == nx);

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 1).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));

      #if 0
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(0,1), Q = " << q11 << "\n";
      }
      
      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      #endif
  
   }
    
   void testSolver3D()
   {
      printMethod(TEST_FUNC);
      Rp::Mixture<3,CPT> mixture;

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
      FFT<3,CPT> fft;
      fft.setup(d);

      WaveList<3,CPT> waveList;
      waveList.allocate(mesh, unitCell);

      mixture.associate(mesh, fft, unitCell, waveList);
      mixture.allocate();

      #if 0
      std::cout << "\n";
      mixture.writeParam(std::cout);
      std::cout << "unitCell  " << unitCell << std::endl;
      std::cout << "mesh      " << mesh.dimensions() << std::endl;
      #endif

      int nMonomer = mixture.nMonomer();
      DArray< RField<3,CPT> > wFields;
      DArray< RField<3,CPT> > cFields;
      wFields.allocate(nMonomer);
      cFields.allocate(nMonomer);
      double nx = (double)mesh.size();
      for (int i = 0; i < nMonomer; ++i) {
         wFields[i].allocate(d);
         cFields[i].allocate(d);
      }

      double cs;
      for (int i = 0; i < nx; ++i) {
         //cs = cos(2.0*Constants::Pi*(double(i)+0.5)/nx);
         //cs = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         cs = cos(2.0*Constants::Pi*double(i)/double(nx));
         wFields[0][i] = 0.5 + cs;
         wFields[1][i] = 0.5 - cs;
      }

      mixture.compute(wFields, cFields);

      // Test if same Q is obtained from different methods
      double q00, q10, q01, q11;
      mixture.polymer(0).propagator(0, 1).computeQ(q00);
      mixture.polymer(0).propagator(0, 1).computeQ(q01);
      mixture.polymer(0).propagator(1, 0).computeQ(q10);
      mixture.polymer(0).propagator(1, 1).computeQ(q11);
      TEST_ASSERT(eq(q01, q10));
      TEST_ASSERT(eq(q01, q00));
      TEST_ASSERT(eq(q01, q11));

      #if 0
      if (verbose() > 0) {
         std::cout << "Propagator(0,0), Q = " << q00 << "\n";
         std::cout << "Propagator(1,0), Q = " << q10 << "\n";
         std::cout << "Propagator(0,1), Q = " << q01 << "\n";
         std::cout << "Propagator(0,1), Q = " << q11 << "\n";
      }
      
      // Test spatial integral of block concentration
      double sum0 = domain.spatialAverage(cFields[0]);
      double sum1 = domain.spatialAverage(cFields[1]);
      std::cout << "Volume fraction of block 0 = " << sum0 << "\n";
      std::cout << "Volume fraction of block 1 = " << sum1 << "\n";
      #endif
 
   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor1D)
TEST_ADD(MixtureTest, testReadParameters1D)
TEST_ADD(MixtureTest, testReadParameters1DBranched)
TEST_ADD(MixtureTest, testReadParameters1D_bead)
TEST_ADD(MixtureTest, testSolver1D)
TEST_ADD(MixtureTest, testSolver1D_bead)
TEST_ADD(MixtureTest, testSolver2D)
TEST_ADD(MixtureTest, testSolver2D_hex)
TEST_ADD(MixtureTest, testSolver3D)
TEST_END(MixtureTest)

#endif
