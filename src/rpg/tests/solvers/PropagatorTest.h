#ifndef RPG_PROPAGATOR_TEST_H
#define RPG_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rpg/solvers/Block.h>
#include <rpg/solvers/Propagator.h>
#include <rpg/solvers/WaveList.h>

#include <prdc/crystal/UnitCell.h>
#include <prdc/cuda/FFT.h>
#include <prdc/cuda/resources.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Rpg;

class PropagatorTest : public UnitTest
{

public:

   void setUp()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   template <int D> 
   void setupBlock(Pscf::Rpg::Block<D>& block)
   {
      block.setId(0);
      if (PolymerModel::isThread()) {
         double length = 2.0;
         block.setLength(length);
      } else {
         int nBead = 20;
         block.setNBead(nBead);
      }
      block.setMonomerId(1);
      double step = sqrt(6.0);
      block.setKuhn(step);
      return;
   }

   template <int D>
   void setupMesh(Mesh<D>& mesh) 
   {
      IntVec<D> d;
      for (int i = 0; i < D; ++i) {
         d[i] = 32;
      }
      mesh.setDimensions(d);
   }

   template <int D>
   void setupUnitCell(UnitCell<D>& unitCell, std::string fname)
   {
      std::ifstream in;
      openInputFile(fname, in);
      in >> unitCell;
      in.close();
   }

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Pscf::Rpg::Block<1> block;
   }

   void testSetup1D() // test allocate and associate methods
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      // Create wavelist
      WaveList<1> wavelist;

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.02));
      TEST_ASSERT(block.ns() == 101);
      TEST_ASSERT(block.mesh().dimensions()[0] == 32);

   }

   void testSetup2D() // test allocate and associate methods
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      // Create wavelist 
      WaveList<2> wavelist;

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);

      // Allocate block
      double ds = 0.26;
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.25));
      TEST_ASSERT(block.ns() == 9);
      TEST_ASSERT(block.mesh().dimensions()[0] == 32);
      TEST_ASSERT(block.mesh().dimensions()[1] == 32);
   }

   void testSetup3D() // test allocate and associate methods
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<3> block;
      setupBlock<3>(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh<3>(mesh);
      FFT<3> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Hexagonal");

      // Create wavelist 
      WaveList<3> wavelist;

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);

      // Allocate block
      double ds = 0.3;
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(block.ns() == 7);
      TEST_ASSERT(eq(block.ds(), 1.0/3.0));
      TEST_ASSERT(block.mesh().dimensions()[0] == 32);
      TEST_ASSERT(block.mesh().dimensions()[1] == 32);
      TEST_ASSERT(block.mesh().dimensions()[2] == 32);
   }
   
   void testSetupSolver1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 4.0));

      // Construct wavelist
      WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<1> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      
      for (int i = 0; i < nx; ++i) {
         w_h[i] = 1.0;
      }

      w = w_h; // copy to device

      block.setupSolver(w);
   }

   void testSetupSolver2D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Construct wavelist
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<2> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      
      for (int i = 0; i < nx; ++i) {
         w_h[i] = 1.0;
      }

      w = w_h; // copy to device

      block.setupSolver(w);
   }

   void testSetupSolver2D_bead()
   {
      printMethod(TEST_FUNC);

      PolymerModel::setModel(PolymerModel::Bead);

      // Create and initialize block
      Pscf::Rpg::Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      // Create wavelist 
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);

      double ds = 0.5;
      block.associate(mesh, fft, unitCell, wavelist);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Setup chemical potential field
      RField<2> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      VecOp::eqS(w, 1.0);
      //for (int i=0; i < w.capacity(); ++i) {
      //   w[i] = 1.0;
      //}

      block.clearUnitCellData();
      block.setupSolver(w);
   }

   void testSetupSolver3D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<3> block;
      setupBlock<3>(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh<3>(mesh);
      FFT<3> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Orthorhombic");

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Construct wavelist
      WaveList<3> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<3> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      
      for (int i = 0; i < nx; ++i) {
         w_h[i] = 1.0;
      }

      w = w_h; // copy to device

      block.setupSolver(w);
   }

   void testSolver1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 4.0));

      // Construct wavelist
      WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<1> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      double wc = 0.3;
      for (int i = 0; i < nx; ++i) {
         w_h[i] = wc;
      }

      w = w_h; // copy to device

      block.setupSolver(w);

      // Setup fields on host and device
      Propagator<1>::QField d_qin, d_qout;
      d_qin.allocate(mesh.dimensions());
      d_qout.allocate(mesh.dimensions());
      HostDArray<cudaReal> qin(nx);
      HostDArray<cudaReal> qout(nx);

      // Run block step
      double twoPi = 2.0*Constants::Pi;
      for (int i = 0; i < nx; ++i) {
         qin[i] = cos(twoPi*double(i)/double(nx));
      }

      d_qin = qin;
      block.stepThread(d_qin, d_qout);
      qout = d_qout;

      // Test block step output against expected output
      double a = 4.0;
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      ds = block.ds();
      double expected = exp(-(wc + r)*ds);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(qout[i], qin[i]*expected));
      }

      // Test propagator solve 
      block.propagator(0).solve();

      // Copy results from propagator solve
      HostDArray<cudaReal> propHead(nx);
      HostDArray<cudaReal> propTail(nx);
      propHead = block.propagator(0).head();
      propTail = block.propagator(0).tail();

      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(propHead[i],1.0));
      }
      
      expected = exp(-wc*block.length());
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(propTail[i], expected));
      }

   }

   void testSolver1D_bead()
   {
      printMethod(TEST_FUNC);

      PolymerModel::setModel(PolymerModel::Bead);

      // Create and initialize block
      Pscf::Rpg::Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      WaveList<1> wavelist;
      wavelist.allocate(mesh, unitCell);

      double ds = 1.00;
      block.associate(mesh, fft, unitCell, wavelist);
      bool own0 = true;
      bool own1 = true;
      block.setVertexOwnership(own0, own1);
      block.allocate(ds);
      TEST_ASSERT(block.ns() == block.nBead());

      // Setup chemical potential field
      RField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wc = 0.3;
      VecOp::eqS(w, wc);
      //for (int i=0; i < nx; ++i) {
      //   w[i] = wc;
      //}

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<1>::QField qin;
      Propagator<1>::QField qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      // Initialize qin as a cosine(kx)
      HostDArray<cudaReal> qin_h(nx);
      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         qin_h[i] = cos(twoPi*double(i)/double(nx));
      }
      qin = qin_h;

      // Take a step
      block.stepBead(qin, qout);
      HostDArray<cudaReal> qout_h(nx);
      qout_h = qout;

      // Test qout_h
      //double a = 4.0;
      double a = unitCell.parameter(0);
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      double expected = exp(-(wc + r));
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(qout_h[i], qin_h[i]*expected));
      }
   
      block.propagator(0).setVertexOwnership(own0, own1);
      block.propagator(1).setVertexOwnership(own1, own0);
      TEST_ASSERT(block.propagator(0).ownsHead() == own0);
      TEST_ASSERT(block.propagator(0).ownsTail() == own1);
      TEST_ASSERT(block.propagator(1).ownsHead() == own1);
      TEST_ASSERT(block.propagator(1).ownsTail() == own0);

      // Test propagator solve, block owns both vertices
      Propagator<1>& q0 = block.propagator(0);
      q0.solve();

      // Check head slice
      HostDArray<cudaReal> qh_h(nx);
      qh_h = q0.head();
      expected = exp(-wc*ds);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(qh_h[i], expected));
      }
      
      // Check tail slice
      HostDArray<cudaReal> qt_h(nx);
      qt_h = q0.tail();
      expected = exp(-wc*ds*block.nBead());
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(qt_h[i], expected));
      }

      Propagator<1>& q1 = block.propagator(1);
      q1.solve();

      double qh = block.averageProductBead(q0.head(), q1.tail());
      double qt = block.averageProductBead(q0.tail(), q1.head());
      TEST_ASSERT(eq( log(qh), log(qt) )) ;
      TEST_ASSERT(eq( log(qh), -wc*ds*block.nBead() )) ;
   }

   void testSolver2D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Construct wavelist
      WaveList<2> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<2> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      double wc = 0.3;
      for (int i = 0; i < nx; ++i) {
         w_h[i] = wc;
      }

      w = w_h; // copy to device

      block.setupSolver(w);

      // Setup fields on host and device
      Propagator<2>::QField d_qin, d_qout;
      d_qin.allocate(mesh.dimensions());
      d_qout.allocate(mesh.dimensions());
      HostDArray<cudaReal> qin(nx);
      HostDArray<cudaReal> qout(nx);

      // Run block step
      MeshIterator<2> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) ) );
      }

      d_qin = qin;
      block.stepThread(d_qin, d_qout);
      qout = d_qout;

      // Test block step output against expected output
      double b = block.kuhn();
      double Gb;
      double expected;
      IntVec<2> temp;
      temp[0] = 1;
      temp[1] = 1;

      ds = block.ds();
      for (iter.begin(); !iter.atEnd(); ++iter){
         Gb = unitCell.ksq(temp);
         double factor = b;
         double r = Gb*factor*factor/6.0;
         expected = exp(-(wc + r)*ds);
         TEST_ASSERT(eq(qout[iter.rank()], qin[iter.rank()]*expected));
      }

      // Test propagator solve 
      block.propagator(0).solve();

      // Copy results from propagator solve
      HostDArray<cudaReal> propHead(nx);
      HostDArray<cudaReal> propTail(nx);
      propHead = block.propagator(0).head();
      propTail = block.propagator(0).tail();

      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(propHead[iter.rank()], 1.0));
      }
      
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(propTail[iter.rank()], expected));
      }
   
   }

   void testSolver3D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Pscf::Rpg::Block<3> block;
      setupBlock<3>(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh<3>(mesh);
      FFT<3> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Orthorhombic");
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Construct wavelist
      WaveList<3> wavelist;
      wavelist.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, wavelist);      

      // Allocate block
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      int nx = mesh.size();
      RField<3> w;
      w.allocate(mesh.dimensions());
      HostDArray<cudaReal> w_h(nx);

      TEST_ASSERT(w.capacity() == mesh.size());
      double wc = 0.3;
      for (int i = 0; i < nx; ++i) {
         w_h[i] = wc;
      }

      w = w_h; // copy to device

      block.setupSolver(w);

      // Setup fields on host and device
      Propagator<3>::QField d_qin, d_qout;
      d_qin.allocate(mesh.dimensions());
      d_qout.allocate(mesh.dimensions());
      HostDArray<cudaReal> qin(nx);
      HostDArray<cudaReal> qout(nx);

      // Run block step
      MeshIterator<3> iter(mesh.dimensions());
      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) + 
                         double(iter.position(2))/double(mesh.dimension(2)) ) );
      }

      d_qin = qin;
      block.stepThread(d_qin, d_qout);
      qout = d_qout;

      // Test block step output against expected output
      double b = block.kuhn();
      double Gb;
      double expected;
      IntVec<3> temp;
      temp[0] = 1;
      temp[1] = 1;
      temp[2] = 1;


      ds = block.ds();
      for (iter.begin(); !iter.atEnd(); ++iter){
         Gb = unitCell.ksq(temp);
         double factor = b;
         double r = Gb*factor*factor/6.0;
         expected = exp(-(wc + r)*ds);
         TEST_ASSERT(eq(qout[iter.rank()], qin[iter.rank()]*expected));
      }

      // Test propagator solve 
      block.propagator(0).solve();

      // Copy results from propagator solve
      HostDArray<cudaReal> propHead(nx);
      HostDArray<cudaReal> propTail(nx);
      propHead = block.propagator(0).head();
      propTail = block.propagator(0).tail();

      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(propHead[iter.rank()], 1.0));
      }
      
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(propTail[iter.rank()], expected));
      }
   
   }
};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor1D)
TEST_ADD(PropagatorTest, testSetup1D)
TEST_ADD(PropagatorTest, testSetup2D)
TEST_ADD(PropagatorTest, testSetup3D)
TEST_ADD(PropagatorTest, testSetupSolver1D)
TEST_ADD(PropagatorTest, testSetupSolver2D)
TEST_ADD(PropagatorTest, testSetupSolver2D_bead)
TEST_ADD(PropagatorTest, testSetupSolver3D)
TEST_ADD(PropagatorTest, testSolver1D)
TEST_ADD(PropagatorTest, testSolver1D_bead)
TEST_ADD(PropagatorTest, testSolver2D)
TEST_ADD(PropagatorTest, testSolver3D)
TEST_END(PropagatorTest)

#endif
