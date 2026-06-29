#ifndef CPC_PROPAGATOR_TEST_H
#define CPC_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cpc/solvers/Block.h>
#include <cpc/solvers/Propagator.h>
#include <cpc/field/Domain.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/WaveList.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/chem/PolymerModel.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/cpu/complex.h>

#include <util/param/BracketPolicy.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Cpc;
using namespace Pscf::Prdc;

class PropagatorTest : public UnitTest
{

public:

   void setUp()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   void tearDown()
   {  PolymerModel::setModel(PolymerModel::Thread); }

   // Utility functions (used in unit tests)

   template <int D> 
   void setupBlock(Block<D>& block)
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

   /*
   * Open and read file header to initialize Domain<D> system.
   */
   template <int D>
   int readHeader(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      int nMonomer;
      domain.readFieldHeader(in, nMonomer);
      in.close();
      return nMonomer;
   }

   // Unit test functions

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Block<1> block;
   }

   void testSetup1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      WaveList<1> waveList(false);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.02));
      TEST_ASSERT(block.ns() == 101);
      TEST_ASSERT(mesh.dimensions()[0] == 32);

   }

   void testSetup2D()
   {
      printMethod(TEST_FUNC);

      //Create and initialize block
      Block<2> block;
      setupBlock<2>(block);

      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      WaveList<2> waveList(false);
      waveList.allocate(mesh, unitCell);

      block.associate(mesh, fft, unitCell, waveList);

      double ds = 0.26;
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.25));
      TEST_ASSERT(block.ns() == 9);
      TEST_ASSERT(mesh.dimensions()[0] == 32);
      TEST_ASSERT(mesh.dimensions()[1] == 32);
   }

   void testSetup3D()
   {
      printMethod(TEST_FUNC);

      //Create and initialize block
      Block<3> block;
      setupBlock<3>(block);

      Mesh<3> mesh;
      setupMesh<3>(mesh);
      FFT<3> fft;
      fft.setup(mesh.dimensions());

      // Set up unit cell
      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Hexagonal");

      WaveList<3> waveList(false);
      waveList.allocate(mesh, unitCell);

      // Associate block
      block.associate(mesh, fft, unitCell, waveList);

      // Allocate block
      double ds = 0.3;
      block.allocate(ds);

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(block.ns() == 7);
      TEST_ASSERT(eq(block.ds(), 1.0/3.0));
      TEST_ASSERT(mesh.dimensions()[0] == 32);
      TEST_ASSERT(mesh.dimensions()[1] == 32);
      TEST_ASSERT(mesh.dimensions()[2] == 32);
   }

   void testSetupSolver1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      bool isRealField = false;
      WaveList<1> waveList(isRealField);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 4.0));

      // Setup chemical potential field
      CField<1> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i][0] = 1.0;
         w[i][1] = 0.0;
      }

      block.clearUnitCellData();
      block.setupSolver(w);
   }
   
   void testSetupSolver2D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      bool isRealField = false;
      WaveList<2> waveList(isRealField);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Setup chemical potential field
      CField<2> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i][0] = 0.0;
         w[i][1] = 1.0;
      }

      block.clearUnitCellData();
      block.setupSolver(w);
   }

   void testSetupSolver2D_bead()
   {
      printMethod(TEST_FUNC);

      PolymerModel::setModel(PolymerModel::Bead);

      // Create and initialize block
      Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);
      FFT<2> fft;
      fft.setup(mesh.dimensions());

      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      bool isRealField = false;
      WaveList<2> waveList(isRealField);
      waveList.allocate(mesh, unitCell);

      double ds = 0.5;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Setup chemical potential field
      CField<2> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i][0] = 1.5;
         w[i][0] = 0.5;
      }

      block.clearUnitCellData();
      block.setupSolver(w);
   }

   void testSetupSolver3D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<3> block;
      setupBlock<3>(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh<3>(mesh);

      FFT<3> fft;
      fft.setup(mesh.dimensions());

      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Orthorhombic");

      bool isRealField = false;
      WaveList<3> waveList(isRealField);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Setup chemical potential field
      CField<3> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i][0] = 1.5;
         w[i][0] = 0.5;
      }

      block.clearUnitCellData();
      block.setupSolver(w);
   }

   void testSetupSolver3D_domain()
   {
      printMethod(TEST_FUNC);


      Domain<3> domain;
      int nMonomer = readHeader("in/Cubic.hdr", domain);
      TEST_ASSERT(nMonomer == 2);
      Mesh<3>& mesh = domain.mesh();
      UnitCell<3>& unitCell = domain.unitCell();
      FFT<3>& fft = domain.fft();
      WaveList<3>& waveList = domain.waveList();

      // Create and initialize block
      Block<3> block;
      block.associate(mesh, fft, unitCell, waveList);
      setupBlock<3>(block);

      double ds = 0.02;
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 1.8));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 1.8));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 1.8));

      // Setup chemical potential field
      CField<3> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i][0] = 1.5;
         w[i][0] = 0.5;
      }

      block.clearUnitCellData();
      block.setupSolver(w);
   }

   void testSolver1D()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(PolymerModel::isThread());

      // Create and initialize block
      Block<1> block;
      setupBlock<1>(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      bool isRealField = false;
      WaveList<1> waveList(isRealField);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      // Setup chemical potential field
      CField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wcr = 0.3;
      double wci = 0.2;
      fftw_complex wc;
      assign(wc, wcr, wci);
      for (int i=0; i < nx; ++i) {
         w[i][0] = wcr;
         w[i][1] = wci;
      }

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<1>::FieldT qin;
      Propagator<1>::FieldT qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         qin[i][0] = cos(twoPi*double(i)/double(nx));
         qin[i][1] = 0.5 * cos(twoPi*double(i)/double(nx));
         //qin[i][0] = 1.0;
         //qin[i][1] = 0.0;
      }

      block.stepThread(qin, qout);

      double a = unitCell.parameter(0);
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      ds = block.ds();
      double ratioQ = exp(- r * ds);
      fftw_complex ratioW;
      fftw_complex wArg;
      mul(wArg, wc, -1.0*ds);
      assignExp(ratioW, wArg);
      fftw_complex ratio;
      mul(ratio, ratioW, ratioQ);
      fftw_complex expected;
      //setVerbose(1);
      for (int i = 0; i < nx; ++i) {
         mul(expected, qin[i], ratio);
         if (verbose() > 0) {
            Log::file() << "\n" << Dbl(qin[i][0])
                        << "  " << Dbl(qin[i][1])
                        << "  " << Dbl(qout[i][0])
                        << "  " << Dbl(qout[i][1])
                        << "  " << Dbl(expected[0])
                        << "  " << Dbl(expected[1]);
         }
         TEST_ASSERT( eq( qout[i][0], expected[0]) );
         TEST_ASSERT( eq( qout[i][1], expected[1]) );
      }
    
      // Test propagator solve , homogeneous initial condition
      block.propagator(0).solve();

      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(block.propagator(0).head()[i][0], 1.0));
         TEST_ASSERT(eq(block.propagator(0).head()[i][1], 0.0));
      }
     
      ds = block.length();
      ratioQ = exp(- r * ds);
      mul(wArg, wc, -1.0*ds);
      assignExp(ratioW, wArg);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(block.propagator(0).tail()[i][0], ratioW[0]));
         TEST_ASSERT(eq(block.propagator(0).tail()[i][1], ratioW[1]));
      }

   }

   void testSolver1D_domain()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(PolymerModel::isThread());

      Block<1> block;

      Domain<1> domain;
      Mesh<1>& mesh = domain.mesh();
      UnitCell<1>& unitCell = domain.unitCell();
      FFT<1>& fft = domain.fft();
      WaveList<1>& waveList = domain.waveList();
      block.associate(mesh, fft, unitCell, waveList);
      int nMonomer = readHeader("in/Lamellar.hdr", domain);
      TEST_ASSERT(nMonomer == 2);

      // Initialize block
      setupBlock<1>(block);
      double ds = 0.02;
      block.allocate(ds);

      // Setup chemical potential field
      CField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wcr = 0.3;
      double wci = 0.2;
      fftw_complex wc;
      assign(wc, wcr, wci);
      for (int i=0; i < nx; ++i) {
         w[i][0] = wcr;
         w[i][1] = wci;
      }

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<1>::FieldT qin;
      Propagator<1>::FieldT qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         qin[i][0] = cos(twoPi*double(i)/double(nx));
         qin[i][1] = 0.5 * cos(twoPi*double(i)/double(nx));
         //qin[i][0] = 1.0;
         //qin[i][1] = 0.0;
      }

      block.stepThread(qin, qout);

      double a = unitCell.parameter(0);
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      ds = block.ds();
      double ratioQ = exp(- r * ds);
      fftw_complex ratioW;
      fftw_complex wArg;
      mul(wArg, wc, -1.0*ds);
      assignExp(ratioW, wArg);
      fftw_complex ratio;
      mul(ratio, ratioW, ratioQ);
      fftw_complex expected;
      //setVerbose(1);
      for (int i = 0; i < nx; ++i) {
         mul(expected, qin[i], ratio);
         if (verbose() > 0) {
            Log::file() << "\n" << Dbl(qin[i][0])
                        << "  " << Dbl(qin[i][1])
                        << "  " << Dbl(qout[i][0])
                        << "  " << Dbl(qout[i][1])
                        << "  " << Dbl(expected[0])
                        << "  " << Dbl(expected[1]);
         }
         TEST_ASSERT( eq( qout[i][0], expected[0]) );
         TEST_ASSERT( eq( qout[i][1], expected[1]) );
      }
    
      // Test propagator solve , homogeneous initial condition
      block.propagator(0).solve();

      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(block.propagator(0).head()[i][0], 1.0));
         TEST_ASSERT(eq(block.propagator(0).head()[i][1], 0.0));
      }
     
      ds = block.length();
      ratioQ = exp(- r * ds);
      mul(wArg, wc, -1.0*ds);
      assignExp(ratioW, wArg);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(block.propagator(0).tail()[i][0], ratioW[0]));
         TEST_ASSERT(eq(block.propagator(0).tail()[i][1], ratioW[1]));
      }

   }

   #if 0
   void testSolver1D_bead()
   {
      printMethod(TEST_FUNC);

      PolymerModel::setModel(PolymerModel::Bead);

      // Create and initialize block
      Block<1> block;
      setupBlock<1>(block);
      int nBead = block.nBead();

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh<1>(mesh);
      FFT<1> fft;
      fft.setup(mesh.dimensions());

      UnitCell<1> unitCell;
      setupUnitCell<1>(unitCell, "in/Lamellar");

      WaveList<1> waveList(false);
      waveList.allocate(mesh, unitCell);

      double ds = 1.00;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);
      TEST_ASSERT(block.ns() == block.nBead() + 2);
      bool isEnd0 = true;
      bool isEnd1 = true;
      block.propagator(0).setEndFlags(isEnd0, isEnd1);
      block.propagator(1).setEndFlags(isEnd1, isEnd0);

      // Setup chemical potential field
      RField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wc = 0.3;
      for (int i=0; i < nx; ++i) {
         w[i] = wc;
      }

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<1>::FieldT qin;
      Propagator<1>::FieldT qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         qin[i] = cos(twoPi*double(i)/double(nx));
      }

      block.stepBead(qin, qout);
      //double a = 4.0;
      double a = unitCell.parameter(0);
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      double expected = exp(-(wc + r));
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(qout[i], qin[i]*expected));
      }
  
      // Test propagator solve, block owns both vertices
      Propagator<1>& p0 = block.propagator(0);
      p0.solve();

      // Check head slice
      expected = 1.0;
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(p0.head()[i], 1.0));
      }

      // Check first bead, j=1 
      expected = exp(-wc);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(p0.q(1)[i], expected));
      }
      
      // Check last bead, j=nBead
      expected = exp(-wc*nBead);
      for (int i = 0; i < nx; ++i) {
         TEST_ASSERT(eq(p0.q(nBead)[i], expected));
      }

      Propagator<1>& p1 = block.propagator(1);
      p1.solve();

   }

   void testSolver2D()
   {

      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<2> block;
      setupBlock<2>(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh<2>(mesh);

      FFT<2> fft;
      fft.setup(mesh.dimensions());

      UnitCell<2> unitCell;
      setupUnitCell<2>(unitCell, "in/Rectangular");

      WaveList<2> waveList(false);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Setup chemical potential field
      RField<2> w;
      w.allocate(mesh.dimensions());
      MeshIterator<2> iter(mesh.dimensions());

      TEST_ASSERT(w.capacity() == mesh.size());
      double wc = 0.3;
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = wc;
      }

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<2>::FieldT qin;
      Propagator<2>::FieldT qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) ) );
      }
      
      block.stepThread(qin, qout);
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

         TEST_ASSERT(eq(qout[iter.rank()], 
                        qin[iter.rank()]*expected));
      }
    
      // Test propagator solve 
      block.propagator(0).solve();

      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).head()[iter.rank()], 1.0));
      }
      
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).tail()[iter.rank()], expected));
      }

   }

   void testSolver3D()
   {

      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<3> block;
      setupBlock<3>(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh<3>(mesh);

      FFT<3> fft;
      fft.setup(mesh.dimensions());

      UnitCell<3> unitCell;
      setupUnitCell<3>(unitCell, "in/Orthorhombic");

      WaveList<3> waveList(false);
      waveList.allocate(mesh, unitCell);

      double ds = 0.02;
      block.associate(mesh, fft, unitCell, waveList);
      block.allocate(ds);

      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Setup chemical potential field
      RField<3> w;
      w.allocate(mesh.dimensions());
      MeshIterator<3> iter(mesh.dimensions());

      TEST_ASSERT(w.capacity() == mesh.size());
      double wc = 0.3;
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = wc;
      }

      block.clearUnitCellData();
      block.setupSolver(w);

      // Test step
      Propagator<3>::FieldT qin;
      Propagator<3>::FieldT qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) + 
                         double(iter.position(2))/double(mesh.dimension(2)) ) );
      }
      
      block.stepThread(qin, qout);
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

         TEST_ASSERT(eq(qout[iter.rank()], 
                        qin[iter.rank()]*expected));
      }
    
      // Test propagator solve 
      block.propagator(0).solve();

      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).head()[iter.rank()], 1.0));
      }
      
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).tail()[iter.rank()], expected));
      }
   }
   #endif

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
TEST_ADD(PropagatorTest, testSetupSolver3D_domain)
TEST_ADD(PropagatorTest, testSolver1D)
TEST_ADD(PropagatorTest, testSolver1D_domain)
TEST_END(PropagatorTest)

#endif
