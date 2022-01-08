#ifndef PSPC_PROPAGATOR_TEST_H
#define PSPC_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/solvers/Block.h>
#include <pscf/mesh/MeshIterator.h>
#include <pspc/solvers/Propagator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class PropagatorTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      Block<1> block;
   }

   void setupBlock1D(Block<1>& block) 
   {
      block.setId(0);
      double length = 2.0;
      block.setLength(length);
      block.setMonomerId(1);
      double step = sqrt(6.0);
      block.setKuhn(step);
   }

   void setupBlock2D(Block<2>& block)
   {
      block.setId(0);
      double length = 2.0;
      block.setLength(length);
      block.setMonomerId(1);
      double step = sqrt(6.0);
      block.setKuhn(step);
   }

   void setupBlock3D(Block<3>& block) 
   {
      block.setId(0);
      double length = 2.0;
      block.setLength(length);
      block.setMonomerId(1);
      double step = sqrt(6.0);
      //double step = 0;
      block.setKuhn(step);
   }

   void setupMesh1D(Mesh<1>& mesh) {
      IntVec<1> d;
      d[0] = 10;
      mesh.setDimensions(d);
   }

   void setupMesh2D(Mesh<2>& mesh) {
      IntVec<2> d;
      d[0] = 10;
      d[1] = 10;
      mesh.setDimensions(d);
   }

   void setupMesh3D(Mesh<3>& mesh) {
      IntVec<3> d;
      d[0] = 10;
      d[1] = 10;
      d[2] = 10;
      mesh.setDimensions(d);
   }

   void setupUnitCell1D(UnitCell<1>& unitCell) 
   {
      std::ifstream in;
      openInputFile("in/Lamellar", in);
      in >> unitCell;
      in.close();
   }

   void setupUnitCell2D(UnitCell<2>& unitCell)
   {
      std::ifstream in;
      openInputFile("in/Rectangular", in);
      in >> unitCell;
      in.close();
   }

   void setupUnitCell3D(UnitCell<3>& unitCell) 
   {
      std::ifstream in;
      openInputFile("in/Orthorhombic", in);
      in >> unitCell;
      in.close();
   }

   void testSetDiscretization1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setupBlock1D(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh1D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);
      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.02));
      TEST_ASSERT(block.ns() == 101);
      TEST_ASSERT(block.mesh().dimensions()[0] == 10);

      // std::cout << std::endl;
      // std::cout << "ns   = " << block.ns() << std::endl;
      // std::cout << "mesh = " 
      //          << block.mesh().dimensions() << std::endl;
   }

   void testSetDiscretization2D()
   {
      printMethod(TEST_FUNC);

      //Create and initialize block
      Block<2> block;
      setupBlock2D(block);

      Mesh<2> mesh;
      setupMesh2D(mesh);

      double ds = 0.26;
      block.setDiscretization(ds, mesh);
      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.25));
      TEST_ASSERT(block.ns() == 9);
      TEST_ASSERT(block.mesh().dimensions()[0] == 10);
      TEST_ASSERT(block.mesh().dimensions()[1] == 10);
   }

   void testSetDiscretization3D()
   {
      printMethod(TEST_FUNC);

      //Create and initialize block
      Block<3> block;
      setupBlock3D(block);

      Mesh<3> mesh;
      setupMesh3D(mesh);

      double ds = 0.3;
      block.setDiscretization(ds, mesh);
      // std::cout << "block len = " << block.length() << "\n";
      // std::cout << "block ns  = " << block.ns() << "\n";
      // std::cout << "block ds  = " << block.ds() << "\n";

      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(block.ns() == 7);
      TEST_ASSERT(eq(block.ds(), 1.0/3.0));
      TEST_ASSERT(block.mesh().dimensions()[0] == 10);
      TEST_ASSERT(block.mesh().dimensions()[1] == 10);
      TEST_ASSERT(block.mesh().dimensions()[2] == 10);
   }

   void testSetupSolver1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setupBlock1D(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh1D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<1> unitCell;
      setupUnitCell1D(unitCell);
      // // std::cout << std::endl;
      // // std::cout << "unit cell = " << unitCell << std::endl;
      // TEST_ASSERT(eq(unitCell.rBasis(0)[0], 4.0));

      // // Setup chemical potential field
      // RField<1> w;
      // w.allocate(mesh.dimensions());
      // TEST_ASSERT(w.capacity() == mesh.size());
      // for (int i=0; i < w.capacity(); ++i) {
      //    w[i] = 1.0;
      // }

      // block.setupUnitCell(unitCell);
      // block.setupSolver(w);
   }
   
   void testSetupSolver2D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<2> block;
      setupBlock2D(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh2D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<2> unitCell;
      setupUnitCell2D(unitCell);
      // std::cout << std::endl;
      // std::cout << "unit cell = " << unitCell << std::endl;
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));

      // Setup chemical potential field
      RField<2> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = 1.0;
      }

      block.setupUnitCell(unitCell);
      block.setupSolver(w);
   }

   void testSetupSolver3D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<3> block;
      setupBlock3D(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh3D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<3> unitCell;
      setupUnitCell3D(unitCell);
      // std::cout << std::endl;
      // std::cout << "unit cell = " << unitCell << std::endl;
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Setup chemical potential field
      RField<3> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = 1.0;
      }

      block.setupUnitCell(unitCell);
      block.setupSolver(w);
   }

   void testSolver1D()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setupBlock1D(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setupMesh1D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<1> unitCell;
      setupUnitCell1D(unitCell);
      // std::cout << std::endl;
      // std::cout << "unit cell = " << unitCell << std::endl;

      // Setup chemical potential field
      RField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wc = 0.3;
      for (int i=0; i < nx; ++i) {
         w[i] = wc;
      }

      block.setupUnitCell(unitCell);
      block.setupSolver(w);

      // Test step
      Propagator<1>::QField qin;
      Propagator<1>::QField qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         qin[i] = cos(twoPi*double(i)/double(nx));
      }
      block.step(qin, qout);
      double a = 4.0;
      double b = block.kuhn();
      double Gb = twoPi*b/a;
      double r = Gb*Gb/6.0;
      ds = block.ds();
      double expected = exp(-(wc + r)*ds);
      for (int i = 0; i < nx; ++i) {
         // std::cout << "  " << qout[i]
         //           << "  " << qin[i]
         //           << "  " << qin[i]*expected
         //           << "  " << expected << std::endl;
         TEST_ASSERT(eq(qout[i], qin[i]*expected));
      }
      // std::cout << "\n";
      // std::cout << "expected ratio = " << expected << "\n";
    
      // Test propagator solve 
      block.propagator(0).solve();

      // std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         // std::cout << "  " << block.propagator(0).head()[i];
         TEST_ASSERT(eq(block.propagator(0).head()[i],1.0));
      }
      // std::cout << "\n";
      

      // std::cout << "\n Tail:\n";
      expected = exp(-wc*block.length());
      for (int i = 0; i < nx; ++i) {
         // std::cout << "  " << block.propagator(0).tail()[i];
         TEST_ASSERT(eq(block.propagator(0).tail()[i], expected));
      }
      // std::cout << "\n";
      // std::cout << exp(-wc*block.length()) << "\n";

   }

   void testSolver2D()
   {

      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<2> block;
      setupBlock2D(block);

      // Create and initialize mesh
      Mesh<2> mesh;
      setupMesh2D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<2> unitCell;
      setupUnitCell2D(unitCell);
      // std::cout << std::endl;
      // std::cout << "unit cell = " << unitCell << std::endl;
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

      block.setupUnitCell(unitCell);
      block.setupSolver(w);

      // Test step
      Propagator<2>::QField qin;
      Propagator<2>::QField qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) ) );
      }
      
      // std::cout<<std::endl;
      block.step(qin, qout);
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

         // std::cout << "  " << qout[iter.rank()]
         //          << "  " << qin[iter.rank()]
         //         << "  " << qin[iter.rank()]*expected
         //          << "  " << expected << std::endl;
         //TEST_ASSERT(eq(qout[iter.rank()], 
         //               qin[iter.rank()]*expected));
      }

      // std::cout << "\n";
      // std::cout << "expected ratio = " << expected << "\n";
    
      // Test propagator solve 
      block.propagator(0).solve();

      // std::cout << "\n Head:\n";
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).head()[iter.rank()], 1.0));
      }
      // std::cout << "\n";
      
      // std::cout << "\n Tail:\n";
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).tail()[iter.rank()], expected));
      }
      // std::cout << "\n";
      // std::cout << exp(-wc*block.length()) << "\n";

   }

   void testSolver3D()
   {

      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<3> block;
      setupBlock3D(block);

      // Create and initialize mesh
      Mesh<3> mesh;
      setupMesh3D(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<3> unitCell;
      setupUnitCell3D(unitCell);
      // std::cout << std::endl;
      // std::cout << "unit cell = " << unitCell << std::endl;
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

      block.setupUnitCell(unitCell);
      block.setupSolver(w);

      // Test step
      Propagator<3>::QField qin;
      Propagator<3>::QField qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());

      double twoPi = 2.0*Constants::Pi;
      for (iter.begin(); !iter.atEnd(); ++iter){
         qin[iter.rank()] = cos(twoPi * 
                        (double(iter.position(0))/double(mesh.dimension(0)) + 
                         double(iter.position(1))/double(mesh.dimension(1)) + 
                         double(iter.position(2))/double(mesh.dimension(2)) ) );
      }
      
      // std::cout<<std::endl;
      block.step(qin, qout);
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

         // std::cout << "  " << qout[iter.rank()]
         //          << "  " << qin[iter.rank()]
         //         << "  " << qin[iter.rank()]*expected
         //          << "  " << expected << std::endl;
         // TEST_ASSERT(eq(qout[iter.rank()], 
         //               qin[iter.rank()]*expected));
      }

      // std::cout << "\n";
      // std::cout << "expected ratio = " << expected << "\n";
    
      // Test propagator solve 
      block.propagator(0).solve();

      // std::cout << "\n Head:\n";
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).head()[iter.rank()], 1.0));
      }
      // std::cout << "\n";
      
      // std::cout << "\n Tail:\n";
      expected = exp(-wc*block.length());
      for (iter.begin(); !iter.atEnd(); ++iter){
         TEST_ASSERT(eq(block.propagator(0).tail()[iter.rank()], expected));
      }
      // std::cout << "\n";
      // std::cout << exp(-wc*block.length()) << "\n";

   }

};

TEST_BEGIN(PropagatorTest)
// TEST_ADD(PropagatorTest, testConstructor1D)
// TEST_ADD(PropagatorTest, testSetDiscretization1D)
// TEST_ADD(PropagatorTest, testSetDiscretization2D)
// TEST_ADD(PropagatorTest, testSetDiscretization3D)
TEST_ADD(PropagatorTest, testSetupSolver1D)
// TEST_ADD(PropagatorTest, testSetupSolver3D)
// TEST_ADD(PropagatorTest, testSolver1D)
// TEST_ADD(PropagatorTest, testSolver2D)
// TEST_ADD(PropagatorTest, testSolver3D)
TEST_END(PropagatorTest)

#endif
