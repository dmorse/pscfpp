#ifndef PSSP_PROPAGATOR_TEST_H
#define PSSP_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/solvers/Block.h>
#include <pssp/solvers/Propagator.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/math/IntVec.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

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


   void setupUnitCell3D(UnitCell<3>& unitCell) 
   {
      std::ifstream in;
      openInputFile("in/Cubic", in);
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

      //std::cout << std::endl;
      //std::cout << "ns   = " << block.ns() << std::endl;
      //std::cout << "mesh = " 
      //          << block.mesh().dimensions() << std::endl;
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
      TEST_ASSERT(eq(block.length(), 2.0));
      TEST_ASSERT(eq(block.ds(), 0.25));
      TEST_ASSERT(block.ns() == 9);
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
      std::cout << std::endl;
      std::cout << "unit cell = " << unitCell << std::endl;
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 4.0));

      // Setup chemical potential field
      RField<1> w;
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
      std::cout << std::endl;
      std::cout << "unit cell = " << unitCell << std::endl;
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
         //std::cout << "  " << qout[i]
         //          << "  " << qin[i]
         //          << "  " << qin[i]*expected
         //          << "  " << expected << std::endl;
         TEST_ASSERT(eq(qout[i], qin[i]*expected));
      }
      //std::cout << "\n";
      //std::cout << "expected ratio = " << expected << "\n";
    
      // Test propagator solve 
      block.propagator(0).solve();

      //std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         //std::cout << "  " << block.propagator(0).head()[i];
         TEST_ASSERT(eq(block.propagator(0).head()[i],1.0));
      }
      std::cout << "\n";
      

      // std::cout << "\n Tail:\n";
      expected = exp(-wc*block.length());
      for (int i = 0; i < nx; ++i) {
         //std::cout << "  " << block.propagator(0).tail()[i];
         TEST_ASSERT(eq(block.propagator(0).tail()[i], expected));
      }
      // std::cout << "\n";
      //std::cout << exp(-wc*block.length()) << "\n";

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
      //std::cout << std::endl;
      //std::cout << "unit cell = " << unitCell << std::endl;
      TEST_ASSERT(eq(unitCell.rBasis(0)[0], 3.0));
      TEST_ASSERT(eq(unitCell.rBasis(1)[1], 4.0));
      TEST_ASSERT(eq(unitCell.rBasis(2)[2], 5.0));

      // Setup chemical potential field
      RField<3> w;
      w.allocate(mesh.dimensions());
      IntVec<3> dim = mesh.dimensions();
      int nx = dim[0];
      int ny = dim[1];
      int nz = dim[2];
      

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
      for (int i=0; i < nx; ++i) {
         for (int j=0; j < ny; ++j){
            for (int k=0; k < nz; ++k){
               dim[0] = i;
               dim[1] = j;
               dim[2] = k;
               qin[mesh.rank(dim)] = cos(twoPi*(double(i)/double(nx)
                                     + double(j)/double(ny)
                                     + double(k)/double(nz)));
               //qin[mesh.rank(dim)] = 1;
            }
         }
      }
      
      std::cout<<std::endl;
      block.step(qin, qout);
      double b = block.kuhn();
      double Gb;
      double expected;
      IntVec<3> temp;
      ds = block.ds();
      for (int i=0; i < nx; ++i) {
         for (int j=0; j < ny; ++j){
            for (int k=0; k < nz; ++k){
               dim[0] = i;
               dim[1] = j;
               dim[2] = k;
               temp[0] = 1;
               temp[1] = 1;
               temp[2] = 1;
               Gb = unitCell.ksq(temp);
              // Gb = 1;
               //double factor = twoPi*b/4.0;
               double factor = b;
               double r = Gb*factor*factor/6.0;
               expected = exp(-(wc + r)*ds);

               //std::cout << "  " << qout[mesh.rank(dim)]
               //          << "  " << qin[mesh.rank(dim)]
               //         << "  " << qin[mesh.rank(dim)]*expected
               //          << "  " << expected << std::endl;
               TEST_ASSERT(eq(qout[mesh.rank(dim)], 
                            qin[mesh.rank(dim)]*expected));
            }
         }
      }

      //std::cout << "\n";
      //std::cout << "expected ratio = " << expected << "\n";
    
      // Test propagator solve 
      block.propagator(0).solve();

      //std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k){
               //std::cout << "  " << block.propagator(0).head()[i];
               dim[0] = i;
               dim[1] = j;
               dim[2] = k;
               TEST_ASSERT(eq(block.propagator(0).head()[mesh.rank(dim)],1.0)); 
            }
         }
      }
      std::cout << "\n";
      

      // std::cout << "\n Tail:\n";
      expected = exp(-wc*block.length());
      for (int i = 0; i < nx; ++i) {
         for (int j = 0; j < ny; ++j){
            for (int k = 0; k < nz; ++k){
               dim[0] = i;
               dim[1] = j;
               dim[2] = k;
               //std::cout << "  " << block.propagator(0).tail()[i];
               TEST_ASSERT(eq(block.propagator(0).tail()[i], expected));
            }
         }
      }
      // std::cout << "\n";
      //std::cout << exp(-wc*block.length()) << "\n";

   }

};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor1D)
TEST_ADD(PropagatorTest, testSetDiscretization1D)
TEST_ADD(PropagatorTest, testSetDiscretization3D)
TEST_ADD(PropagatorTest, testSetupSolver1D)
TEST_ADD(PropagatorTest, testSetupSolver3D)
TEST_ADD(PropagatorTest, testSolver1D)
TEST_ADD(PropagatorTest, testSolver3D)
TEST_END(PropagatorTest)

#endif
