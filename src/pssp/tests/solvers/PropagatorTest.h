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

   void setupMesh1D(Mesh<1>& mesh) {
      IntVec<1> d;
      d[0] = 10;
      mesh.setDimensions(d);
   }

   void setupUnitCell1D(UnitCell<1>& unitCell) 
   {
      std::ifstream in;
      openInputFile("in/Lamellar", in);
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

      std::cout << std::endl;
      std::cout << "ns   = " << block.ns() << std::endl;
      std::cout << "mesh = " 
                << block.mesh().dimensions() << std::endl;
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

      // Setup chemical potential field
      RField<1> w;
      w.allocate(mesh.dimensions());
      TEST_ASSERT(w.capacity() == mesh.size());
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = 1.0;
      }

      block.setupSolver(w, unitCell);

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
      std::cout << std::endl;
      std::cout << "unit cell = " << unitCell << std::endl;

      // Setup chemical potential field
      RField<1> w;
      w.allocate(mesh.dimensions());
      int nx = mesh.size();
      TEST_ASSERT(w.capacity() == nx);
      double wc = 0.3;
      for (int i=0; i < w.capacity(); ++i) {
         w[i] = wc;
      }

      block.setupSolver(w, unitCell);

      #if 0
      Propagator<1>::QField qin;
      Propagator<1>::QField qout;
      qin.allocate(mesh.dimensions());
      qout.allocate(mesh.dimensions());
      for (int i=0; i < w.capacity(); ++i) {
         qin[i] = 1.0;
      }
      block.step(qin, qout);
      std::cout << "\n qout:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << qout[i];
      }
      std::cout << "\n";
      std::cout << exp(-wc*block.ds()) << "\n";
      #endif
     
      block.propagator(0).solve();

      std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << block.propagator(0).head()[i];
      }
      std::cout << "\n";

      std::cout << "\n Tail:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << block.propagator(0).tail()[i];
      }
      std::cout << "\n";
      std::cout << exp(-wc*block.length()) << "\n";

      #if 0 
      #endif

   }

};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor1D)
TEST_ADD(PropagatorTest, testSetDiscretization1D)
TEST_ADD(PropagatorTest, testSetupSolver1D)
TEST_ADD(PropagatorTest, testSolver1D)
TEST_END(PropagatorTest)

#endif
