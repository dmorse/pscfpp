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

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Block<1> block;
   }

   void setup1DBlock(Block<1>& block) 
   {
      block.setId(0);
      double length = 2.0;
      block.setLength(length);
      block.setMonomerId(1);
      double step = sqrt(6.0);
      block.setKuhn(step);
   }

   void setup1DMesh(Mesh<1>& mesh) {
      IntVec<1> d;
      d[0] = 10;
      mesh.setDimensions(d);
   }

   void setup1DUnitCell(UnitCell<1>& unitCell) 
   {
      std::ifstream in;
      openInputFile("in/Lamellar", in);
      in >> unitCell;
      in.close();
   }

   void testSetDiscretization()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setup1DBlock(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setup1DMesh(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      std::cout << std::endl;
      std::cout << "ns   = " << block.ns() << std::endl;
      std::cout << "mesh = " 
                << block.mesh().dimensions() << std::endl;
   }

   void testSetupSolver()
   {
      printMethod(TEST_FUNC);

      // Create and initialize block
      Block<1> block;
      setup1DBlock(block);

      // Create and initialize mesh
      Mesh<1> mesh;
      setup1DMesh(mesh);

      double ds = 0.02;
      block.setDiscretization(ds, mesh);

      UnitCell<1> unitCell;
      setup1DUnitCell(unitCell);
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

      #if 0
      #endif
   }

   #if 0
   void testPlanarSolve1()
   {
      printMethod(TEST_FUNC);

      // Create and initialize Domain
      double xMin = 0.0;
      double xMax = 1.0;
      int nx = 11;
      Domain domain;
      domain.setPlanarParameters(xMin, xMax, nx);
      TEST_ASSERT(eq(domain.volume(), xMax - xMin));

      // Create and initialize block
      Block b;
      b.setId(0);
      double length = 2.0;
      double ds = 0.02;
      double step = sqrt(6.0);
      b.setLength(length);
      b.setMonomerId(1);
      b.setKuhn(step);
      b.setDiscretization(domain, ds);

      // Create W field
      DArray<double> w;
      w.allocate(nx);
      double wc = 0.3;
      for (int i = 0; i < nx; ++i) {
         w[i] = wc;
      }

      // Solve
      b.setupSolver(w);
      b.propagator(0).solve();

      std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << b.propagator(0).head()[i];
      }
      std::cout << "\n";

      std::cout << "\n Tail:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << b.propagator(0).tail()[i];
      }
      std::cout << "\n";
      std::cout << exp(-wc*b.length()) << "\n";
   }
   #endif

};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor)
TEST_ADD(PropagatorTest, testSetDiscretization)
TEST_ADD(PropagatorTest, testSetupSolver)
TEST_END(PropagatorTest)

#endif
