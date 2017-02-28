#ifndef PSSP_PROPAGATOR_TEST_H
#define PSSP_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

//#include <pssp/domain/Domain.h>
#include <pssp/solvers/Block.h>
#include <pssp/solvers/Propagator.h>
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
TEST_END(PropagatorTest)

#endif
