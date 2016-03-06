#ifndef FD1D_PROPAGATOR_TEST_H
#define FD1D_PROPAGATOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fd1d/Domain.h>
#include <fd1d/Block.h>
#include <fd1d/Propagator.h>
#include <util/math/Constants.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Fd1d;

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
      Block block;
   }

   void testSolve1()
   {
      printMethod(TEST_FUNC);
      Block b;
      b.setId(0);
      double length = 2.0;
      double ds = 0.02;
      double step = sqrt(6.0);
      b.setLength(length);
      b.setMonomerId(1);
      b.setKuhn(step);

      double xMin = 0.0;
      double xMax = 1.0;
      int nx = 11;
      Domain domain;
      domain.setParameters(xMin, xMax, nx);
      b.setDiscretization(domain, ds);
      DArray<double> w;
      w.allocate(nx);
      double wc = 0.3;
      for (int i = 0; i < nx; ++i) {
         w[i] = wc;
      }

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

   void testSolve2()
   {
      printMethod(TEST_FUNC);
      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);

      double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setParameters(xMin, xMax, nx);
      b.setDiscretization(domain, ds);

      DArray<double> q, w;
      q.allocate(nx);
      w.allocate(nx);
      double wc = 0.5;
      for (int i = 0; i < nx; ++i) {
         q[i] = cos(2.0*Constants::Pi*double(i)/double(nx-1));
         w[i] = wc;
      }

      b.setupSolver(w);
      b.propagator(0).solve(q);

      std::cout << "\n Head:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " << b.propagator(0).head()[i];
      }
      std::cout << "\n";

      std::cout << "\n Tail:\n";
      for (int i = 0; i < nx; ++i) {
         std::cout << "  " 
                   << b.propagator(0).tail()[i]/b.propagator(0).head()[i];
      }
      std::cout << "\n";

      double dx = (xMax - xMin)/double(nx - 1);
      double k = 2.0*sin(Constants::Pi/double(nx-1))/dx;
      double f = k*k*step*step/6.0 + wc;
      std::cout << exp(-f*length) << "\n";
   }

};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor)
TEST_ADD(PropagatorTest, testSolve1)
TEST_ADD(PropagatorTest, testSolve2)
TEST_END(PropagatorTest)

#endif
