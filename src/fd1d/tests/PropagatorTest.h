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

   void testPlanarSolve2()
   {
      printMethod(TEST_FUNC);

      // Setup Domain
      double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setPlanarParameters(xMin, xMax, nx);
      TEST_ASSERT(eq(domain.volume(), xMax - xMin));

      // Setup Block
      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);
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

   void testCylinderSolve1()
   {
      printMethod(TEST_FUNC);

      // Setup Domain
      //double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setCylinderParameters(xMax, nx);
      double volume = Constants::Pi*xMax*xMax;
      TEST_ASSERT(eq(domain.volume(), volume));

      // Setup Block
      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);
      b.setDiscretization(domain, ds);
      int ns = b.ns();

      // Create example W field and initial q field
      DArray<double> q, w;
      q.allocate(nx);
      w.allocate(nx);
      double wc = 0.5;
      for (int i = 0; i < nx; ++i) {
         q[i] = 1.0;
         w[i] = wc;
      }
      b.setupSolver(w);
      b.propagator(0).solve(q);

      std::cout << "\n";
      double final = exp(-length*wc);
      double value;
      for (int i = 0; i < nx; ++i) {
         value = b.propagator(0).tail()[i];
         // std::cout << "  " << value;
         TEST_ASSERT(eq(value, final));
      }
      // std::cout << "\n";

      int m = ns/2;
      double sum0 = domain.spatialAverage( b.propagator(0).tail() );
      double sum1 = domain.innerProduct( b.propagator(0).q(m),
                                         b.propagator(0).q(ns-1-m) );
      TEST_ASSERT(eq(sum0, sum1));
      //std::cout << "Average m eq 0  " << sum0 << "\n";
      //std::cout << "Average m neq 0 " << sum1 << "\n";
   }

   void testCylinderSolve2()
   {
      printMethod(TEST_FUNC);

      // Setup Domain
      //double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setCylinderParameters(xMax, nx);
      double volume = Constants::Pi*xMax*xMax;
      TEST_ASSERT(eq(domain.volume(), volume));

      // Setup Block
      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);
      b.setDiscretization(domain, ds);
      int ns = b.ns();

      // Create W and initial q fields
      DArray<double> q, w;
      q.allocate(nx);
      w.allocate(nx);
      double wc = 0.5;
      for (int i = 0; i < nx; ++i) {
         q[i] = 1.0;
         w[i] = wc*cos(2.0*Constants::Pi*double(i)/double(nx-1));
      }

      b.setupSolver(w);
      b.propagator(0).solve(q);

      int m = ns/2;
      double sum0 = domain.spatialAverage( b.propagator(0).tail() );
      double sum1 = domain.innerProduct( b.propagator(0).q(m),
                                         b.propagator(0).q(ns-1-m) );
      // std::cout << "Average m eq 0  " << sum0 << "\n";
      // std::cout << "Average m neq 0 " << sum1 << "\n";
      TEST_ASSERT(eq(sum0, sum1));
   }

   void testSphereSolve1()
   {
      printMethod(TEST_FUNC);

      // Setup Domain
      //double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setSphereParameters(xMax, nx);
      double volume = 4.0*Constants::Pi*xMax*xMax*xMax/3.0;
      TEST_ASSERT(eq(domain.volume(), volume));

      // Setup Block
      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);
      b.setDiscretization(domain, ds);
      int ns = b.ns();

      // Setup W and Initial Q
      DArray<double> q, w;
      q.allocate(nx);
      w.allocate(nx);
      double wc = 0.5;
      for (int i = 0; i < nx; ++i) {
         q[i] = 1.0;
         w[i] = wc;
      }

      b.setupSolver(w);
      b.propagator(0).solve(q);

      std::cout << "\n";
      double final = exp(-length*wc);
      double value;
      for (int i = 0; i < nx; ++i) {
         value = b.propagator(0).tail()[i];
         // std::cout << "  " << value;
         TEST_ASSERT(eq(value, final));
      }
      // std::cout << "\n";

      int m = ns/2;
      double sum0 = domain.spatialAverage( b.propagator(0).tail() );
      double sum1 = domain.innerProduct( b.propagator(0).q(m),
                                         b.propagator(0).q(ns-1-m) );
      TEST_ASSERT(eq(sum0, sum1));
      //std::cout << "Average m eq 0  " << sum0 << "\n";
      //std::cout << "Average m neq 0 " << sum1 << "\n";
   }

   void testSphereSolve2()
   {
      printMethod(TEST_FUNC);

      // Setup Domain
      //double xMin = 0.0;
      double xMax = 1.0;
      int nx = 33;
      Domain domain;
      domain.setSphereParameters(xMax, nx);
      double volume = 4.0*Constants::Pi*xMax*xMax*xMax/3.0;
      TEST_ASSERT(eq(domain.volume(), volume));

      Block b;
      double length = 0.5;
      double ds = 0.00005;
      double step = 1.0;
      b.setId(0);
      b.setMonomerId(1);
      b.setLength(length);
      b.setKuhn(step);
      b.setDiscretization(domain, ds);
      int ns = b.ns();

      DArray<double> q, w;
      q.allocate(nx);
      w.allocate(nx);
      double wc = 0.5;
      for (int i = 0; i < nx; ++i) {
         q[i] = 1.0;
         w[i] = wc*cos(2.0*Constants::Pi*double(i)/double(nx-1));
      }

      b.setupSolver(w);
      b.propagator(0).solve(q);

      int m = ns/2;
      double sum0 = domain.spatialAverage( b.propagator(0).tail() );
      double sum1 = domain.innerProduct( b.propagator(0).q(m),
                                         b.propagator(0).q(ns-1-m) );
      // std::cout << "Average m eq 0  " << sum0 << "\n";
      // std::cout << "Average m neq 0 " << sum1 << "\n";
      TEST_ASSERT(eq(sum0, sum1));
   }

};

TEST_BEGIN(PropagatorTest)
TEST_ADD(PropagatorTest, testConstructor)
TEST_ADD(PropagatorTest, testPlanarSolve1)
TEST_ADD(PropagatorTest, testPlanarSolve2)
TEST_ADD(PropagatorTest, testCylinderSolve1)
TEST_ADD(PropagatorTest, testCylinderSolve2)
TEST_ADD(PropagatorTest, testSphereSolve1)
TEST_ADD(PropagatorTest, testSphereSolve2)
TEST_END(PropagatorTest)

#endif
