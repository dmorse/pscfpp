#ifndef TRIDIAGONAL_SOLVER_TEST_H
#define TRIDIAGONAL_SOLVER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/TridiagonalSolver.h>

#include <cmath>
#include <fstream>

using namespace Util;
using namespace Pscf;

class TridiagonalSolverTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);
   }

   void testSymmetricDecompose()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u;
      d.allocate(3);
      u.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      solver.computeLU(d, u);
   }

   void testSymmetricMultiply()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u;
      d.allocate(3);
      u.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      solver.computeLU(d, u);
 
      DArray<double> b, x;
      b.allocate(3);
      x.allocate(3);
      // std::cout << "\n";
      b[0] = 1.0;
      b[1] = 0.0;
      b[2] = 0.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 1.0));
      TEST_ASSERT(eq(x[1], 2.0));
      TEST_ASSERT(eq(x[2], 0.0));
      b[0] = 0.0;
      b[1] = 1.0;
      b[2] = 0.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 2.0));
      TEST_ASSERT(eq(x[1], 3.0));
      TEST_ASSERT(eq(x[2], 5.0));
      b[0] = 0.0;
      b[1] = 0.0;
      b[2] = 1.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 0.0));
      TEST_ASSERT(eq(x[1], 5.0));
      TEST_ASSERT(eq(x[2], 4.0));
   }

   void testSymmetricSolve()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u;
      d.allocate(3);
      u.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      solver.computeLU(d, u);
 
      DArray<double> b, x, y;
      b.allocate(3);
      x.allocate(3);
      b[0] = 1.0;
      b[1] = 2.0;
      b[2] = 3.0;
      solver.solve(b, x);
      y.allocate(3);
      solver.multiply(x, y);
      // std::cout << y[0] << "  " << y[1] << "  " << y[2] << "\n";
      TEST_ASSERT(eq(b[0], y[0]));
      TEST_ASSERT(eq(b[1], y[1]));
      TEST_ASSERT(eq(b[2], y[2]));
   }

   void testDecompose()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u, l;
      d.allocate(3);
      u.allocate(2);
      l.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      l[0] = 1.0;
      l[1] = -3.0;
      solver.computeLU(d, u, l);
   }

   void testMultiply()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u, l;
      d.allocate(3);
      u.allocate(2);
      l.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      l[0] = 1.0;
      l[1] = -3.0;
      solver.computeLU(d, u, l);
 
      DArray<double> b, x;
      b.allocate(3);
      x.allocate(3);
      // std::cout << "\n";
      b[0] = 1.0;
      b[1] = 0.0;
      b[2] = 0.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 1.0));
      TEST_ASSERT(eq(x[1], 1.0));
      TEST_ASSERT(eq(x[2], 0.0));
      b[0] = 0.0;
      b[1] = 1.0;
      b[2] = 0.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 2.0));
      TEST_ASSERT(eq(x[1], 3.0));
      TEST_ASSERT(eq(x[2], -3.0));
      b[0] = 0.0;
      b[1] = 0.0;
      b[2] = 1.0;
      solver.multiply(b, x);
      // std::cout << x[0] << "  " << x[1] << "  " << x[2] << "\n";
      TEST_ASSERT(eq(x[0], 0.0));
      TEST_ASSERT(eq(x[1], 5.0));
      TEST_ASSERT(eq(x[2], 4.0));
   }

   void testSolve()
   {
      printMethod(TEST_FUNC);
      TridiagonalSolver solver;
      solver.allocate(3);

      DArray<double> d, u, l;
      d.allocate(3);
      u.allocate(2);
      l.allocate(2);
      d[0] = 1.0;
      d[1] = 3.0;
      d[2] = 4.0;
      u[0] = 2.0;
      u[1] = 5.0;
      l[0] = 1.0;
      l[1] = -3.0;
      solver.computeLU(d, u, l);
 
      DArray<double> b, x, y;
      b.allocate(3);
      x.allocate(3);
      b[0] = 1.0;
      b[1] = 2.0;
      b[2] = 3.0;
      solver.solve(b, x);
      y.allocate(3);
      solver.multiply(x, y);
      // std::cout << y[0] << "  " << y[1] << "  " << y[2] << "\n";
      TEST_ASSERT(eq(b[0], y[0]));
      TEST_ASSERT(eq(b[1], y[1]));
      TEST_ASSERT(eq(b[2], y[2]));
   }

   /*
   * Matrices with a zero or small pivot in Gauss elimination without
   * row swaps: solve() gave NaN or an inaccurate solution.
   */
   void testSolvePivot()
   {
      printMethod(TEST_FUNC);

      double d1[3] = {0.0, 2.0, 2.0};
      double u1[2] = {1.0, 1.0};
      TEST_ASSERT(solvePivot(3, d1, u1, u1));

      double d2[3] = {1.0E-20, 2.0, 2.0};
      TEST_ASSERT(solvePivot(3, d2, u1, u1));

      // Swaps of consecutive rows
      double d3[6] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
      double u3[5] = {0.125, 0.125, 0.125, 0.125, 0.125};
      double l3[5] = {2.0, 2.0, 2.0, 2.0, 2.0};
      TEST_ASSERT(solvePivot(6, d3, u3, l3));

      double d4[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
      double u4[5] = {2.0, 2.0, 2.0, 2.0, 2.0};
      double l4[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
      TEST_ASSERT(solvePivot(6, d4, u4, l4));
   }

   #if 0
   #endif

private:

   /*
   * Solve Ax = b for the tridiagonal matrix A with diagonal d, upper
   * off-diagonal u and lower off-diagonal l, and b = A x0 with x0[i] =
   * i + 1. Return true if solve() gives x0 and multiply() gives b.
   */
   bool solvePivot(int n, double const * d, double const * u,
                   double const * l)
   {
      DArray<double> dA, uA, lA, x0, b, x;
      dA.allocate(n);
      uA.allocate(n-1);
      lA.allocate(n-1);
      x0.allocate(n);
      b.allocate(n);
      x.allocate(n);
      int i;
      for (i = 0; i < n; ++i) {
         dA[i] = d[i];
         x0[i] = double(i + 1);
      }
      for (i = 0; i < n - 1; ++i) {
         uA[i] = u[i];
         lA[i] = l[i];
      }
      for (i = 0; i < n; ++i) {
         b[i] = d[i]*x0[i];
         if (i > 0) b[i] += l[i-1]*x0[i-1];
         if (i < n - 1) b[i] += u[i]*x0[i+1];
      }

      TridiagonalSolver solver;
      solver.allocate(n);
      solver.computeLU(dA, uA, lA);
      solver.solve(b, x);
      bool success = true;
      for (i = 0; i < n; ++i) {
         if (!(std::abs(x[i] - x0[i]) < 1.0E-12*x0[i])) success = false;
      }
      solver.multiply(x0, x);
      for (i = 0; i < n; ++i) {
         if (!(std::abs(x[i] - b[i]) < 1.0E-12)) success = false;
      }
      return success;
   }

};

TEST_BEGIN(TridiagonalSolverTest)
TEST_ADD(TridiagonalSolverTest, testConstructor)
TEST_ADD(TridiagonalSolverTest, testSymmetricDecompose)
TEST_ADD(TridiagonalSolverTest, testSymmetricMultiply)
TEST_ADD(TridiagonalSolverTest, testSymmetricSolve)
TEST_ADD(TridiagonalSolverTest, testDecompose)
TEST_ADD(TridiagonalSolverTest, testMultiply)
TEST_ADD(TridiagonalSolverTest, testSolve)
TEST_ADD(TridiagonalSolverTest, testSolvePivot)
TEST_END(TridiagonalSolverTest)

#endif
