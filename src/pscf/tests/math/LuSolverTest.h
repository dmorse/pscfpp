#ifndef LU_SOLVER_TEST_H
#define LU_SOLVER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/LuSolver.h>
#include <util/containers/DMatrix.h>

#include <cmath>
#include <fstream>

using namespace Util;
using namespace Pscf;

class LuSolverTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      LuSolver solver;
      solver.allocate(3);
   }

   void testDecompose()
   {
      printMethod(TEST_FUNC);

      DMatrix<double> a;
      a.allocate(3, 3);
      a(0,0) = 1.0;
      a(1,1) = 3.0;
      a(2,2) = 4.0;
      a(0,1) = 2.0;
      a(1,0) = 2.0;
      a(1,2) = 5.0;
      a(2,1) = 5.0;
      a(0,2) = 0.0;
      a(2,0) = 0.0;

      LuSolver solver;
      solver.allocate(3);
      solver.computeLU(a);
   }

   void testSolve()
   {
      printMethod(TEST_FUNC);

      DMatrix<double> a;
      a.allocate(3,3);
      a(0,0) = 1.0;
      a(1,1) = 3.0;
      a(2,2) = 4.0;
      a(0,1) = 2.0;
      a(1,0) = 2.0;
      a(1,2) = 5.0;
      a(2,1) = 5.0;
      a(0,2) = 0.0;
      a(2,0) = 0.0;

      DArray<double> b, x, y;
      b.allocate(3);
      x.allocate(3);
      y.allocate(3);

      b[0] = 1.0;
      b[1] = 2.0;
      b[2] = 3.0;

      LuSolver solver;
      solver.allocate(3);
      solver.computeLU(a);
      solver.solve(b, x);

      int i, j;
      for (i = 0; i < 3; ++i) {
         y[i] = 0.0;
         for (j = 0; j < 3; ++j) {
            y[i] += a(i,j)*x[j];
         }
      }
      // std::cout << "y = " << y[0] << "  " << y[1] << "  " << y[2] << "\n";

      TEST_ASSERT(eq(b[0], y[0]));
      TEST_ASSERT(eq(b[1], y[1]));
      TEST_ASSERT(eq(b[2], y[2]));
   }

   /*
   * Matrices larger than 24 x 24, which gsl_linalg_LU_decomp (GSL 2.6
   * and later) decomposes recursively. For a singular matrix, it did
   * not complete the decomposition (GSL 2.6 and 2.7 fill the remaining
   * submatrix with NaN; GSL 2.8 stops at the zero pivot and leaves the
   * decomposition and the permutation incomplete, or crashes).
   * computeLU must complete the decomposition, so that isSingular()
   * reports the zero pivot and solveBasic() solves a consistent system.
   */
   void testSingularLarge()
   {
      printMethod(TEST_FUNC);

      const int n = 100;
      DMatrix<double> a;
      a.allocate(n, n);
      int i, j;
      for (i = 0; i < n; ++i) {
         for (j = 0; j < n; ++j) {
            a(i,j) = std::sin(1.0 + 7.0*i + 3.0*j*j);
         }
         a(i,i) += 20.0;
      }

      // Nonsingular matrix: solve Ax = b with a known solution x
      DArray<double> b, x;
      b.allocate(n);
      x.allocate(n);
      for (i = 0; i < n; ++i) {
         b[i] = 0.0;
         for (j = 0; j < n; ++j) {
            b[i] += a(i,j)*std::cos(0.5 + 1.3*j);
         }
      }
      LuSolver solver;
      solver.allocate(n);
      solver.computeLU(a);
      TEST_ASSERT(!solver.isSingular());
      solver.solve(b, x);
      for (i = 0; i < n; ++i) {
         TEST_ASSERT(std::abs(x[i] - std::cos(0.5 + 1.3*i)) < 1.0E-12);
      }

      // Singular matrix: columns 0 and n/2 are e_{n/2} and -e_{n/2}
      for (i = 0; i < n; ++i) {
         a(i, 0) = 0.0;
         a(i, n/2) = 0.0;
      }
      a(n/2, 0) = 1.0;
      a(n/2, n/2) = -1.0;
      solver.computeLU(a);
      TEST_ASSERT(solver.isSingular());

      // Consistent b = Ax: the solution is unique up to the null
      // vector, which changes x[0] and x[n/2] by the same amount
      for (i = 0; i < n; ++i) {
         b[i] = 0.0;
         for (j = 0; j < n; ++j) {
            b[i] += a(i,j)*std::cos(0.5 + 1.3*j);
         }
      }
      solver.solveBasic(b, x);
      double r;
      for (i = 0; i < n; ++i) {
         r = -b[i];
         for (j = 0; j < n; ++j) {
            r += a(i,j)*x[j];
         }
         TEST_ASSERT(std::abs(r) < 1.0E-12);
         if (i != 0 && i != n/2) {
            TEST_ASSERT(std::abs(x[i] - std::cos(0.5 + 1.3*i)) < 1.0E-12);
         }
      }
      TEST_ASSERT(std::abs(x[0] - x[n/2] - std::cos(0.5)
                           + std::cos(0.5 + 1.3*(n/2))) < 1.0E-12);
   }

};

TEST_BEGIN(LuSolverTest)
TEST_ADD(LuSolverTest, testConstructor)
TEST_ADD(LuSolverTest, testDecompose)
TEST_ADD(LuSolverTest, testSolve)
TEST_ADD(LuSolverTest, testSingularLarge)
TEST_END(LuSolverTest)

#endif
