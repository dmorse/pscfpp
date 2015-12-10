#ifndef TRIDIAGONAL_SOLVER_TEST_H
#define TRIDIAGONAL_SOLVER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pfts/TridiagonalSolver.h>

#include <fstream>

using namespace Util;
using namespace Pfts;

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

   void testDecompose()
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

   void testMultiply()
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

   void testSolve()
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
};

TEST_BEGIN(TridiagonalSolverTest)
TEST_ADD(TridiagonalSolverTest, testConstructor)
TEST_ADD(TridiagonalSolverTest, testDecompose)
TEST_ADD(TridiagonalSolverTest, testMultiply)
TEST_ADD(TridiagonalSolverTest, testSolve)
TEST_END(TridiagonalSolverTest)

#endif
