#ifndef LU_SOLVER_TEST_H
#define LU_SOLVER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/LuSolver.h>

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
};

TEST_BEGIN(LuSolverTest)
TEST_ADD(LuSolverTest, testConstructor)
TEST_ADD(LuSolverTest, testDecompose)
TEST_ADD(LuSolverTest, testSolve)
TEST_END(LuSolverTest)

#endif
