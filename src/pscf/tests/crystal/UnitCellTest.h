#ifndef PSCF_UNIT_CELL_TEST_H
#define PSCF_UNIT_CELL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <util/math/Constants.h>
#include <util/format/Int.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class UnitCellTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   template <int D>
   bool isValidReciprocal(UnitCell<D> cell)
   {
      double sum;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (i=0; i < D; ++i ) {
         for (j=0; j < D; ++j ) {
            sum = 0.0;
            for (k=0; k < D; ++k ) {
               sum += cell.rBasis(i)[k]*cell.kBasis(j)[k];  
            }
            if (i == j) {
               sum -= twoPi;
            }
            if (std::abs(sum) > 1.0E-8) {
               return false;
            }
         }
      }
      return true;
   }

   template <int D>
   bool isValidDerivative(UnitCell<D> cell)
   {
      double sum;
      double nParams = cell.nParams();
      int i, j, k, m;
      for (k = 0; k < nParams; ++k) {
         for (i = 0; i < D; ++i) {
            for (j = 0; j < D; ++j) {
               sum = 0.0;
               for (m = 0; m < D; ++m) {
                  sum += cell.drBasis(k, i, m)*cell.kBasis(j)[m];
                  sum += cell.dkBasis(k, j, m)*cell.rBasis(i)[m];
               }
            }
            if (std::abs(sum) > 1.0E-8) {
               return false;
            }
         }
      }
      return true;
   }

   void test1DLamellar() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<1> v;
      std::ifstream in;
      openInputFile("in/Lamellar", in);
      in >> v;
      double param = v.params()[0];
      double twoPi = 2.0*Constants::Pi;

      TEST_ASSERT(eq(v.rBasis(0)[0], param));
      TEST_ASSERT(eq(v.kBasis(0)[0], twoPi/param));
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      // Test ksq function
      IntVec<1> x;
      int m = -3;
      x[0] = m;
      double xSq = v.ksq(x);
      double y = (twoPi*m)/param;
      TEST_ASSERT(eq(xSq, y*y));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;
      std::cout << v.rBasis(0) << std::endl;
      std::cout << v.kBasis(0) << std::endl;
      #endif

   }



   void test2DSquare() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Square", in);

      in >> v;
      TEST_ASSERT(v.nParams() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      #endif

      double param = v.params()[0];
      double twoPi = 2.0*Constants::Pi;
      double b, dbb;
      int i, j, k;
      for (k = 0; k < v.nParams(); ++k) {
         for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
               if (i == j) {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
               if (i == j) {
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }
   }

   void test2DHexagonal() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Hexagonal", in);

      in >> v;
      TEST_ASSERT(v.nParams() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      #endif

      IntVec<2> d;
      d[0] = 8;
      d[1] = 8;
      IntVec<2> x;
      x[0] = -4;
      x[1] = +5;
      IntVec<2> y;
     
      std::cout << "Before shift " << x << std::endl;
      y = shiftToMinimum(x, d, v);
      std::cout << "After shift  " << y << std::endl;
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == -3);
      y = shiftToMinimum(y, d, v);
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == -3);
      //std::cout << "After again  " << y << std::endl;

   }

   void test3DOrthorhombic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Orthorhombic", in);
      in >> v;

      TEST_ASSERT(v.nParams() == 3);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "a(2) = " << v.rBasis(2) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      std::cout << "b(2) = " << v.kBasis(2) << std::endl;
      #endif

      double param, b, dbb;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (k = 0; k < v.nParams(); ++k) {
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.params()[i];
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.params()[i];
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }

   }

   void test3DCubic() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      UnitCell<3> v;
      std::ifstream in;
      openInputFile("in/Cubic", in);
      in >> v;

      TEST_ASSERT(v.nParams() == 1);
      TEST_ASSERT(isValidReciprocal(v));
      TEST_ASSERT(isValidDerivative(v));

      #if 0
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      std::cout << "a(0) = " << v.rBasis(0) << std::endl;
      std::cout << "a(1) = " << v.rBasis(1) << std::endl;
      std::cout << "a(2) = " << v.rBasis(2) << std::endl;
      std::cout << "b(0) = " << v.kBasis(0) << std::endl;
      std::cout << "b(1) = " << v.kBasis(1) << std::endl;
      std::cout << "b(2) = " << v.kBasis(2) << std::endl;
      #endif

      #if 0
      double param, b, dbb;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (k = 0; k < v.nParams(); ++k) {
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.params()[i];
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 2.0*param));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
         for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
               if (i == j && i == k) {
                  param = v.params()[i];
                  b = twoPi/param;
                  dbb = -2.0*b*b/param;
                  TEST_ASSERT(eq(v.dkkBasis(k, i, j), dbb));
               } else {
                  TEST_ASSERT(eq(v.drrBasis(k, i, j), 0.0));
               }
            }
         }
      }
      #endif

      IntVec<3> d;
      d[0] = 8;
      d[1] = 8;
      d[2] = 8;
      IntVec<3> x;
      x[0] = -4;
      x[1] = +4;
      x[2] =  7;
      IntVec<3> y;
     
      //std::cout << "Before shift " << x << std::endl;
      y = shiftToMinimum(x, d, v);
      TEST_ASSERT(y[0] == 4);
      TEST_ASSERT(y[1] == 4);
      TEST_ASSERT(y[2] == -1);
      //std::cout << "After shift  " << y << std::endl;
      //y = shiftToMinimum(y, d, v);
      //std::cout << "After again  " << y << std::endl;

   }

};

TEST_BEGIN(UnitCellTest)
TEST_ADD(UnitCellTest, test1DLamellar)
TEST_ADD(UnitCellTest, test2DSquare)
TEST_ADD(UnitCellTest, test2DHexagonal)
TEST_ADD(UnitCellTest, test3DOrthorhombic)
TEST_ADD(UnitCellTest, test3DCubic)
TEST_END(UnitCellTest)

#endif
