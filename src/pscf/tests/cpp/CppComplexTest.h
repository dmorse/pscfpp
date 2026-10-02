#ifndef PSCF_CPP_COMPLEX_TEST_H
#define PSCF_CPP_COMPLEX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/complex.h>

using namespace Util;
using namespace Pscf;

class CppComplexTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testAddCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 0.25;
         b[1] = 1.5;
         add(z, a, b);
         TEST_ASSERT(eq(z[0], 2.25));
         TEST_ASSERT(eq(z[1], 2.00));
         TEST_ASSERT(eq(z[0], a[0] + b[0]));
         TEST_ASSERT(eq(z[1], a[1] + b[1]));
      }
   }

   void testAddCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         add(z, a, b);
         TEST_ASSERT(eq(z[0], 2.25));
         TEST_ASSERT(eq(z[1], 0.50));
         TEST_ASSERT(eq(z[0], a[0] + b));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testAddEqCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 0.25;
         b[1] = 1.5;
         z[0] = a[0];
         z[1] = a[1];
         addEq(z, b);
         TEST_ASSERT(eq(z[0], 2.25));
         TEST_ASSERT(eq(z[1], 2.00));
         TEST_ASSERT(eq(z[0], a[0] + b[0]));
         TEST_ASSERT(eq(z[1], a[1] + b[1]));
      }
   }

   void testAddEqCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         z[0] = a[0];
         z[1] = a[1];
         addEq(z, b);
         TEST_ASSERT(eq(z[0], 2.25));
         TEST_ASSERT(eq(z[1], 0.50));
         TEST_ASSERT(eq(z[0], a[0] + b));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testSubCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 0.25;
         b[1] = 1.5;
         sub(z, a, b);
         TEST_ASSERT(eq(z[0], 1.75));
         TEST_ASSERT(eq(z[1], -1.00));
         TEST_ASSERT(eq(z[0], a[0] - b[0]));
         TEST_ASSERT(eq(z[1], a[1] - b[1]));
      }
   }

   void testSubCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex a;
         double b;
         fftw_complex z;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         sub(z, a, b);
         TEST_ASSERT(eq(z[0], 1.75));
         TEST_ASSERT(eq(z[1], 0.50));
         TEST_ASSERT(eq(z[0], a[0] - b));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testSubEqCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 0.25;
         b[1] = 1.5;
         z[0] = a[0];
         z[1] = a[1];
         subEq(z, b);
         TEST_ASSERT(eq(z[0],  1.75));
         TEST_ASSERT(eq(z[1], -1.00));
         TEST_ASSERT(eq(z[0], a[0] - b[0]));
         TEST_ASSERT(eq(z[1], a[1] - b[1]));
      }
   }

   void testSubEqCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         z[0] = a[0];
         z[1] = a[1];
         subEq(z, b);
         TEST_ASSERT(eq(z[0], 1.75));
         TEST_ASSERT(eq(z[1], 0.50));
         TEST_ASSERT(eq(z[0], a[0] - b));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testMulCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 3.0;
         b[1] = 2.0;
         mul(z, a, b);
         TEST_ASSERT(eq(z[0], 5.0));
         TEST_ASSERT(eq(z[1], 5.5));
         TEST_ASSERT(eq(z[0], a[0]*b[0] - a[1]*b[1]));
         TEST_ASSERT(eq(z[1], a[1]*b[0] + a[0]*b[1]));

         mulEq(a, b);
         TEST_ASSERT(eq(z[0], a[0]));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testMulCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         mul(z, a, b);
         TEST_ASSERT(eq(z[0], 0.50));
         TEST_ASSERT(eq(z[1], 0.125));
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         mulEq(a, b);
         TEST_ASSERT(eq(z[0], a[0]));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testMulEqCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 3.0;
         b[1] = 2.0;
         z[0] = a[0];
         z[1] = a[1];
         mulEq(z, b);
         TEST_ASSERT(eq(z[0], 5.0));
         TEST_ASSERT(eq(z[1], 5.5));
         TEST_ASSERT(eq(z[0], a[0]*b[0] - a[1]*b[1]));
         TEST_ASSERT(eq(z[1], a[1]*b[0] + a[0]*b[1]));
      }
   }

   void testMulEqCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         z[0] = a[0];
         z[1] = a[1];
         mulEq(z, b);
         TEST_ASSERT(eq(z[0], 0.50));
         TEST_ASSERT(eq(z[1], 0.125));
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));
      }
   }

   void testDivCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 3.0;
         b[1] = 2.0;
         mul(z, a, b);
         TEST_ASSERT(eq(z[0], a[0]*b[0] - a[1]*b[1]));
         TEST_ASSERT(eq(z[1], a[1]*b[0] + a[0]*b[1]));
 
         fftw_complex x;
         div(x, z, b);
         TEST_ASSERT(eq(x[0], a[0]));
         TEST_ASSERT(eq(x[1], a[1]));
      }
   }

   void testDivCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         mul(z, a, b);
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         fftw_complex x;
         div(x, z, b);
         TEST_ASSERT(eq(x[0], a[0]));
         TEST_ASSERT(eq(x[1], a[1]));
      }
   }

   void testDivEqCc()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         fftw_complex b;
         a[0] = 2.0;
         a[1] = 0.5;
         b[0] = 3.0;
         b[1] = 2.0;
         z[0] = a[0];
         z[1] = a[1];
         mulEq(z, b);
         TEST_ASSERT(eq(z[0], 5.0));
         TEST_ASSERT(eq(z[1], 5.5));

         divEq(z, b);
         TEST_ASSERT(eq(z[0], a[0]));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

   void testDivEqCr()
   {
      printMethod(TEST_FUNC);
      {
         fftw_complex z;
         fftw_complex a;
         double b;
         a[0] = 2.0;
         a[1] = 0.5;
         b    = 0.25;
         z[0] = a[0];
         z[1] = a[1];
         mulEq(z, b);
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         divEq(z, b);
         TEST_ASSERT(eq(z[0], a[0]));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

};

TEST_BEGIN(CppComplexTest)
TEST_ADD(CppComplexTest, testAddCc)
TEST_ADD(CppComplexTest, testAddCr)
TEST_ADD(CppComplexTest, testAddEqCc)
TEST_ADD(CppComplexTest, testAddEqCr)
TEST_ADD(CppComplexTest, testSubCc)
TEST_ADD(CppComplexTest, testSubCr)
TEST_ADD(CppComplexTest, testSubEqCc)
TEST_ADD(CppComplexTest, testSubEqCr)
TEST_ADD(CppComplexTest, testMulCc)
TEST_ADD(CppComplexTest, testMulCr)
TEST_ADD(CppComplexTest, testMulEqCc)
TEST_ADD(CppComplexTest, testMulEqCr)
TEST_ADD(CppComplexTest, testDivCc)
TEST_ADD(CppComplexTest, testDivCr)
TEST_ADD(CppComplexTest, testDivEqCc)
TEST_ADD(CppComplexTest, testDivEqCr)
TEST_END(CppComplexTest)

#endif
