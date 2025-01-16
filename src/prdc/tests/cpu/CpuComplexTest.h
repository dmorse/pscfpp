#ifndef PSCF_CPU_COMPLEX_TEST_H
#define PSCF_CPU_COMPLEX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/complex.h>

using namespace Util;
using namespace Pscf;

class CpuComplexTest : public UnitTest
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
         Cpu::add(z, a, b);
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
         Cpu::add(z, a, b);
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
         Cpu::addEq(z, b);
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
         Cpu::addEq(z, b);
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
         Cpu::sub(z, a, b);
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
         Cpu::sub(z, a, b);
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
         Cpu::subEq(z, b);
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
         Cpu::subEq(z, b);
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
         Cpu::mul(z, a, b);
         TEST_ASSERT(eq(z[0], 5.0));
         TEST_ASSERT(eq(z[1], 5.5));
         TEST_ASSERT(eq(z[0], a[0]*b[0] - a[1]*b[1]));
         TEST_ASSERT(eq(z[1], a[1]*b[0] + a[0]*b[1]));

         Cpu::mulEq(a, b);
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
         Cpu::mul(z, a, b);
         TEST_ASSERT(eq(z[0], 0.50));
         TEST_ASSERT(eq(z[1], 0.125));
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         Cpu::mulEq(a, b);
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
         Cpu::mulEq(z, b);
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
         Cpu::mulEq(z, b);
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
         Cpu::mul(z, a, b);
         TEST_ASSERT(eq(z[0], a[0]*b[0] - a[1]*b[1]));
         TEST_ASSERT(eq(z[1], a[1]*b[0] + a[0]*b[1]));
 
         fftw_complex x;
         Cpu::div(x, z, b);
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
         Cpu::mul(z, a, b);
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         fftw_complex x;
         Cpu::div(x, z, b);
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
         Cpu::mulEq(z, b);
         TEST_ASSERT(eq(z[0], 5.0));
         TEST_ASSERT(eq(z[1], 5.5));

         Cpu::divEq(z, b);
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
         Cpu::mulEq(z, b);
         TEST_ASSERT(eq(z[0], a[0]*b));
         TEST_ASSERT(eq(z[1], a[1]*b));

         Cpu::divEq(z, b);
         TEST_ASSERT(eq(z[0], a[0]));
         TEST_ASSERT(eq(z[1], a[1]));
      }
   }

};

TEST_BEGIN(CpuComplexTest)
TEST_ADD(CpuComplexTest, testAddCc)
TEST_ADD(CpuComplexTest, testAddCr)
TEST_ADD(CpuComplexTest, testAddEqCc)
TEST_ADD(CpuComplexTest, testAddEqCr)
TEST_ADD(CpuComplexTest, testSubCc)
TEST_ADD(CpuComplexTest, testSubCr)
TEST_ADD(CpuComplexTest, testSubEqCc)
TEST_ADD(CpuComplexTest, testSubEqCr)
TEST_ADD(CpuComplexTest, testMulCc)
TEST_ADD(CpuComplexTest, testMulCr)
TEST_ADD(CpuComplexTest, testMulEqCc)
TEST_ADD(CpuComplexTest, testMulEqCr)
TEST_ADD(CpuComplexTest, testDivCc)
TEST_ADD(CpuComplexTest, testDivCr)
TEST_ADD(CpuComplexTest, testDivEqCc)
TEST_ADD(CpuComplexTest, testDivEqCr)
TEST_END(CpuComplexTest)

#endif
