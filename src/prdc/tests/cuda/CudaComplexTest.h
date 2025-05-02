#ifndef PSCF_CUDA_COMPLEX_TEST_H
#define PSCF_CUDA_COMPLEX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/complex.h>

using namespace Util;
using namespace Pscf;
using namespace Prdc;

class CudaComplexTest : public UnitTest
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
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 0.25;
         b.y = 1.5;
         add(z, a, b);
         TEST_ASSERT(eq(z.x, 2.25));
         TEST_ASSERT(eq(z.y, 2.00));
         TEST_ASSERT(eq(z.x, a.x + b.x));
         TEST_ASSERT(eq(z.y, a.y + b.y));
      }
   }

   void testAddCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         add(z, a, b);
         TEST_ASSERT(eq(z.x, 2.25));
         TEST_ASSERT(eq(z.y, 0.50));
         TEST_ASSERT(eq(z.x, a.x + b));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testAddEqCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 0.25;
         b.y = 1.5;
         z.x = a.x;
         z.y = a.y;
         addEq(z, b);
         TEST_ASSERT(eq(z.x, 2.25));
         TEST_ASSERT(eq(z.y, 2.00));
         TEST_ASSERT(eq(z.x, a.x + b.x));
         TEST_ASSERT(eq(z.y, a.y + b.y));
      }
   }

   void testAddEqCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         z.x = a.x;
         z.y = a.y;
         addEq(z, b);
         TEST_ASSERT(eq(z.x, 2.25));
         TEST_ASSERT(eq(z.y, 0.50));
         TEST_ASSERT(eq(z.x, a.x + b));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testSubCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 0.25;
         b.y = 1.5;
         sub(z, a, b);
         TEST_ASSERT(eq(z.x, 1.75));
         TEST_ASSERT(eq(z.y, -1.00));
         TEST_ASSERT(eq(z.x, a.x - b.x));
         TEST_ASSERT(eq(z.y, a.y - b.y));
      }
   }

   void testSubCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex a;
         cudaReal b;
         cudaComplex z;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         sub(z, a, b);
         TEST_ASSERT(eq(z.x, 1.75));
         TEST_ASSERT(eq(z.y, 0.50));
         TEST_ASSERT(eq(z.x, a.x - b));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testSubEqCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 0.25;
         b.y = 1.5;
         z.x = a.x;
         z.y = a.y;
         subEq(z, b);
         TEST_ASSERT(eq(z.x,  1.75));
         TEST_ASSERT(eq(z.y, -1.00));
         TEST_ASSERT(eq(z.x, a.x - b.x));
         TEST_ASSERT(eq(z.y, a.y - b.y));
      }
   }

   void testSubEqCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         z.x = a.x;
         z.y = a.y;
         subEq(z, b);
         TEST_ASSERT(eq(z.x, 1.75));
         TEST_ASSERT(eq(z.y, 0.50));
         TEST_ASSERT(eq(z.x, a.x - b));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testMulCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 3.0;
         b.y = 2.0;
         mul(z, a, b);
         TEST_ASSERT(eq(z.x, 5.0));
         TEST_ASSERT(eq(z.y, 5.5));
         TEST_ASSERT(eq(z.x, a.x * b.x - a.y * b.y));
         TEST_ASSERT(eq(z.y, a.y * b.x + a.x * b.y));

         mulEq(a, b);
         TEST_ASSERT(eq(z.x, a.x));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testMulCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         mul(z, a, b);
         TEST_ASSERT(eq(z.x, 0.50));
         TEST_ASSERT(eq(z.y, 0.125));
         TEST_ASSERT(eq(z.x, a.x * b));
         TEST_ASSERT(eq(z.y, a.y * b));

         mulEq(a, b);
         TEST_ASSERT(eq(z.x, a.x));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testMulEqCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 3.0;
         b.y = 2.0;
         z.x = a.x;
         z.y = a.y;
         mulEq(z, b);
         TEST_ASSERT(eq(z.x, 5.0));
         TEST_ASSERT(eq(z.y, 5.5));
         TEST_ASSERT(eq(z.x, a.x * b.x - a.y * b.y));
         TEST_ASSERT(eq(z.y, a.y * b.x + a.x * b.y));
      }
   }

   void testMulEqCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         z.x = a.x;
         z.y = a.y;
         mulEq(z, b);
         TEST_ASSERT(eq(z.x, 0.50));
         TEST_ASSERT(eq(z.y, 0.125));
         TEST_ASSERT(eq(z.x, a.x * b));
         TEST_ASSERT(eq(z.y, a.y * b));
      }
   }

   void testDivCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 3.0;
         b.y = 2.0;
         mul(z, a, b);
         TEST_ASSERT(eq(z.x, a.x * b.x - a.y * b.y));
         TEST_ASSERT(eq(z.y, a.y * b.x + a.x * b.y));
 
         cudaComplex x;
         div(x, z, b);
         TEST_ASSERT(eq(x.x, a.x));
         TEST_ASSERT(eq(x.y, a.y));
      }
   }

   void testDivCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         mul(z, a, b);
         TEST_ASSERT(eq(z.x, a.x * b));
         TEST_ASSERT(eq(z.y, a.y * b));

         cudaComplex x;
         div(x, z, b);
         TEST_ASSERT(eq(x.x, a.x));
         TEST_ASSERT(eq(x.y, a.y));
      }
   }

   void testDivEqCc()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaComplex b;
         a.x = 2.0;
         a.y = 0.5;
         b.x = 3.0;
         b.y = 2.0;
         z.x = a.x;
         z.y = a.y;
         mulEq(z, b);
         TEST_ASSERT(eq(z.x, 5.0));
         TEST_ASSERT(eq(z.y, 5.5));

         divEq(z, b);
         TEST_ASSERT(eq(z.x, a.x));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

   void testDivEqCr()
   {
      printMethod(TEST_FUNC);
      {
         cudaComplex z;
         cudaComplex a;
         cudaReal b;
         a.x = 2.0;
         a.y = 0.5;
         b = 0.25;
         z.x = a.x;
         z.y = a.y;
         mulEq(z, b);
         TEST_ASSERT(eq(z.x, a.x * b));
         TEST_ASSERT(eq(z.y, a.y * b));

         divEq(z, b);
         TEST_ASSERT(eq(z.x, a.x));
         TEST_ASSERT(eq(z.y, a.y));
      }
   }

};

TEST_BEGIN(CudaComplexTest)
TEST_ADD(CudaComplexTest, testAddCc)
TEST_ADD(CudaComplexTest, testAddCr)
TEST_ADD(CudaComplexTest, testAddEqCc)
TEST_ADD(CudaComplexTest, testAddEqCr)
TEST_ADD(CudaComplexTest, testSubCc)
TEST_ADD(CudaComplexTest, testSubCr)
TEST_ADD(CudaComplexTest, testSubEqCc)
TEST_ADD(CudaComplexTest, testSubEqCr)
TEST_ADD(CudaComplexTest, testMulCc)
TEST_ADD(CudaComplexTest, testMulCr)
TEST_ADD(CudaComplexTest, testMulEqCc)
TEST_ADD(CudaComplexTest, testMulEqCr)
TEST_ADD(CudaComplexTest, testDivCc)
TEST_ADD(CudaComplexTest, testDivCr)
TEST_ADD(CudaComplexTest, testDivEqCc)
TEST_ADD(CudaComplexTest, testDivEqCr)
TEST_END(CudaComplexTest)

#endif
