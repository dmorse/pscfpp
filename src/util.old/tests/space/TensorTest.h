#ifndef TENSOR_TEST_H
#define TENSOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/space/Tensor.h>

#include <fstream>

using namespace Util;

class TensorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      Tensor v(0.5);
      Tensor u(0.5);
      TEST_ASSERT(u == v);
      for (int i=0; i < Dimension; ++i) {
         for (int j=0; j < Dimension; ++j) {
            TEST_ASSERT(eq(v(i,j), 0.5));
            TEST_ASSERT(eq(u(i,j), 0.5));
         }
      }
   } 

   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      Tensor v(0.0);
      v(0,1) = 1.0;
      v(1,2) = 0.5;
      v(2,0) = 2.0;
      double a[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}, {2.0, 0.0, 0.0}};
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
   }

   void testConstructor3()
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}, {2.0, 0.0, 0.0}};
      Tensor v(a);
      v(0, 0) = 0.0; 
      v(0, 1) = 1.0; 
      v(0, 2) = 0.0; 
      v(0, 0) = 0.0; 
      v(1, 1) = 0.0; 
      v(1, 2) = 0.5; 
      v(2, 0) = 2.0; 
      v(2, 1) = 0.0; 
      v(2, 2) = 0.0; 
      TEST_ASSERT(v == a);
   }
 
   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      Tensor v(0.0);
      v(0,1) = 1.0;
      v(1,2) = 0.5;
      v(2,0) = 2.0;
      Tensor u(v);
      u(0, 0) = 0.0; 
      u(0, 1) = 1.0; 
      u(0, 2) = 0.0; 
      u(0, 0) = 0.0; 
      u(1, 1) = 0.0; 
      u(1, 2) = 0.5; 
      u(2, 0) = 2.0; 
      u(2, 1) = 0.0; 
      u(2, 2) = 0.0; 
      TEST_ASSERT(u == v);
   }
 
   void testEquality()
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}, {2.0, 0.0, 0.0}};
      Tensor v(a);
      const Tensor u(v);
      TEST_ASSERT(eq(u(0, 0), 0.0)); 
      TEST_ASSERT(eq(u(0, 1), 1.0)); 
      TEST_ASSERT(eq(u(0, 2), 0.0)); 
      TEST_ASSERT(eq(u(0, 0), 0.0)); 
      TEST_ASSERT(eq(u(1, 1), 0.0)); 
      TEST_ASSERT(eq(u(1, 2), 0.5)); 
      TEST_ASSERT(eq(u(2, 0), 2.0)); 
      TEST_ASSERT(eq(u(2, 1), 0.0)); 
      TEST_ASSERT(eq(u(2, 2), 0.0)); 
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
      TEST_ASSERT(v == u);
      v(0,1) = v(0, 1) + 0.1;
      TEST_ASSERT(eq(v(0, 1), 1.1));
      TEST_ASSERT(v != u);
      TEST_ASSERT(a != v);
      TEST_ASSERT(v != a);
   }
 
   void testAssignment()
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.5}, {2.0, 0.0, 0.0}};
      Tensor v(a);
      Tensor u = v;
      TEST_ASSERT(v == u);
      v(1,2) = v(1,2) + 0.1;
      TEST_ASSERT(v != u);
   }
 
   void testSetRow()
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{-0.2, 1.0, 0.3}, {-0.7, 0.8, 0.0}, {2.0, 0.0, 0.0}};
      Tensor v(a);
      Tensor u = v;
      double b[3] = {3.0, 4.0, 5.0};
      Vector t(b);
      u.setRow(1, t);
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(1,i), t[i]));
      }
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(0,i), v(0,i)));
      }
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(2,i), v(2,i)));
      }
   }
 
   void testSetColumn()
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{-0.2, 1.0, 0.3}, {-0.7, 0.8, 0.0}, {2.0, 0.0, 0.0}};
      Tensor v(a);
      Tensor u = v;
      double b[3] = {3.0, 4.0, 5.0};
      Vector t(b);
      u.setColumn(1, t);
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(i, 1), t[i]));
      }
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(i, 0), v(i, 0)));
      }
      for (int i = 0; i < 3; ++i) {
         TEST_ASSERT(eq(u(i, 2), v(i, 2)));
      }
   }

   void testAdd() 
   {
      printMethod(TEST_FUNC);
      double a1[3][3] = {{0.0, 1.0, 0.0}, {0.3, 0.0, 0.5}, { 2.0, 0.5, 0.0}};
      double a2[3][3] = {{0.2, 0.0, 2.0}, {1.0, 3.2, 0.5}, {-2.0, 0.0, 1.0}};
      double a3[3][3] = {{0.2, 1.0, 2.0}, {1.3, 3.2, 1.0}, { 0.0, 0.5, 1.0}};
      Tensor u(a1);
      Tensor v(a2);
      Tensor e(a3);
      Tensor r;
      r.add(u,v);
      TEST_ASSERT(r == e);
   }

   void testSubtract() 
   {
      printMethod(TEST_FUNC);
      double a1[3][3] = {{0.0, 1.0, 0.0}, {0.3, 1.0, 0.5}, { 2.0, 0.5, 4.0}};
      double a2[3][3] = {{0.2, 0.0, 2.0}, {1.0, 3.2, 0.5}, {-2.0, 0.0, 1.0}};
      double a3[3][3] = {{-0.2, 1.0, -2.0}, {-0.7, -2.2, 0.0}, {4.0, 0.5, 3.0}};
      Tensor u(a1);
      Tensor v(a2);
      Tensor e(a3);
      Tensor r;
      r.subtract(u, v);
      TEST_ASSERT(r == e);
   }

   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      double a1[3][3] = {{0.1, 1.0,-2.0}, {0.3, 1.0, 0.5}, {2.0, -0.5, 0.4}};
      double a3[3][3] = {{0.3, 3.0, -6.0}, {0.9, 3.0, 1.5}, {6.0, -1.5, 1.2}};
      Tensor u(a1);
      Tensor e(a3);
      Tensor r;
      r.multiply(u, 3.0);
      TEST_ASSERT(r == e);
   }

   void testDivide() 
   {
      printMethod(TEST_FUNC);
      double a1[3][3] = {{0.1, 1.0,-2.0}, {0.3, 1.0, 0.5}, {2.0, -0.5, 0.4}};
      double a3[3][3] = {{0.3, 3.0, -6.0}, {0.9, 3.0, 1.5}, {6.0, -1.5, 1.2}};
      Tensor u(a1);
      Tensor e(a3);
      Tensor r;
      r.divide(e, 3.0);
      TEST_ASSERT(r == u);
   }

   void testTrace() 
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.1, 1.0,-2.0}, {0.3, 1.0, 0.5}, {2.0, -0.5, 0.4}};
      Tensor u(a);
      double trace = u.trace();
      TEST_ASSERT( eq(trace, 1.5) );
   }

   void testDyad() 
   {
      printMethod(TEST_FUNC);
      double a[3] = {1.0, -2.0, 3.0};
      double b[3] = {2.0, 3.0, 4.0};
      double c[3][3] = {{2.0, 3.0, 4.0}, {-4.0, -6.0, -8.0}, {6.0, 9.0, 12.0}};
      Vector x(a);
      Vector y(b);
      Tensor u;
      Tensor v(c);

      u.dyad(x, y);
      TEST_ASSERT(u == v);
   }



   void testSymmetrizeSelf() 
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.1, 1.0,-2.0}, {0.4, 1.0, 0.3}, {2.0, -0.5, 0.4}};
      double b[3][3] = {{0.1, 0.7, 0.0}, {0.7, 1.0,-0.1}, {0.0, -0.1, 0.4}};
      Tensor u(a);
      Tensor v(b);
      u.symmetrize();
      TEST_ASSERT( u == v );
   }

   void testSymmetrizeOther() 
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.1, 1.0,-2.0}, {0.4, 1.0, 0.3}, {2.0, -0.5, 0.4}};
      double b[3][3] = {{0.1, 0.7, 0.0}, {0.7, 1.0,-0.1}, {0.0, -0.1, 0.4}};
      Tensor u(a);
      Tensor v(b);
      Tensor w;
      w.symmetrize(u);
      TEST_ASSERT( w == v );
   }

   void testTransposeSelf() 
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.1, 1.0,-2.0}, {0.8, 3.0,  0.3}, { 2.0, -0.5, 0.4}};
      double b[3][3] = {{0.1, 0.8, 2.0}, {1.0, 3.0, -0.5}, {-2.0,  0.3, 0.4}};
      Tensor u(a);
      Tensor v(b);
      u.transpose();
      TEST_ASSERT( u == v );
   }

   void testTransposeOther() 
   {
      printMethod(TEST_FUNC);
      double a[3][3] = {{0.1, 1.0,-2.0}, {0.8, 3.0,  0.3}, { 2.0, -0.5, 0.4}};
      double b[3][3] = {{0.1, 0.8, 2.0}, {1.0, 3.0, -0.5}, {-2.0,  0.3, 0.4}};
      Tensor u(a);
      Tensor v(b);
      Tensor w;
      w.transpose(u);
      TEST_ASSERT( w == v );
   }

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Tensor v;
      std::ifstream in;
      openInputFile("in/Tensor", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(TensorTest)
TEST_ADD(TensorTest, testConstructor1)
TEST_ADD(TensorTest, testConstructor2)
TEST_ADD(TensorTest, testConstructor3)
TEST_ADD(TensorTest, testCopyConstructor)
TEST_ADD(TensorTest, testEquality)
TEST_ADD(TensorTest, testAssignment)
TEST_ADD(TensorTest, testSetRow)
TEST_ADD(TensorTest, testSetColumn)
TEST_ADD(TensorTest, testAdd)
TEST_ADD(TensorTest, testSubtract)
TEST_ADD(TensorTest, testMultiply)
TEST_ADD(TensorTest, testDivide)
TEST_ADD(TensorTest, testTrace)
TEST_ADD(TensorTest, testDyad)
TEST_ADD(TensorTest, testSymmetrizeSelf)
TEST_ADD(TensorTest, testSymmetrizeOther)
TEST_ADD(TensorTest, testTransposeSelf)
TEST_ADD(TensorTest, testTransposeOther)
TEST_ADD(TensorTest, testReadWrite)
TEST_END(TensorTest)

#endif
