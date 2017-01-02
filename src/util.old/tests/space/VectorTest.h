#ifndef VECTOR_TEST_H
#define VECTOR_TEST

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/space/Vector.h>

#include <fstream>

using namespace Util;

class VectorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      Vector v(0.0);
      Vector u(0.0, 0.0, 0.0);
      TEST_ASSERT(u == v);
   } 

   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      double a[3]  = { 1.0, 0.5, 2.0 };
      Vector v(a);
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
   }
 
   void testConstructor3()
   {
      printMethod(TEST_FUNC);
      double a[3]  = { 1.0, 0.5, 2.0 };
      Vector v(a);
      Vector u(1.0, 0.5, 2.0);
      TEST_ASSERT(u == v);
   }
 
   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      Vector v(1.0, 0.5, 2.0);
      Vector u(v);
      TEST_ASSERT(u == v);
   }
 
   void testEquality()
   {
      printMethod(TEST_FUNC);
      double a[3]  = { 1.0, 0.5, 2.0 };
      Vector v(a);
      Vector u(v);
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
      TEST_ASSERT(v == u);
      v[0] = v[0] + 0.1;
      TEST_ASSERT(v != u);
      TEST_ASSERT(a != v);
      TEST_ASSERT(v != a);
   }
 
   void testAssignment()
   {
      printMethod(TEST_FUNC);
      double a[3]  = { 1.0, 0.5, 2.0 };
      Vector v(a);
      Vector u = v;
      TEST_ASSERT(v == u);
      v[0] = v[0] + 0.1;
      TEST_ASSERT(v != u);
   }
 
   void testAdd() 
   {
      printMethod(TEST_FUNC);
      Vector u(1.0, 0.4,  2.0);
      Vector v(0.2, 1.0, -2.1);
      Vector e(1.2, 1.4, -0.1);
      Vector r;
      r.add(u,v);
      TEST_ASSERT(r == e);
   }

   void testSubtract() 
   {
      printMethod(TEST_FUNC);
      Vector u(1.0, 0.4,  2.0);
      Vector v(0.2, 1.0, -2.1);
      Vector e(0.8,-0.6,  4.1);
      Vector r;
      r.subtract(u,v);
      TEST_ASSERT(r == e);
   }

   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      Vector u(1.0, 0.4,  2.0);
      Vector e(3.0, 1.2,  6.0);
      Vector r;
      r.multiply(u, 3.0);
      TEST_ASSERT(r == e);
   }

   void testDivide() 
   {
      printMethod(TEST_FUNC);
      Vector u(1.0, 0.4,  2.0);
      Vector e(0.5, 0.2,  1.0);
      Vector r;
      r.divide(u, 2.0);
      TEST_ASSERT(r == e);
   }

   void testDot() 
   {
      printMethod(TEST_FUNC);
      double a[3] = { 1.0, 0.5, 2.0 };
      double b[3] = { 0.5, 1.0, 2.0 };
      double d;
      Vector v(a);
      Vector u(b);
      d = u.dot(v);
      TEST_ASSERT(eq(d, 5.0));
   }

   void testCross() 
   {
      printMethod(TEST_FUNC);
      Vector u(0.5, 1.0, 2.0);
      Vector v(1.0, 0.0, 0.0);
      Vector w;
      w.cross(u,v);
      TEST_ASSERT(w == Vector(0.0, 2.0, -1.0));
      TEST_ASSERT(Vector(0.0, 2.0, -1.0) == w);
   }

   void testSquare() 
   {
      printMethod(TEST_FUNC);
      Vector v(1.0, 0.5, 2.0);
      double c;
      c = v.square();
      TEST_ASSERT(eq(c , 5.25));
   }

   void testAbs() 
   {
      Vector v(1.0, 0.5, 2.0);
      double c;
      printMethod(TEST_FUNC);
      c = v.abs();
      TEST_ASSERT(eq(c*c , 5.25));
   }

   void testVersor() 
   {
      Vector v(1.0, 0.5, 2.0);
      Vector u;
      double d;
      printMethod(TEST_FUNC);
      u.versor(v);
      d = v.abs();
      u *= d;
      TEST_ASSERT(u == v);
   }

   void testTransverse() {
      printMethod(TEST_FUNC);
      Vector v(1.0, 0.5, 2.0);
      Vector u(0.5, 0.5, 0.0);
      Vector t;
      Vector p;
      Vector w;
      p.parallel(v,u);
      t.transverse(v,u);
      w.add(p,t);
      TEST_ASSERT(w == v);
      TEST_ASSERT(eq(u.dot(t), 0.0));
   }

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Vector v;
      std::ifstream in;
      openInputFile("in/Vector", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(VectorTest)
TEST_ADD(VectorTest, testConstructor1)
TEST_ADD(VectorTest, testConstructor2)
TEST_ADD(VectorTest, testConstructor3)
TEST_ADD(VectorTest, testCopyConstructor)
TEST_ADD(VectorTest, testEquality)
TEST_ADD(VectorTest, testAssignment)
TEST_ADD(VectorTest, testAdd)
TEST_ADD(VectorTest, testSubtract)
TEST_ADD(VectorTest, testMultiply)
TEST_ADD(VectorTest, testDivide)
TEST_ADD(VectorTest, testDot)
TEST_ADD(VectorTest, testCross)
TEST_ADD(VectorTest, testSquare)
TEST_ADD(VectorTest, testAbs)
TEST_ADD(VectorTest, testVersor)
TEST_ADD(VectorTest, testTransverse)
TEST_ADD(VectorTest, testReadWrite)
TEST_END(VectorTest)

#endif
