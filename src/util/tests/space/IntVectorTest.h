#ifndef INT_VECTOR_TEST_H
#define INT_VECTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/space/IntVector.h>

#include <iostream>
#include <fstream>

using namespace Util;

class IntVectorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      IntVector v(0);
      IntVector u(0, 0, 0);
      TEST_ASSERT(u == v);
   } 

   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      int a[3]  = { 1, 3, 2 };
      IntVector v(a);
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
   }
 
   void testConstructor3()
   {
      printMethod(TEST_FUNC);
      int a[3]  = { 1, -1, 2 };
      IntVector v(a);
      IntVector u(1, -1, 2);
      TEST_ASSERT(u == v);
   }
 
   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      IntVector v(1, -3, 2);
      IntVector u(v);
      TEST_ASSERT(u == v);
   }
 
   void testEquality()
   {
      printMethod(TEST_FUNC);
      int a[3]  = { 1, -4, 2 };
      IntVector v(a);
      IntVector u(v);
      TEST_ASSERT(v == a);
      TEST_ASSERT(a == v);
      TEST_ASSERT(v == u);
      v[0] = v[0] + 5;
      TEST_ASSERT(v != u);
      TEST_ASSERT(a != v);
      TEST_ASSERT(v != a);
   }
 
   void testAssignment()
   {
      printMethod(TEST_FUNC);
      int a[3]  = { 1, -3, 2 };
      IntVector v(a);
      IntVector u = v;
      TEST_ASSERT(v == u);
      v[0] = v[0] + 1;
      TEST_ASSERT(v != u);
   }
 
   void testAdd() 
   {
      printMethod(TEST_FUNC);
      IntVector u(1, 4,  2);
      IntVector v(2, 1, -2);
      IntVector e(3, 5,  0);
      IntVector r;
      r.add(u,v);
      TEST_ASSERT(r == e);
   }

   void testSubtract() 
   {
      printMethod(TEST_FUNC);
      IntVector u(1, -4,  2);
      IntVector v(3,  1, -6);
      IntVector e(-2, -5, 8);
      IntVector r;
      r.subtract(u,v);
      TEST_ASSERT(r == e);
   }

   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      IntVector u(1,  -4,  2);
      IntVector e(3, -12,  6);
      IntVector r;
      r.multiply(u, 3);
      TEST_ASSERT(r == e);
   }

   void testDot() 
   {
      printMethod(TEST_FUNC);
      int a[3] = { 1, -2, 2 };
      int b[3] = { 2,  3, 4 };
      int d;
      IntVector v(a);
      IntVector u(b);
      d = u.dot(v);
      TEST_ASSERT(d == 4);
   }

   void testCross() 
   {
      printMethod(TEST_FUNC);
      IntVector u(2, 1, 3);
      IntVector v(1, 0, -1);
      IntVector w;
      w.cross(u,v);
      TEST_ASSERT(w == IntVector(-1, 5, -1) );
      TEST_ASSERT(IntVector(-1, 5, -1) == w);
   }

   void testSquare() 
   {
      printMethod(TEST_FUNC);
      IntVector v(1, 3, 2);
      int c;
      c = v.square();
      TEST_ASSERT(c == 14);
   }

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();
      IntVector v;
      std::ifstream in;
      openInputFile("in/IntVector", in);
      //std::ifstream in("in/IntVector");

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

   }

};

TEST_BEGIN(IntVectorTest)
TEST_ADD(IntVectorTest, testConstructor1)
TEST_ADD(IntVectorTest, testConstructor2)
TEST_ADD(IntVectorTest, testConstructor3)
TEST_ADD(IntVectorTest, testCopyConstructor)
TEST_ADD(IntVectorTest, testEquality)
TEST_ADD(IntVectorTest, testAssignment)
TEST_ADD(IntVectorTest, testAdd)
TEST_ADD(IntVectorTest, testSubtract)
TEST_ADD(IntVectorTest, testMultiply)
TEST_ADD(IntVectorTest, testDot)
TEST_ADD(IntVectorTest, testCross)
TEST_ADD(IntVectorTest, testSquare)
TEST_ADD(IntVectorTest, testReadWrite)
TEST_END(IntVectorTest)

#endif
