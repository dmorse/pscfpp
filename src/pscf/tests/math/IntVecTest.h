#ifndef PSCF_INT_VEC_TEST_H
#define PSCF_INT_VEC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/IntVec.h>

#include <iostream>
#include <fstream>

using namespace Pscf;

class IntVecTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      int a[3]  = {1, 3, 2};
      IntVec<3> v(a);
      TEST_ASSERT(1 == v[0]);
      TEST_ASSERT(3 == v[1]);
      TEST_ASSERT(2 == v[2]);
   }
 
   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      IntVec<3> v(3);
      TEST_ASSERT(3 == v[0]);
      TEST_ASSERT(3 == v[1]);
      TEST_ASSERT(3 == v[2]);
      int a[3] = {3, 3, 3};
      IntVec<3> u(a);
      TEST_ASSERT(u == v);
   }

   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      int va[3] = {1, -3, 2};
      IntVec<3> v(va);
      IntVec<3> u(v);
      TEST_ASSERT(u == v);
   }
 
   void testEquality()
   {
      printMethod(TEST_FUNC);
      int a[3]  = {1, -4, 2 };
      IntVec<3> v(a);
      IntVec<3> u(v);
      TEST_ASSERT(v == u);
      v[0] = v[0] + 5;
      TEST_ASSERT(v != u);
   }
 
   void testComparison()
   {
      printMethod(TEST_FUNC);
      int a[3]  = {2, -3, 2};
      int b[3]  = {2, -3, 1};
      int c[3]  = {1, -3, 1};
      IntVec<3> va(a);
      IntVec<3> vb(b);
      IntVec<3> vc(c);
      TEST_ASSERT(vc < va);
      TEST_ASSERT(va > vc);
      TEST_ASSERT(!(vc > va));
      TEST_ASSERT(!(va < vc));
      TEST_ASSERT(!(vc < vc));
      TEST_ASSERT(!(vc > vc));
      TEST_ASSERT(vc < vb);
      TEST_ASSERT(vb > vc);
      TEST_ASSERT(!(vc > vb));
      TEST_ASSERT(!(vb < vc));
      TEST_ASSERT(!(vb < vb));
      TEST_ASSERT(!(vb > vb));
      TEST_ASSERT(vb < va);
      TEST_ASSERT(va > vb);
      TEST_ASSERT(!(vb > va));
      TEST_ASSERT(!(va < vb));
   }
 
   void testAssignment()
   {
      printMethod(TEST_FUNC);
      int a[3]  = { 1, -3, 2 };
      IntVec<3> v(a);
      IntVec<3> u = v;
      TEST_ASSERT(v == u);
      v[0] = v[0] + 1;
      TEST_ASSERT(v != u);
   }
 
   void testAdd() 
   {
      printMethod(TEST_FUNC);
      int ua[3] = {1, 4, 2};
      int va[3] = {2, 1, -2};
      int ea[3] = {3, 5, 0};
      IntVec<3> u(ua);
      IntVec<3> v(va);
      IntVec<3> e(ea);
      IntVec<3> r;
      r.add(u,v);
      TEST_ASSERT(r == e);
   }

   void testSubtract() 
   {
      printMethod(TEST_FUNC);
      int ua[3] = {1, -4,  2};
      int va[3] = {3,  1, -6};
      int ea[3] = {-2, -5, 8};
      IntVec<3> u(ua);
      IntVec<3> v(va);
      IntVec<3> e(ea);
      IntVec<3> r;
      r.subtract(u,v);
      TEST_ASSERT(r == e);
   }

   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      int ua[3] = {1,  -4,  2};
      int ea[3] = {3, -12,  6};
      IntVec<3> u(ua);
      IntVec<3> e(ea);
      IntVec<3> r;
      r.multiply(u, 3);
      TEST_ASSERT(r == e);
   }

   void testDot() 
   {
      printMethod(TEST_FUNC);
      int a[3] = { 1, -2, 2 };
      int b[3] = { 2,  3, 4 };
      int d;
      IntVec<3> v(a);
      IntVec<3> u(b);
      d = dot(u, v);
      TEST_ASSERT(d == 4);
   }

   #if 0
   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();
      IntVec<3> v;
      std::ifstream in;
      openInputFile("in/IntVec", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;
   }
   #endif

};

TEST_BEGIN(IntVecTest)
TEST_ADD(IntVecTest, testConstructor1)
TEST_ADD(IntVecTest, testConstructor2)
TEST_ADD(IntVecTest, testCopyConstructor)
TEST_ADD(IntVecTest, testEquality)
TEST_ADD(IntVecTest, testComparison)
TEST_ADD(IntVecTest, testAssignment)
TEST_ADD(IntVecTest, testAdd)
TEST_ADD(IntVecTest, testSubtract)
TEST_ADD(IntVecTest, testMultiply)
TEST_ADD(IntVecTest, testDot)
//TEST_ADD(IntVecTest, testReadWrite)
TEST_END(IntVecTest)

#endif
