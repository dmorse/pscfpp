#ifndef PSCF_REAL_VEC_TEST_H
#define PSCF_REAL_VEC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/math/RealVec.h>

#include <iostream>
#include <fstream>

using namespace Pscf;

class RealVecTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      double a[3]  = {1.1, 3.0, -2.2};
      RealVec<3> v(a);
      TEST_ASSERT(1.1 == v[0]);
      TEST_ASSERT(3.0 == v[1]);
      TEST_ASSERT(-2.2 == v[2]);
   }
 
   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      RealVec<3> v(3.0);
      TEST_ASSERT(3 == v[0]);
      TEST_ASSERT(3 == v[1]);
      TEST_ASSERT(3 == v[2]);
      double a[3] = {3.0, 3.0, 3.0};
      RealVec<3> u(a);
      TEST_ASSERT(u == v);
   }

   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      double va[3] = {1.0, -3.0, 2.2};
      RealVec<3> v(va);
      RealVec<3> u(v);
      TEST_ASSERT(u == v);
   }
 
   void testEquality()
   {
      printMethod(TEST_FUNC);
      double a[3]  = {1.3, -4.8, 2.0};
      RealVec<3> v(a);
      RealVec<3> u(v);
      TEST_ASSERT(v == u);
      v[0] = v[0] + 5.0;
      TEST_ASSERT(v != u);
   }
 
   void testAssignment()
   {
      printMethod(TEST_FUNC);
      double a[3]  = {1.0, -3.2, 2.0 };
      RealVec<3> v(a);
      RealVec<3> u = v;
      TEST_ASSERT(v == u);
      v[0] = v[0] + 1.9;
      TEST_ASSERT(v != u);
   }
 
   void testAdd() 
   {
      printMethod(TEST_FUNC);
      double ua[3] = {1.4, 4.1, 2.0};
      double va[3] = {2.3, 1.0, -2.0};
      double ea[3] = {3.7, 5.1, 0.0};
      RealVec<3> u(ua);
      RealVec<3> v(va);
      RealVec<3> e(ea);
      RealVec<3> r;
      r.add(u,v);
      TEST_ASSERT(r == e);
   }

   void testSubtract() 
   {
      printMethod(TEST_FUNC);
      double ua[3] = {1.0, -4.0,  2.0};
      double va[3] = {3.2,  1.0, -6.2};
      double ea[3] = {-2.2, -5.0, 8.2};
      RealVec<3> u(ua);
      RealVec<3> v(va);
      RealVec<3> e(ea);
      RealVec<3> r;
      r.subtract(u,v);
      TEST_ASSERT(r == e);
   }

   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      double ua[3] = {1.0,  -4.0, 2.2};
      double ea[3] = {3.0, -12.0, 6.6};
      RealVec<3> u(ua);
      RealVec<3> e(ea);
      RealVec<3> r;
      r.multiply(u, 3);
      TEST_ASSERT(r == e);
   }

   void testDot() 
   {
      printMethod(TEST_FUNC);
      double a[3] = {1.0, -2.0, 2.0};
      double b[3] = {2.0, 3.0, 4.0 };
      double d;
      RealVec<3> v(a);
      RealVec<3> u(b);
      d = dot(u, v);
      TEST_ASSERT(eq(d,4.0));
   }

   #if 0
   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();
      RealVec<3> v;
      std::ifstream in;
      openInputFile("in/RealVec", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;
   }
   #endif

};

TEST_BEGIN(RealVecTest)
TEST_ADD(RealVecTest, testConstructor1)
TEST_ADD(RealVecTest, testConstructor2)
TEST_ADD(RealVecTest, testCopyConstructor)
TEST_ADD(RealVecTest, testEquality)
TEST_ADD(RealVecTest, testAssignment)
TEST_ADD(RealVecTest, testAdd)
TEST_ADD(RealVecTest, testSubtract)
TEST_ADD(RealVecTest, testMultiply)
TEST_ADD(RealVecTest, testDot)
//TEST_ADD(RealVecTest, testReadWrite)
TEST_END(RealVecTest)

#endif
