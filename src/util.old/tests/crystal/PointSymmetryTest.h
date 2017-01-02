#ifndef POINT_SYMMETRY_TEST_H
#define POINT_SYMMETRY_TEST_H

#include <util/crystal/PointSymmetry.h>
#include <util/space/IntVector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fstream>

using namespace Util;

class PointSymmetryTest : public UnitTest 
{

public:

   void testConstructor()
   {
      printMethod(TEST_FUNC);
      PointSymmetry a();
   } 

   void testCopyConstructor()
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      a.R(0,1) = 1;
      a.R(1,0) = 1;
      a.R(2,2) = 1;

      PointSymmetry b(a);
      TEST_ASSERT(b == a);

      //std::cout << std::endl;
      //std::cout << b;

   }

   void testAssignment()
   {
      printMethod(TEST_FUNC);

      PointSymmetry a;
      PointSymmetry b;

      a.R(0,1) = 1;
      a.R(1,0) = 1;
      a.R(2,2) = 1;
      TEST_ASSERT(b != a);

      b = a;
      TEST_ASSERT(b == a);

      //std::cout << std::endl;
      //std::cout << b;

   }

   void testEquality()
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      a.R(2,1) = 1;
      a.R(2,1) = 1;
      a.R(0,0) = 1;

      PointSymmetry b;
      b.R(2,1) = 1;
      b.R(2,1) = 1;
      b.R(0,0) = 1;

      TEST_ASSERT(b == a);
      TEST_ASSERT(a == b);
   }
 
   void testInEquality()
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      a.R(2,1) = 1;
      a.R(2,1) = 1;
      a.R(0,0) = 1;

      PointSymmetry b;
      b.R(0,1) = 1;
      b.R(1,0) = 1;
      b.R(2,2) = 1;

      TEST_ASSERT(!(b == a));
      TEST_ASSERT(b != a);
      TEST_ASSERT(a != b);
   }
 
   void testMultiply() 
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      a.R(2,1) =  1;
      a.R(1,2) =  1;
      a.R(0,0) = -1;

      PointSymmetry b;
      b.R(0,1) = 1;
      b.R(1,0) = 1;
      b.R(2,2) = 1;

      PointSymmetry c;
      c = a*b;

      TEST_ASSERT(c.R(0,0) == 0);
      TEST_ASSERT(c.R(0,1) == -1);
      TEST_ASSERT(c.R(0,2) ==  0);
      TEST_ASSERT(c.R(1,0) ==  0);
      TEST_ASSERT(c.R(1,1) ==  0);
      TEST_ASSERT(c.R(1,2) ==  1);
      TEST_ASSERT(c.R(2,0) ==  1);
      TEST_ASSERT(c.R(2,1) ==  0);
      TEST_ASSERT(c.R(2,2) ==  0);

      //std::cout << std::endl;
      //std::cout << c;

   }
 
   void testIdentity() 
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      PointSymmetry b;
   
      a = PointSymmetry::identity();
      TEST_ASSERT(a != b);

      //std::cout << std::endl;
      //std::cout << a;

      b = PointSymmetry::identity();
      TEST_ASSERT(a == b);

      //std::cout << std::endl;
      //std::cout << b;
 
      
   }

   void testInverse() 
   {
      printMethod(TEST_FUNC);
      PointSymmetry a;
      a.R(2,1) =  1;
      a.R(1,2) =  1;
      a.R(0,0) = -1;

      PointSymmetry b;
      b = a.inverse();

      PointSymmetry c;
      c = a*b;

      TEST_ASSERT(c == PointSymmetry::identity());

      //std::cout << std::endl;
      //std::cout << c;

   }
 
   void testMultiplyVector() 
   {
      printMethod(TEST_FUNC);
      PointSymmetry s;
      s.R(2,1) =  1;
      s.R(1,2) =  1;
      s.R(0,0) = -1;

      IntVector v;
      v[0] =  1;
      v[1] = -2;
      v[2] =  3;

      IntVector u;
      u = s*v;

      TEST_ASSERT(u[0] == -1);
      TEST_ASSERT(u[1] ==  3);
      TEST_ASSERT(u[2] == -2);

   }
 
   void testVectorMultiply() 
   {
      printMethod(TEST_FUNC);
      PointSymmetry s;
      s.R(2,1) =  1;
      s.R(1,2) = -1;
      s.R(0,0) = -1;

      IntVector v;
      v[0] =  1;
      v[1] = -2;
      v[2] =  3;

      IntVector u;
      u = v*s;

      TEST_ASSERT(u[0] == -1);
      TEST_ASSERT(u[1] ==  3);
      TEST_ASSERT(u[2] ==  2);

   }
 
};

TEST_BEGIN(PointSymmetryTest)
TEST_ADD(PointSymmetryTest, testConstructor)
TEST_ADD(PointSymmetryTest, testCopyConstructor)
TEST_ADD(PointSymmetryTest, testAssignment)
TEST_ADD(PointSymmetryTest, testEquality)
TEST_ADD(PointSymmetryTest, testInEquality)
TEST_ADD(PointSymmetryTest, testMultiply)
TEST_ADD(PointSymmetryTest, testIdentity)
TEST_ADD(PointSymmetryTest, testInverse)
TEST_ADD(PointSymmetryTest, testMultiplyVector)
TEST_ADD(PointSymmetryTest, testVectorMultiply)
TEST_END(PointSymmetryTest)

#endif
