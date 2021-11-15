#ifndef PSCF_SPACE_SYMMETRY_TEST_H
#define PSCF_SPACE_SYMMETRY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/SpaceSymmetry.h>
#include <util/math/Constants.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class SpaceSymmetryTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   template <int D>
   bool isValid(SpaceSymmetry<D> cell)
   {
      return true;
   }

   void test2DIdentity() 
   {
      printMethod(TEST_FUNC);
      printEndl();
      SpaceSymmetry<2> E = SpaceSymmetry<2>::identity();
      TEST_ASSERT(E.R(0,0) == 1);
      TEST_ASSERT(E.R(1,0) == 0);
      TEST_ASSERT(E.R(0,1) == 0);
      TEST_ASSERT(E.R(1,1) == 1);
      TEST_ASSERT(E.t(0) == 0);
      TEST_ASSERT(E.t(1) == 0);
      //std::cout << SpaceSymmetry<2>::identity() << std::endl;
   }

   void test2DConstruct() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      SpaceSymmetry<2> A;
      A.R(0,0) = 0;
      A.R(1,0) = 1;
      A.R(0,1) = -1;
      A.R(1,1) = 0;
      A.t(0) = 0;
      A.t(1) = Rational(1, 2);

      std::cout << A << std::endl;
   }

   void test2DRead() 
   {
      printMethod(TEST_FUNC);
      printEndl();
     
      std::ifstream in;
      openInputFile("in/Symmetry2D", in);

      SpaceSymmetry<2> A;
      in >> A;
      std::cout << A;
   }

   void test2DEquality() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      SpaceSymmetry<2> A;
      A.R(0,0) = 0;
      A.R(1,0) = 1;
      A.R(0,1) = -1;
      A.R(1,1) = 0;
      A.t(0) = 0;
      A.t(1) = Rational(1, 2);
      //std::cout << A << std::endl;

      SpaceSymmetry<2> B = A;
      //std::cout << B << std::endl;

      TEST_ASSERT(A == B);

      B.t(0) = Rational(1, 2);
      TEST_ASSERT(A != B);

      A.t(0) = Rational(1, 2);
      TEST_ASSERT(A == B);

      A.R(0,0) = 2;
      TEST_ASSERT(A != B);

   }

   void test2DInvertMultiply() 
   {
      printMethod(TEST_FUNC);
      // printEndl();

      SpaceSymmetry<2> A;
      A.R(0,0) = 0;
      A.R(1,0) = 1;
      A.R(0,1) = -1;
      A.R(1,1) = 0;
      A.t(0) = 0;
      A.t(1) = Rational(1, 2);
      //std::cout << A << std::endl;

      SpaceSymmetry<2> B;
      B = A.inverse();
      //std::cout << B << std::endl;
      
      SpaceSymmetry<2> C;
      C = A*B;
      //std::cout << C << std::endl;
      TEST_ASSERT(C == A*B);
      TEST_ASSERT(C == SpaceSymmetry<2>::identity());

      SpaceSymmetry<2> D;
      D = A*B;
      TEST_ASSERT(D == C);

      SpaceSymmetry<2> E;
      E = A*A;
      //std::cout << E << std::endl;
      TEST_ASSERT(A == B*E);
      TEST_ASSERT(A == E*B);
   
   }
  
   void test3DInvertMultiply() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceSymmetry<3> A;
      A.R(0,0) =  0;
      A.R(0,1) = -1;
      A.R(0,2) =  0;
      A.R(1,0) =  1;
      A.R(1,1) =  0;
      A.R(1,2) =  0;
      A.R(2,0) =  0;
      A.R(2,1) =  0;
      A.R(2,2) =  1;
      A.t(0) = 0;
      A.t(1) = Rational(1,2);
      A.t(2) = Rational(-1,4);
      //std::cout << A << std::endl;
      SpaceSymmetry<3> B = A.inverse();
      //std::cout << B << std::endl;
      SpaceSymmetry<3> C = A*B;
      // std::cout << C << std::endl;
      TEST_ASSERT(C == SpaceSymmetry<3>::identity());

      A.R(0,0) =  1;
      A.R(0,1) = -1;
      A.R(0,2) =  0;
      A.R(1,0) =  0;
      A.R(1,1) = -1;
      A.R(1,2) =  0;
      A.R(2,0) =  0;
      A.R(2,1) =  0;
      A.R(2,2) =  1;
      A.t(0) = 0;
      A.t(1) = Rational(1,2);
      A.t(2) = Rational(-1,4);
      //std::cout << A << std::endl;
      B = A.inverse();
      //std::cout << B << std::endl;
      C = A*B;
      //std::cout << C << std::endl;
      TEST_ASSERT(C == SpaceSymmetry<3>::identity());

      A.R(0,0) =  0;
      A.R(0,1) =  1;
      A.R(0,2) =  0;
      A.R(1,0) =  0;
      A.R(1,1) =  0;
      A.R(1,2) =  1;
      A.R(2,0) =  1;
      A.R(2,1) =  0;
      A.R(2,2) =  0;
      A.t(0) = 0;
      A.t(1) = Rational(-3,2);
      A.t(2) = Rational(-1,4);
      B = A.inverse();
      //std::cout << B << std::endl;
      C = A*B;
      //std::cout << C << std::endl;
      TEST_ASSERT(C == SpaceSymmetry<3>::identity());
   }

   void testShiftOrigin() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceSymmetry<3> A;
      A.R(0,0) =  0;
      A.R(0,1) = -1;
      A.R(0,2) =  0;
      A.R(1,0) =  1;
      A.R(1,1) =  0;
      A.R(1,2) =  0;
      A.R(2,0) =  0;
      A.R(2,1) =  0;
      A.R(2,2) =  1;
      A.t(0) = 0;
      A.t(1) = Rational(1,2);
      A.t(2) = Rational(-1,4);
      A.normalize();

      SpaceSymmetry<3> B = A;
      TEST_ASSERT(B == A);
      //TEST_ASSERT(A == B);
 
      SpaceSymmetry<3>::Translation origin;
      origin[0] = Rational(1, 8);
      origin[1] = Rational(-3, 4);
      origin[2] = Rational(0, 1);

      B.shiftOrigin(origin);
      TEST_ASSERT(A != B);

      origin[0] = Rational(-1, 8);
      origin[1] = Rational(3, 4);
      origin[2] = Rational(0, 1);
      B.shiftOrigin(origin);

      TEST_ASSERT(A == B);
   }
};

TEST_BEGIN(SpaceSymmetryTest)
TEST_ADD(SpaceSymmetryTest, test2DIdentity)
TEST_ADD(SpaceSymmetryTest, test2DConstruct)
TEST_ADD(SpaceSymmetryTest, test2DRead)
TEST_ADD(SpaceSymmetryTest, test2DEquality)
TEST_ADD(SpaceSymmetryTest, test2DInvertMultiply)
TEST_ADD(SpaceSymmetryTest, test3DInvertMultiply)
TEST_ADD(SpaceSymmetryTest, testShiftOrigin)
TEST_END(SpaceSymmetryTest)

#endif
