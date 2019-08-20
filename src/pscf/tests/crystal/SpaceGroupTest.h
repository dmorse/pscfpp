#ifndef PSCF_SPACE_GROUP_TEST_H
#define PSCF_SPACE_GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/SpaceGroup.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class SpaceGroupTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   void test2DmakeIdentity() 
   {
      printMethod(TEST_FUNC);
      //printEndl();
      SpaceGroup<2> G;
      G.makeCompleteGroup();
      TEST_ASSERT(1 == G.size());
      //std::cout << G << std::endl;
   }

   void test2Dmake1() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceGroup<2> G;
      SpaceSymmetry<2> A;
      A.R(0,0) = 0;
      A.R(1,0) = 1;
      A.R(0,1) = -1;
      A.R(1,1) = 0;
      A.t(0) = 0;
      A.t(1) = Rational(1, 2);

      G.add(A);
      G.makeCompleteGroup();
      TEST_ASSERT(4 == G.size());
      //std::cout << G << std::endl;
   }

   void test2Dmake2() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceGroup<2> G;

      // Below: Add generating elements of simple Hexagonal group

      SpaceSymmetry<2> A;
      A.R(0,0) = -1;
      A.R(1,0) =  0;
      A.R(0,1) =  0;
      A.R(1,1) = -1;
      A.t(0) = 0;
      A.t(1) = 0;
      G.add(A);

      A.R(0,0) = -1;
      A.R(1,0) =  1;
      A.R(0,1) =  0;
      A.R(1,1) =  1;
      A.t(0) = 0;
      A.t(1) = 0;
      G.add(A);

      A.R(0,0) =  0;
      A.R(1,0) =  1;
      A.R(0,1) = -1;
      A.R(1,1) =  1;
      A.t(0) = 0;
      A.t(1) = 0;
      G.add(A);

      G.makeCompleteGroup();
      //std::cout << G << std::endl;
      TEST_ASSERT(12 == G.size());
   }

   void test2Dmake3() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceGroup<2> G;
      SpaceSymmetry<2> A;

      A.R(0,0) = 0;
      A.R(1,0) = 1;
      A.R(0,1) = 1;
      A.R(1,1) = 0;
      A.t(0) = 0;
      A.t(1) = 0;
      G.add(A);

      A.R(0,0) = -1;
      A.R(1,0) = 0;
      A.R(0,1) = 0;
      A.R(1,1) = 1;
      A.t(0) = 0;
      A.t(1) = 0;
      G.add(A);

      A.R(0,0) = 1;
      A.R(1,0) = 0;
      A.R(0,1) = 0;
      A.R(1,1) = 1;
      A.t(0) = Rational(1, 2);
      A.t(1) = A.t(0);
      G.add(A);

      G.makeCompleteGroup();
      //std::cout << G << std::endl;
      TEST_ASSERT(16 == G.size());
   }

   void test3Dmake() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      SpaceGroup<3> G;

      // Below: Add generating elements of BCC group

      SpaceSymmetry<3> A;
      A.R(0,0) = -1;
      A.R(0,1) = 0;
      A.R(0,2) = 0;
      A.R(1,0) = 0;
      A.R(1,1) = 1;
      A.R(1,2) = 0;
      A.R(2,0) = 0;
      A.R(2,1) = 0;
      A.R(2,2) = 1;
      A.t(0) = 0;
      A.t(1) = 0;
      A.t(2) = 0;
      //std::cout << A << std::endl;

      G.add(A);
      //G.makeCompleteGroup();
      //std::cout << G << std::endl;

      A.R(0,0) = 0;
      A.R(0,1) = 1;
      A.R(0,2) = 0;
      A.R(1,0) = 1;
      A.R(1,1) = 0;
      A.R(1,2) = 0;
      A.R(2,0) = 0;
      A.R(2,1) = 0;
      A.R(2,2) = 1;
      A.t(0) = 0;
      A.t(1) = 0;
      A.t(2) = 0;

      G.add(A);
      //G.makeCompleteGroup();
      //std::cout << G << std::endl;

      A.R(0,0) = 0;
      A.R(0,1) = 0;
      A.R(0,2) = 1;
      A.R(1,0) = 0;
      A.R(1,1) = 1;
      A.R(1,2) = 0;
      A.R(2,0) = 1;
      A.R(2,1) = 0;
      A.R(2,2) = 0;
      A.t(0) = 0;
      A.t(1) = 0;
      A.t(2) = 0;
      G.add(A);
      //G.makeCompleteGroup();
      //std::cout << G << std::endl;
      //std::cout << "size =" << G.size() << std::endl;

      A.R(0,0) = 1;
      A.R(0,1) = 0;
      A.R(0,2) = 0;
      A.R(1,0) = 0;
      A.R(1,1) = 1;
      A.R(1,2) = 0;
      A.R(2,0) = 0;
      A.R(2,1) = 0;
      A.R(2,2) = 1;
      A.t(0) = Rational(1,2);
      A.t(1) = Rational(1,2);
      A.t(2) = Rational(1,2);
      G.add(A);

      G.makeCompleteGroup();
      TEST_ASSERT(96 == G.size());

      //std::cout << G << std::endl;
   }

   void test2Dread() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      std::ifstream in;
      openInputFile("in/p_6_m_m", in);

      SpaceGroup<2> g;
      in >> g;
      TEST_ASSERT(12 == g.size());
      TEST_ASSERT(g.isValid());

      // std::cout << std::endl;
      //std::cout << g << std::endl;
   }

};

TEST_BEGIN(SpaceGroupTest)
TEST_ADD(SpaceGroupTest, test2DmakeIdentity)
TEST_ADD(SpaceGroupTest, test2Dmake1)
TEST_ADD(SpaceGroupTest, test2Dmake2)
TEST_ADD(SpaceGroupTest, test2Dmake3)
TEST_ADD(SpaceGroupTest, test3Dmake)
TEST_ADD(SpaceGroupTest, test2Dread)
TEST_END(SpaceGroupTest)

#endif
