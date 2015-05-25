#ifndef POINT_GROUP_TEST_H
#define POINT_GROUP_TEST_H

#include <util/crystal/PointGroup.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fstream>

using namespace Util;

class PointGroupTest : public UnitTest 
{

public:

   void testConstructor()
   {
      printMethod(TEST_FUNC);

      PointGroup g;
      TEST_ASSERT(g.size() == 1);
   } 

   void testAdd()
   {
      printMethod(TEST_FUNC);

      PointGroup g;
      TEST_ASSERT(g.size() == 1);

      PointSymmetry        a, b, c;
      const PointSymmetry* p;
      bool result;

      a.R(0,1) = 1;
      a.R(1,0) = 1;
      a.R(2,2) = 1;
      b.R(0,1) = 1;
      b.R(1,0) = 1;
      b.R(2,2) = 1;
      c.R(0,1) =  1;
      c.R(1,0) = -1;
      c.R(2,2) =  1;
      TEST_ASSERT(g.find(a) == 0);
      TEST_ASSERT(g.find(c) == 0);

      result = g.add(a);
      TEST_ASSERT(result);
      TEST_ASSERT(g.size() == 2);

      p = g.find(a);
      TEST_ASSERT(p != 0);
      TEST_ASSERT(*p == a);

      result = g.add(b);
      TEST_ASSERT(!result);
      TEST_ASSERT(g.size() == 2);
      TEST_ASSERT(g.find(b) != 0);

      p = g.find(b);
      TEST_ASSERT(p != 0);
      TEST_ASSERT(*p == b);

      result = g.add(c);
      TEST_ASSERT(result);
      TEST_ASSERT(g.size() == 3);
      TEST_ASSERT(g.find(c) != 0);

      p = g.find(c);
      TEST_ASSERT(p != 0);
      TEST_ASSERT(*p == c);
      
      //std::cout << std::endl;
      //std::cout << g << std::endl;
   }

   void testMakeComplete()
   {
      printMethod(TEST_FUNC);

      PointGroup g;
      TEST_ASSERT(g.size() == 1);

      PointSymmetry a, b, c;

      a.R(0,1) =  1;
      a.R(1,0) =  1;
      a.R(2,2) =  1;

      b.R(0,0) = -1;
      b.R(1,1) =  1;
      b.R(2,2) =  1;

      c.R(0,1) =  1;
      c.R(1,2) =  1;
      c.R(2,0) =  1;

      g.add(c);
      g.makeCompleteGroup();
      //std::cout << std::endl;
      //std::cout << g;
      TEST_ASSERT(g.size() == 3);
      TEST_ASSERT(g.isValid());

      g.add(b);
      g.makeCompleteGroup();
      //std::cout << std::endl;
      //std::cout << g;
      TEST_ASSERT(g.size() == 24);
      TEST_ASSERT(g.isValid());

      g.add(a);
      g.makeCompleteGroup();
      //std::cout << std::endl;
      //std::cout << g;
      TEST_ASSERT(g.size() == 48);
      TEST_ASSERT(g.isValid());

   }

   void testMakeStar()
   {
      printMethod(TEST_FUNC);

      PointGroup g;
      TEST_ASSERT(g.size() == 1);

      PointSymmetry a, b, c;

      a.R(0,1) =  1;
      a.R(1,0) =  1;
      a.R(2,2) =  1;

      b.R(0,0) = -1;
      b.R(1,1) =  1;
      b.R(2,2) =  1;

      c.R(0,1) =  1;
      c.R(1,2) =  1;
      c.R(2,0) =  1;

      g.add(c);
      g.add(b);
      g.add(a);
      g.makeCompleteGroup();
      TEST_ASSERT(g.size() == 48);
      TEST_ASSERT(g.isValid());

      IntVector root;
      FSArray<IntVector, 48> star;

      root[0] = -1;
      root[1] =  3;
      root[2] =  2;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 48);

      root[0] = -1;
      root[1] = -1;
      root[2] =  2;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 24);

      root[0] = -1;
      root[1] = -1;
      root[2] = -1;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 8);

      root[0] = -1;
      root[1] =  0;
      root[2] =  2;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 24);

      root[0] =  3;
      root[1] =  0;
      root[2] =  3;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 12);

      root[0] =  0;
      root[1] =  0;
      root[2] = -1;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 6);

      root[0] =  0;
      root[1] =  0;
      root[2] =  0;
      g.makeStar(root, star);
      TEST_ASSERT(star.size() == 1);

      //std::cout << "size = " << star.size() << std::endl;
      //for (int i = 0; i < star.size(); ++i) {
      //   std::cout << i << "   " << star[i] << std::endl;
      //}
   }

};

TEST_BEGIN(PointGroupTest)
TEST_ADD(PointGroupTest, testConstructor)
TEST_ADD(PointGroupTest, testAdd)
TEST_ADD(PointGroupTest, testMakeComplete)
TEST_ADD(PointGroupTest, testMakeStar)
TEST_END(PointGroupTest)

#endif
