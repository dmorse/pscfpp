#ifndef GRID_TEST_H
#define GRID_TEST

#include <util/space/Grid.h>
#include <util/space/IntVector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <fstream>

using namespace Util;

class GridTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor1()
   {
      printMethod(TEST_FUNC);
      Grid g;
      TEST_ASSERT(g.dimension(2) == 1);
      TEST_ASSERT(g.dimension(1) == 1);
      TEST_ASSERT(g.dimension(0) == 1);
      TEST_ASSERT(g.size() == 1);
      IntVector u = g.dimensions();
      TEST_ASSERT(u[2] == 1);
      TEST_ASSERT(u[1] == 1);
      TEST_ASSERT(u[0] == 1);
   } 

   void testConstructor2()
   {
      printMethod(TEST_FUNC);
      IntVector v;
      v[0] = 2;
      v[1] = 3;
      v[2] = 2;
      Grid g(v);
      TEST_ASSERT(g.dimension(2) == 2);
      TEST_ASSERT(g.dimension(1) == 3);
      TEST_ASSERT(g.dimension(0) == 2);
      TEST_ASSERT(g.size() == 12);
      IntVector u = g.dimensions();
      TEST_ASSERT(u[2] == 2);
      TEST_ASSERT(u[1] == 3);
      TEST_ASSERT(u[0] == 2);
   }

   void testRankPosition()
   {
      printMethod(TEST_FUNC);
      IntVector v;
      v[0] = 2;
      v[1] = 3;
      v[2] = 2;
      Grid g(v);
      TEST_ASSERT(g.dimension(2) == 2);
      TEST_ASSERT(g.dimension(1) == 3);
      TEST_ASSERT(g.dimension(0) == 2);
      TEST_ASSERT(g.size() == 12);

      IntVector u;
      u[0] = 1;
      u[1] = 1;
      u[2] = 1;
      TEST_ASSERT(g.rank(u) == 9);
      TEST_ASSERT(u == g.position(9));

      u[0] = 1;
      u[1] = 2;
      u[2] = 1;
      TEST_ASSERT(g.rank(u) == 11);
      TEST_ASSERT(u == g.position(g.rank(u)));

      u[0] = 0;
      u[1] = 2;
      u[2] = 1;
      TEST_ASSERT(u == g.position(g.rank(u)));

      TEST_ASSERT(7 == g.rank(g.position(7)));
   }

   void testIsInGrid()
   {
      printMethod(TEST_FUNC);

      IntVector v;
      v[0] = 2;
      v[1] = 3;
      v[2] = 2;
      Grid g(v);
      TEST_ASSERT(g.size() == 12);

      TEST_ASSERT(g.isInGrid(0, 1));
      TEST_ASSERT(g.isInGrid(2, 1));
      TEST_ASSERT(!g.isInGrid(-1, 1));
      TEST_ASSERT(!g.isInGrid(-2, 1));
      TEST_ASSERT(!g.isInGrid(3, 1));
      TEST_ASSERT(!g.isInGrid(4, 1));

      IntVector u;
      u[0] =  1;
      u[1] =  2;
      u[2] =  0;
      TEST_ASSERT(g.isInGrid(u));

      u[0] =  1;
      u[1] =  3;
      u[2] =  0;
      TEST_ASSERT(!g.isInGrid(u));

      u[0] =  1;
      u[1] =  2;
      u[2] =  -1;
      TEST_ASSERT(!g.isInGrid(u));
   }

   void testShift()
   {
      printMethod(TEST_FUNC);

      IntVector v;
      v[0] = 2;
      v[1] = 3;
      v[2] = 2;
      Grid g(v);

      int coordinate = 2;

      coordinate = 0;
      TEST_ASSERT(g.shift(coordinate, 1) == 0);
      TEST_ASSERT(coordinate == 0);

      coordinate = 2;
      TEST_ASSERT(g.shift(coordinate, 1) == 0);
      TEST_ASSERT(coordinate == 2);

      coordinate = 3;
      TEST_ASSERT(g.shift(coordinate, 1) == 1);
      TEST_ASSERT(coordinate == 0);

      coordinate = 4;
      TEST_ASSERT(g.shift(coordinate, 1) == 1);
      TEST_ASSERT(coordinate == 1);

      coordinate = 5;
      TEST_ASSERT(g.shift(coordinate, 1) == 1);
      TEST_ASSERT(coordinate == 2);

      coordinate = -2;
      TEST_ASSERT(g.shift(coordinate, 1) == -1);
      TEST_ASSERT(coordinate == 1);

      coordinate = -3;
      TEST_ASSERT(g.shift(coordinate, 1) == -1);
      TEST_ASSERT(coordinate == 0);

      coordinate = -4;
      TEST_ASSERT(g.shift(coordinate, 1) == -2);
      TEST_ASSERT(coordinate == 2);

      IntVector u, s;
      u[0] =  1;
      u[1] =  4;
      u[2] = -1;
      s = g.shift(u);
      TEST_ASSERT(u[0] == 1);
      TEST_ASSERT(u[1] == 1);
      TEST_ASSERT(u[2] == 1);
      TEST_ASSERT(s[0] ==  0);
      TEST_ASSERT(s[1] ==  1);
      TEST_ASSERT(s[2] == -1);

      u[0] =   2;
      u[1] =  -3;
      u[2] =   5;
      s = g.shift(u);
      TEST_ASSERT(u[0] == 0);
      TEST_ASSERT(u[1] == 0);
      TEST_ASSERT(u[2] == 1);
      TEST_ASSERT(s[0] ==  1);
      TEST_ASSERT(s[1] == -1);
      TEST_ASSERT(s[2] ==  2);
   }

};

TEST_BEGIN(GridTest)
TEST_ADD(GridTest, testConstructor1)
TEST_ADD(GridTest, testConstructor2)
TEST_ADD(GridTest, testRankPosition)
TEST_ADD(GridTest, testIsInGrid)
TEST_ADD(GridTest, testShift)
TEST_END(GridTest)

#endif
