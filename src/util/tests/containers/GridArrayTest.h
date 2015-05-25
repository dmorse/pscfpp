#ifndef GRID_ARRAY_TEST_H
#define GRID_ARRAY_TEST_H

#include <util/containers/GridArray.h>
#include <util/space/IntVector.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class GridArrayTest : public UnitTest 
{

   int memory_;

public:

   void setUp() { memory_ = Memory::total(); }

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testRankPosition();

   void testIsInGrid();

   void testShift();

   void testSubscript();

   void testCopyConstructor();

   void testAssignment();

   #if 0
   void testSerializeFile1();

   void testSerializeFile2();
   #endif

};

void GridArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   TEST_ASSERT(!v.isAllocated() );
}

void GridArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   TEST_ASSERT(!v.isAllocated());

   IntVector dimensions;
   dimensions[0] = 4;
   dimensions[1] = 3;
   dimensions[2] = 2;
   v.allocate(dimensions);

   TEST_ASSERT(v.size() == 24);
   TEST_ASSERT(v.dimension(0) == 4);
   TEST_ASSERT(v.dimension(1) == 3);
   TEST_ASSERT(v.dimension(2) == 2);
   TEST_ASSERT(v.dimensions() == dimensions);
   TEST_ASSERT(v.isAllocated());

   IntVector position;
   position[0] = 2;
   position[1] = 1;
   position[2] = 1;
   TEST_ASSERT(v.rank(position) == 15);
   TEST_ASSERT(v.position(15) == position);

   position[0] = 3;
   position[1] = 2;
   position[2] = 1;
   TEST_ASSERT(v.rank(position) == 23);
   TEST_ASSERT(v.position(23) == position);

} 

void GridArrayTest::testRankPosition()
{
   printMethod(TEST_FUNC);

   GridArray<int> g; 
   IntVector v;
   v[0] = 2;
   v[1] = 3;
   v[2] = 2;
   g.allocate(v);
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

void GridArrayTest::testIsInGrid()
{
   printMethod(TEST_FUNC);

   GridArray<int> g; 
   IntVector v;
   v[0] = 2;
   v[1] = 3;
   v[2] = 2;
   g.allocate(v);
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

void GridArrayTest::testShift()
{
   printMethod(TEST_FUNC);

   GridArray<int> g; 
   IntVector v;
   v[0] = 2;
   v[1] = 3;
   v[2] = 2;
   g.allocate(v);

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

void GridArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;

   IntVector dimensions;
   dimensions[0] = 4;
   dimensions[1] = 3;
   dimensions[2] = 2;
   v.allocate(dimensions);

   TEST_ASSERT(v.size() == 24);
   TEST_ASSERT(v.dimension(0) == 4);
   TEST_ASSERT(v.dimension(1) == 3);
   TEST_ASSERT(v.dimension(2) == 2);
   TEST_ASSERT(v.dimensions() == dimensions);
   TEST_ASSERT(v.isAllocated());

   IntVector position;
   position[0] = 2;
   position[1] = 1;
   position[2] = 1;
   v(position) = 38;
   TEST_ASSERT(v.rank(position) == 15);
   TEST_ASSERT(v.position(15) == position);
   TEST_ASSERT(v[15] == 38);
   TEST_ASSERT(v(position) == 38);
} 

void GridArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;

   IntVector dimensions;
   dimensions[0] = 3;
   dimensions[1] = 1;
   dimensions[2] = 2;
   v.allocate(dimensions);
   TEST_ASSERT(v.isAllocated() );
   TEST_ASSERT(v.size() == 6);

   IntVector position1;
   position1[0] = 0;
   position1[1] = 0;
   position1[2] = 0;
   v(position1) = 38;

   IntVector position2;
   position2[0] = 0;
   position2[1] = 0;
   position2[2] = 1;
   v(position2) = 15;

   IntVector position3;
   position3[0] = 1;
   position3[1] = 0;
   position3[2] = 0;
   v(position3) = 13;

   IntVector position4;
   position4[0] = 1;
   position4[1] = 0;
   position4[2] = 1;
   v(position4) = 9;

   IntVector position5;
   position5[0] = 2;
   position5[1] = 0;
   position5[2] = 0;
   v(position5) = 7;

   IntVector position6;
   position6[0] = 2;
   position6[1] = 0;
   position6[2] = 1;
   v(position6) = 5;

   GridArray<int> u(v);

   TEST_ASSERT(u.isAllocated());
   TEST_ASSERT(u.dimensions() == v.dimensions());
   TEST_ASSERT(u.size() == v.size() );
   TEST_ASSERT(u(position1) == 38 );
   TEST_ASSERT(u(position2) == 15 );
   TEST_ASSERT(u(position3) == 13 );
   TEST_ASSERT(u(position4) == 9);
   TEST_ASSERT(u(position5) == 7);
   TEST_ASSERT(u(position6) == 5);

} 

void GridArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;

   IntVector dimensions;
   dimensions[0] = 3;
   dimensions[1] = 1;
   dimensions[2] = 2;
   v.allocate(dimensions);
   TEST_ASSERT(v.isAllocated() );
   TEST_ASSERT(v.size() == 6);

   IntVector position1;
   position1[0] = 0;
   position1[1] = 0;
   position1[2] = 0;
   v(position1) = 38;

   IntVector position2;
   position2[0] = 0;
   position2[1] = 0;
   position2[2] = 1;
   v(position2) = 15;

   IntVector position3;
   position3[0] = 1;
   position3[1] = 0;
   position3[2] = 0;
   v(position3) = 13;

   IntVector position4;
   position4[0] = 1;
   position4[1] = 0;
   position4[2] = 1;
   v(position4) = 9;

   IntVector position5;
   position5[0] = 2;
   position5[1] = 0;
   position5[2] = 0;
   v(position5) = 7;

   IntVector position6;
   position6[0] = 2;
   position6[1] = 0;
   position6[2] = 1;
   v(position6) = 5;

   GridArray<int> u;
   u = v;

   TEST_ASSERT(u.isAllocated());
   TEST_ASSERT(u.dimensions() == v.dimensions());
   TEST_ASSERT(u.size() == v.size() );
   TEST_ASSERT(u(position1) == 38 );
   TEST_ASSERT(u(position2) == 15 );
   TEST_ASSERT(u(position3) == 13 );
   TEST_ASSERT(u(position4) == 9);
   TEST_ASSERT(u(position5) == 7);
   TEST_ASSERT(u(position6) == 5);

} 

#if 0
void GridArrayTest::testSerializeFile1()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   int i1 = 13;
   int i2;

   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   GridArray<int> u;
   u.allocate(2, 2);

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(u.capacity1() == 2);
   TEST_ASSERT(u.capacity2() == 2);
   TEST_ASSERT(u(0,0) == 3);
   TEST_ASSERT(u(1,0) == 4);
   TEST_ASSERT(u(0,1) == 5);
   TEST_ASSERT(u(1,1) == 6);
   TEST_ASSERT(i2 == i1);
   TEST_ASSERT(i2 == 13);

   #if 0
   // Clear values of u and i2
   for (int i=0; i < capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);
   #endif

} 

void GridArrayTest::testSerializeFile2()
{
   printMethod(TEST_FUNC);
   GridArray<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   int i1 = 13;
   int i2;

   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   GridArray<int> u;

   //u.allocate(2, 2);  -> Note allocation, different from previous

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(u.capacity1() == 2);
   TEST_ASSERT(u.capacity2() == 2);
   TEST_ASSERT(u(0,0) == 3);
   TEST_ASSERT(u(1,0) == 4);
   TEST_ASSERT(u(0,1) == 5);
   TEST_ASSERT(u(1,1) == 6);
   TEST_ASSERT(i2 == i1);
   TEST_ASSERT(i2 == 13);

   #if 0
   // Clear values of u and i2
   for (int i=0; i < capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);
   #endif

}
#endif // if 0
 
TEST_BEGIN(GridArrayTest)
TEST_ADD(GridArrayTest, testConstructor)
TEST_ADD(GridArrayTest, testAllocate)
TEST_ADD(GridArrayTest, testRankPosition)
TEST_ADD(GridArrayTest, testIsInGrid)
TEST_ADD(GridArrayTest, testShift)
TEST_ADD(GridArrayTest, testSubscript)
TEST_ADD(GridArrayTest, testCopyConstructor)
TEST_ADD(GridArrayTest, testAssignment)
//TEST_ADD(GridArrayTest, testSerializeFile1)
//TEST_ADD(GridArrayTest, testSerializeFile2)
TEST_END(GridArrayTest)

#endif
