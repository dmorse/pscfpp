#ifndef D_RAGGED_MATRIX_TEST_H
#define D_RAGGED_MATRIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DRaggedMatrix.h>
#include <util/containers/DArray.h>
#include <util/containers/RaggedMatrix.h>

using namespace Util;

class DRaggedMatrixTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testSubscript();

   void testCopyConstructor();

   void testAssignment();

   void testBaseClassReference();

};


void DRaggedMatrixTest::testConstructor()
{
   printMethod(TEST_FUNC);
   DRaggedMatrix<int> v;
   TEST_ASSERT(v.capacity1() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void DRaggedMatrixTest::testAllocate()
{
   printMethod(TEST_FUNC);
   DRaggedMatrix<int> v;
   DArray<int> rowSizes;
   rowSizes.allocate(2);
   rowSizes[0] = 1;
   rowSizes[1] = 3;
   v.allocate(rowSizes);
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2(0) == 1 );
   TEST_ASSERT(v.capacity2(1) == 3 );
   TEST_ASSERT(v.isAllocated() );
} 

void DRaggedMatrixTest::testSubscript()
{
   printMethod(TEST_FUNC);
   DRaggedMatrix<int> v;
   DArray<int> rowSizes;
   rowSizes.allocate(2);
   rowSizes[0] = 1;
   rowSizes[1] = 3;
   v.allocate(rowSizes);
   v(0,0) = 3;
   v(1,0) = 4;
   v(1,1) = 5;
   v(1,2) = 6;
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(1,1) == 5 );
   TEST_ASSERT(v(1,2) == 6 );
} 

#if 0
void DRaggedMatrixTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   DRaggedMatrix<int> v;
   TEST_ASSERT(v.capacity1() == 0 );
   TEST_ASSERT(v.capacity2() == 0 );
   TEST_ASSERT(!v.isAllocated() );

   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   TEST_ASSERT(v.isAllocated() );
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   DRaggedMatrix<int> u(v);
   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u.isAllocated() );
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );
   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );
} 

void DRaggedMatrixTest::testAssignment()
{
   printMethod(TEST_FUNC);

   DRaggedMatrix<int> v;
   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   DRaggedMatrix<int> u;
   u.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   u  = v;

   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u.isAllocated() );

   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );
} 

void DRaggedMatrixTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   DRaggedMatrix<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   Matrix<int>& u = v;

   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );
   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );

}
#endif 

TEST_BEGIN(DRaggedMatrixTest)
TEST_ADD(DRaggedMatrixTest, testConstructor)
TEST_ADD(DRaggedMatrixTest, testAllocate)
TEST_ADD(DRaggedMatrixTest, testSubscript)
   //TEST_ADD(DRaggedMatrixTest, testCopyConstructor);
   //TEST_ADD(DRaggedMatrixTest, testAssignment);
   //TEST_ADD(DRaggedMatrixTest, testBaseClassReference);
TEST_END(DRaggedMatrixTest)

#endif
