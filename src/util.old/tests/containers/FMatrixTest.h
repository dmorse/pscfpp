#ifndef F_MATRIX_TEST_H
#define F_MATRIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/FMatrix.h>
#include <util/containers/Matrix.h>

using namespace Util;

class FMatrixTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testConstructor();

   void testSubscript();

   void testCopyConstructor();

   void testAssignment();

   void testBaseClassReference();

};


void FMatrixTest::testConstructor()
{
   printMethod(TEST_FUNC);
   FMatrix<int,2,2> v;
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
} 

void FMatrixTest::testSubscript()
{
   printMethod(TEST_FUNC);
   FMatrix<int,2,2> v;
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );
} 

void FMatrixTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   FMatrix<int,2,2> v;
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 2 );
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   FMatrix<int,2,2> u(v);
   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );

   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );
} 

void FMatrixTest::testAssignment()
{
   printMethod(TEST_FUNC);

   FMatrix<int,2,2> v;
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;

   FMatrix<int,2,2> u;
   u  = v;

   TEST_ASSERT(u.capacity1() == 2 );
   TEST_ASSERT(u.capacity2() == 2 );

   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );

   TEST_ASSERT(u(0,0) == 3 );
   TEST_ASSERT(u(1,0) == 4 );
   TEST_ASSERT(u(0,1) == 5 );
   TEST_ASSERT(u(1,1) == 6 );

} 

void FMatrixTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   FMatrix<int,2,2> v;
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

TEST_BEGIN(FMatrixTest)
TEST_ADD(FMatrixTest, testConstructor)
TEST_ADD(FMatrixTest, testSubscript)
TEST_ADD(FMatrixTest, testCopyConstructor)
TEST_ADD(FMatrixTest, testAssignment)
TEST_ADD(FMatrixTest, testBaseClassReference)
TEST_END(FMatrixTest)

#endif
