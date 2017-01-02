#ifndef DMATRIX_TEST_H
#define DMATRIX_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DMatrix.h>
#include <util/containers/Matrix.h>

using namespace Util;

class DMatrixTest : public UnitTest 
{

   int memory_;

public:

   void setUp() 
   {  memory_ = Memory::total(); }

   void tearDown() {}
  
   void testConstructor();

   void testAllocate();

   void testSubscript();

   void testCopyConstructor();

   void testAssignment();

   void testBaseClassReference();

   void testSerializeFile1();

   void testSerializeFile2();

};


void DMatrixTest::testConstructor()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   TEST_ASSERT(v.capacity1() == 0 );
   TEST_ASSERT(v.capacity2() == 0 );
   TEST_ASSERT(!v.isAllocated() );
} 

void DMatrixTest::testAllocate()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   v.allocate(2,3);
   TEST_ASSERT(v.capacity1() == 2 );
   TEST_ASSERT(v.capacity2() == 3 );
   TEST_ASSERT(v.isAllocated() );
} 

void DMatrixTest::testSubscript()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
   v.allocate(2,2);
   v(0,0) = 3;
   v(1,0) = 4;
   v(0,1) = 5;
   v(1,1) = 6;
   TEST_ASSERT(v(0,0) == 3 );
   TEST_ASSERT(v(1,0) == 4 );
   TEST_ASSERT(v(0,1) == 5 );
   TEST_ASSERT(v(1,1) == 6 );
} 

void DMatrixTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
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

   DMatrix<int> u(v);
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

void DMatrixTest::testAssignment()
{
   printMethod(TEST_FUNC);

   DMatrix<int> v;
   v.allocate(2,2);
   TEST_ASSERT(v.capacity1() == 2);
   TEST_ASSERT(v.capacity2() == 2);
   TEST_ASSERT(v.isAllocated() );

   DMatrix<int> u;
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

void DMatrixTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
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

void DMatrixTest::testSerializeFile1()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
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

   DMatrix<int> u;
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

void DMatrixTest::testSerializeFile2()
{
   printMethod(TEST_FUNC);
   DMatrix<int> v;
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

   DMatrix<int> u;

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
TEST_BEGIN(DMatrixTest)
TEST_ADD(DMatrixTest, testConstructor)
TEST_ADD(DMatrixTest, testAllocate)
TEST_ADD(DMatrixTest, testSubscript)
TEST_ADD(DMatrixTest, testCopyConstructor)
TEST_ADD(DMatrixTest, testAssignment)
TEST_ADD(DMatrixTest, testBaseClassReference)
TEST_ADD(DMatrixTest, testSerializeFile1)
TEST_ADD(DMatrixTest, testSerializeFile2)
TEST_END(DMatrixTest)

#endif
