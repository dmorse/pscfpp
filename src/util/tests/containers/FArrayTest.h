#ifndef F_ARRAY_TEST_H
#define F_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/FArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>

#include <complex>

using namespace Util;

class FArrayTest : public UnitTest 
{

public:

   void setUp(){}
   void tearDown(){}
   void testConstructor();
   void testCopyConstructor();
   void testAssignment();
   void testSubscript();
   void testIterator();
   void testSerializeMemory();
   void testSerializeFile();

};

void FArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   FArray<int,3> v;
   TEST_ASSERT(v.capacity() == 3);
} 

void FArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   FArray<int,3 > v;
   TEST_ASSERT(v.capacity() == 3 );
   v[0] = 3;
   v[1] = 4;
   v[2] = 5;

   FArray<int,3> u(v);
   TEST_ASSERT(u.capacity() == 3 );
   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );
   TEST_ASSERT(u[0] == 3 );
   TEST_ASSERT(u[1] == 4 );
   TEST_ASSERT(u[2] == 5 );
} 

void FArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);

   FArray<int,3 > v;
   TEST_ASSERT(v.capacity() == 3 );

   FArray<int,3 > u;
   TEST_ASSERT(u.capacity() == 3 );

   v[0] = 3;
   v[1] = 4;
   v[2] = 5;

   u  = v;

   TEST_ASSERT(v[0] == 3 );
   TEST_ASSERT(v[1] == 4 );
   TEST_ASSERT(v[2] == 5 );
   TEST_ASSERT(u[0] == 3 );
   TEST_ASSERT(u[1] == 4 );
   TEST_ASSERT(u[2] == 5 );
}
 
void FArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   FArray<int, 3> v;
   v[0] = 3;
   v[1] = 4;
   v[2] = 5;
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
}
 
void FArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   FArray<int, 3> v;
   v[0] = 3;
   v[1] = 4;
   v[2] = 5;

   ArrayIterator<int> it;
   v.begin(it);
   TEST_ASSERT(*it == 3 );
   ++it;
   TEST_ASSERT(it.notEnd());
   TEST_ASSERT(*it == 4 );
   ++it;
   TEST_ASSERT(it.notEnd());
   TEST_ASSERT(*it == 5 );
   ++it;
   TEST_ASSERT(it.isEnd());
} 

void FArrayTest::testSerializeMemory()
{
   printMethod(TEST_FUNC);

   const int Capacity = 3;
   typedef std::complex<double> Data;

   FArray<Data, Capacity> v;

   for (int i=0; i < Capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
   int size = memorySize(v);
  
   int i1 = 13;
   int i2;

   MemoryOArchive oArchive;
   oArchive.allocate(size + 12);

   oArchive << v;
   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   oArchive << i1;

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == 3);

   FArray<Data, Capacity> u;

   MemoryIArchive iArchive;
   iArchive = oArchive;
   TEST_ASSERT(iArchive.begin()  == oArchive.begin());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin());

   // Load into u and i2
   iArchive >> u;
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);

   iArchive >> i2;
   TEST_ASSERT(iArchive.cursor() == iArchive.end());
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);

   // Release
   iArchive.release();
   TEST_ASSERT(!iArchive.isAllocated());
   TEST_ASSERT(iArchive.begin() == 0);
   TEST_ASSERT(iArchive.cursor() == 0);
   TEST_ASSERT(iArchive.end() == 0);
   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size + sizeof(int));

   // Clear values of u and i2
   for (int i=0; i < Capacity; i++ ) {
      real(u[i]) = 0.0;
      imag(u[i]) = 0.0;
   }
   i2 = 0;

   // Reload into u and i2
   iArchive = oArchive;
   iArchive >> u;
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());
   TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);

   iArchive >> i2;
   TEST_ASSERT(iArchive.cursor() == iArchive.end());
   TEST_ASSERT(iArchive.begin() == oArchive.begin());
   TEST_ASSERT(iArchive.end() == oArchive.cursor());

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);
}


void FArrayTest::testSerializeFile()
{
   printMethod(TEST_FUNC);

   const int Capacity = 3;
   typedef std::complex<double> Data;

   FArray<Data, Capacity> v;

   for (int i=0; i < Capacity; i++ ) {
      real(v[i]) = (i+1)*10 ;
      imag(v[i]) = (i+1)*10 + 0.1;
   }
  
   int i1 = 13;
   int i2;

   BinaryFileOArchive oArchive;
   openOutputFile("binary", oArchive.file());
   oArchive << v;
   oArchive << i1;
   oArchive.file().close();

   // Show that v is unchanged by packing
   TEST_ASSERT(imag(v[0])==10.1);
   TEST_ASSERT(real(v[1])==20.0);
   TEST_ASSERT(imag(v[2])==30.1);
   TEST_ASSERT(v.capacity() == 3);

   FArray<Data, Capacity> u;

   BinaryFileIArchive iArchive;
   openInputFile("binary", iArchive.file());
   iArchive >> u;
   iArchive >> i2;
   iArchive.file().close();

   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[2]) == 30.1);
   TEST_ASSERT(i2 == 13);
   TEST_ASSERT(u.capacity() == 3);

   // Clear values of u and i2
   for (int i=0; i < Capacity; i++ ) {
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

}


TEST_BEGIN(FArrayTest)
TEST_ADD(FArrayTest, testConstructor)
TEST_ADD(FArrayTest, testCopyConstructor)
TEST_ADD(FArrayTest, testAssignment)
TEST_ADD(FArrayTest, testIterator)
TEST_ADD(FArrayTest, testSerializeMemory)
TEST_ADD(FArrayTest, testSerializeFile)
TEST_END(FArrayTest)

#endif
