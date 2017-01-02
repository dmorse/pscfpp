#ifndef FS_ARRAY_TEST_H
#define FS_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/FSArray.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryCounter.h>

using namespace Util;

class FSArrayTest : public UnitTest 
{

public:

   void setUp(){}
   void tearDown(){}
   void testConstructor();
   void testSubscript();
   void testIterator();
   void testSerialize();

};


void FSArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   FSArray<int,3> v;
   TEST_ASSERT(v.capacity() == 3);
   TEST_ASSERT(v.size() == 0);
} 

void FSArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   FSArray<int,3> v;
   TEST_ASSERT(v.size() == 0);
   v.append(3);
   v.append(4);
   v.append(5);
   TEST_ASSERT(v.size() == 3);
   TEST_ASSERT(v[0] == 3);
   TEST_ASSERT(v[1] == 4);
   TEST_ASSERT(v[2] == 5);
   v.clear();
   TEST_ASSERT(v.size() == 0);
   TEST_ASSERT(v.capacity() == 3);
}
 
void FSArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   FSArray<int, 3> v;
   v.append(3);
   v.append(4);
   v.append(5);

   ArrayIterator<int> it;
   v.begin(it);
   TEST_ASSERT(it.notEnd());
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

void FSArrayTest::testSerialize()
{
   printMethod(TEST_FUNC);

   typedef std::complex<double> Data;
   const int Capacity = 5;
   FSArray<Data, Capacity> v;

   Data a(10, 10.1), b(20,20.1), c(30,30.1);
   Data  d(40,40.1), e(50,50.1), f(60,60.1);
  
   // Fill FSArray v
   v.append(a);
   v.append(b);
   v.append(d);
   TEST_ASSERT(v.capacity() == Capacity);
   TEST_ASSERT(v.size() == 3);
 
   int size = memorySize(v);
    
   MemoryOArchive oArchive;
   oArchive.allocate(size + 2);

   // Save to archive
   oArchive << v;

   TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);

   TEST_ASSERT(real(v[0]) == 10.0);
   TEST_ASSERT(imag(v[0]) == 10.1);
   TEST_ASSERT(real(v[1]) == 20.0);
   TEST_ASSERT(imag(v[1]) == 20.1);
   TEST_ASSERT(real(v[2]) == 40.0);
   TEST_ASSERT(imag(v[2]) == 40.1);
   TEST_ASSERT(v.capacity() == Capacity);
   TEST_ASSERT(v.size() == 3);
  
   MemoryIArchive iArchive;
   iArchive = oArchive;

   // Load from archive
   FSArray<Data, Capacity> u;
   iArchive >> u;
  
   TEST_ASSERT(real(u[0]) == 10.0);
   TEST_ASSERT(imag(u[0]) == 10.1);
   TEST_ASSERT(real(u[1]) == 20.0);
   TEST_ASSERT(imag(u[1]) == 20.1);
   TEST_ASSERT(real(u[2]) == 40.0);
   TEST_ASSERT(imag(u[2]) == 40.1);
   TEST_ASSERT(u.capacity() == 5);
   TEST_ASSERT(u.size() == 3);

}

TEST_BEGIN(FSArrayTest)
TEST_ADD(FSArrayTest, testConstructor)
TEST_ADD(FSArrayTest, testSubscript)
TEST_ADD(FSArrayTest, testIterator)
TEST_ADD(FSArrayTest, testSerialize)
TEST_END(FSArrayTest)

#endif
