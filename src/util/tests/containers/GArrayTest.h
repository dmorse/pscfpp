#ifndef G_ARRAY_TEST_H
#define G_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/GArray.h>
#include <util/containers/Array.h>

using namespace Util;

class GArrayTest : public UnitTest 
{

   int memory_;

public:

   void setUp();
   void testReserve();
   void testConstructor();
   void testSubscript();
   void testDefaultReserve();
   void testAppendResize();
   void testResize1();
   void testResize2();
   void testCopyConstructor();
   void testSerialize1File();

};


void GArrayTest::setUp()
{  memory_ = Memory::total(); } 

void GArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      TEST_ASSERT(v.capacity() == 0);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testReserve()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      v.reserve(3);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.size() == 0 );
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void GArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      v.reserve(3);
      v.append(3);
      v.append(4);
      v.append(5);
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
      TEST_ASSERT(v[2] == 5);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testDefaultReserve()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      v.append(3);
      v.append(4);
      v.append(5);
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 64);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
      TEST_ASSERT(v[2] == 5);
      v.append(6);
      v.append(7);
      TEST_ASSERT(v.size() == 5);
      TEST_ASSERT(v.capacity() == 64);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testAppendResize()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      v.reserve(3);
      v.append(3);
      v.append(4);
      v.append(5);
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 3);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
      TEST_ASSERT(v[2] == 5);
      v.append(6);
      v.append(7);
      TEST_ASSERT(v.size() == 5);
      TEST_ASSERT(v.capacity() == 6);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testResize1()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      GArray<int> v;
      v.reserve(3);
      v.append(3);
      v.append(4);
      TEST_ASSERT(v.size() == 2);
      TEST_ASSERT(v.capacity() == 3);
      v.resize(7);
      TEST_ASSERT(v.size() == 7);
      TEST_ASSERT(v.capacity() == 12);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testResize2()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      v.resize(3);
      v[0]=3;
      v[1]=4;
      v[2]=5;
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 3);
      v.resize(8);
      TEST_ASSERT(v.size() == 8);
      TEST_ASSERT(v.capacity() == 12);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
      TEST_ASSERT(v[2] == 5);
      v[3]=6;
      v[4]=7;
      v[5]=8;
      v[6]=9;
      v[7]=10;
      v.append(11);
      TEST_ASSERT(v.size() == 9);
      TEST_ASSERT(v.capacity() == 12);
      TEST_ASSERT(v[3] == 6);
      TEST_ASSERT(v[4] == 7);
      TEST_ASSERT(v[5] == 8);
      TEST_ASSERT(v[6] == 9);
      TEST_ASSERT(v[7] == 10);
      TEST_ASSERT(v[8] == 11);
      v.resize(2);
      TEST_ASSERT(v.size() == 2);
      TEST_ASSERT(v.capacity() == 12);
      v.resize(5);
      TEST_ASSERT(v.size() == 5);
      TEST_ASSERT(v.capacity() == 12);
      TEST_ASSERT(v[2] == 0);
      TEST_ASSERT(v[3] == 0);
      TEST_ASSERT(v[4] == 0);
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void GArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      GArray<int> v;
      TEST_ASSERT(v.capacity() == 0 );
   
      v.reserve(3);
      TEST_ASSERT(v.capacity() == 3 );
      v.append(3);
      v.append(4);
      v.append(5);
   
      GArray<int> u(v);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.size() == 3);
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v[0] == 3 );
      TEST_ASSERT(v[1] == 4 );
      TEST_ASSERT(v[2] == 5 );
      TEST_ASSERT(u[0] == 3 );
      TEST_ASSERT(u[1] == 4 );
      TEST_ASSERT(u[2] == 5 );
   }
   TEST_ASSERT(Memory::total() == memory_);
}


void GArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      GArray<int> v;
      int i1 = 13;
      TEST_ASSERT(v.capacity() == 0 );
      v.reserve(4);
      TEST_ASSERT(v.capacity() == 4);
      v.append(3);
      v.append(4);
      v.append(5);
   
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 4);
      TEST_ASSERT(v[0] == 3 );
      TEST_ASSERT(v[1] == 4 );
      TEST_ASSERT(v[2] == 5 );
   
      BinaryFileOArchive oArchive;
      openOutputFile("binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 4);
      TEST_ASSERT(v[0] == 3 );
      TEST_ASSERT(v[1] == 4 );
      TEST_ASSERT(v[2] == 5 );
   
      GArray<int> u;
      int i2;
      u.reserve(4);
   
      BinaryFileIArchive iArchive;
      openInputFile("binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(u.size() == 3);
      TEST_ASSERT(u.capacity() == 4);
      TEST_ASSERT(u[0] == 3);
      TEST_ASSERT(u[1] == 4);
      TEST_ASSERT(u[2] == 5);
      TEST_ASSERT(i2 == 13);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(GArrayTest)
TEST_ADD(GArrayTest, testReserve)
TEST_ADD(GArrayTest, testConstructor)
TEST_ADD(GArrayTest, testSubscript)
TEST_ADD(GArrayTest, testDefaultReserve)
TEST_ADD(GArrayTest, testAppendResize)
TEST_ADD(GArrayTest, testResize1)
TEST_ADD(GArrayTest, testResize2)
TEST_ADD(GArrayTest, testCopyConstructor)
TEST_ADD(GArrayTest, testSerialize1File)
TEST_END(GArrayTest)

#endif
