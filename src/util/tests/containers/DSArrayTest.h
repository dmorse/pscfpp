#ifndef DS_ARRAY_TEST_H
#define DS_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DSArray.h>
#include <util/containers/Array.h>

using namespace Util;

class DSArrayTest : public UnitTest 
{

   int memory_;

public:

   void setUp(){ memory_ = Memory::total(); }
   void tearDown(){};
   void testAllocate();
   void testConstructor();
   void testSubscript();
   void testCopyConstructor();
   void testAssignment();
   void testSerialize1File();

};


void DSArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      DSArray<int> v;
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated() );
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DSArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      DSArray<int> v;
      v.allocate(3);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.isAllocated() );
      TEST_ASSERT(v.size() == 0 );
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void DSArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      DSArray<int> v;
      v.allocate(3);
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

void DSArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   {
      DSArray<int> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   
      v.allocate(3);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.isAllocated() );
      v.append(3);
      v.append(4);
      v.append(5);
   
      DSArray<int> u(v);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
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

void DSArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);
   {
      DSArray<int> v;
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
   
      v.allocate(5);
      TEST_ASSERT(v.capacity() == 5);
      TEST_ASSERT(v.isAllocated());
      v.append(3);
      v.append(4);
      v.append(5);
   
      DSArray<int> u;
      u = v;
      TEST_ASSERT(u.capacity() == 5);
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u.size() == 3);
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v[0] == 3);
      TEST_ASSERT(v[1] == 4);
      TEST_ASSERT(v[2] == 5);
      TEST_ASSERT(u[0] == 3);
      TEST_ASSERT(u[1] == 4);
      TEST_ASSERT(u[2] == 5);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void DSArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);

   {
      DSArray<int> v;
      int i1 = 13;
      
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   
      v.allocate(4);
      TEST_ASSERT(v.capacity() == 4);
      TEST_ASSERT(v.isAllocated());
      v.append(3);
      v.append(4);
      v.append(5);
   
      TEST_ASSERT(v.size() == 3);
      TEST_ASSERT(v.capacity() == 4);
      TEST_ASSERT(v.isAllocated());
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
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v[0] == 3 );
      TEST_ASSERT(v[1] == 4 );
      TEST_ASSERT(v[2] == 5 );
   
      DSArray<int> u;
      int i2;
      u.allocate(3);
   
      BinaryFileIArchive iArchive;
      openInputFile("binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(u.size() == 3);
      TEST_ASSERT(u.capacity() == 3);
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u[0] == 3);
      TEST_ASSERT(u[1] == 4);
      TEST_ASSERT(u[2] == 5);
      TEST_ASSERT(i2 == 13);
   
      #if 0
      // Clear values of u and i2
   
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
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(DSArrayTest)
TEST_ADD(DSArrayTest, testAllocate)
TEST_ADD(DSArrayTest, testConstructor)
TEST_ADD(DSArrayTest, testSubscript)
TEST_ADD(DSArrayTest, testCopyConstructor)
TEST_ADD(DSArrayTest, testAssignment)
TEST_ADD(DSArrayTest, testSerialize1File)
TEST_END(DSArrayTest)

#endif
