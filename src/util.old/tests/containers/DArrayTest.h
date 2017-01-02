#ifndef DARRAY_TEST_H
#define DARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;

class DArrayTest : public UnitTest 
{
private:

   const static int capacity = 3;

   typedef std::complex<double> Data;

   int memory_;


public:

   void setUp() 
   {  memory_ = Memory::total(); }

   void tearDown() {}
   void testConstructor();
   void testAllocate();
   void testSubscript();
   void testIterator();
   void testCopyConstructor();
   void testAssignment();
   void testBaseClassReference();
   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();

};


void DArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      DArray<Data> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT((int)Memory::total() == capacity*sizeof(Data));
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      DArray<Data> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
   
      TEST_ASSERT(real(v[0]) == 10 );
      TEST_ASSERT(imag(v[1]) == 20.1 );
      TEST_ASSERT(real(v[2]) == 30 );
      TEST_ASSERT((int)Memory::total() == capacity*sizeof(Data));
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT((int)Memory::total() == 0);
   {
      DArray<Data> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
   
      ArrayIterator<Data> it;
      v.begin(it);
      TEST_ASSERT(imag(*it) == 10.1 );
      TEST_ASSERT(!it.isEnd());
      TEST_ASSERT(it.notEnd());
      ++it;
      TEST_ASSERT(real(*it) == 20 );
      TEST_ASSERT(!it.isEnd());
      TEST_ASSERT(it.notEnd());
      ++it;
      TEST_ASSERT(imag(*it) == 30.1 );
      ++it;
      TEST_ASSERT(it.isEnd());
      TEST_ASSERT(!it.notEnd());
      TEST_ASSERT((int)Memory::total() == capacity*sizeof(Data));
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated() );
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
   
      DArray<Data> u(v);
      TEST_ASSERT(u.capacity() == capacity );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(real(v[0]) == 10 );
      TEST_ASSERT(imag(v[1]) == 20.1 );
      TEST_ASSERT(real(v[2]) == 30 );
      TEST_ASSERT(real(u[0]) == 10 );
      TEST_ASSERT(imag(u[1]) == 20.1 );
      TEST_ASSERT(real(u[2]) == 30 );
      TEST_ASSERT((int)Memory::total() == 2*capacity*sizeof(Data));
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);

   {
      DArray<Data> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.isAllocated() );
   
      DArray<Data> u;
      u.allocate(3);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
   
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
   
      u  = v;
   
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(real(v[0]) == 10 );
      TEST_ASSERT(imag(v[1]) == 20.1 );
      TEST_ASSERT(real(v[2]) == 30 );
      TEST_ASSERT(real(u[0]) == 10 );
      TEST_ASSERT(imag(u[1]) == 20.1 );
      TEST_ASSERT(real(u[2]) == 30 );
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void DArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
      
      Array<Data>& u = v;
      TEST_ASSERT(real(u[0]) == 10 );
      TEST_ASSERT(imag(u[1]) == 20.1 );
      TEST_ASSERT(real(u[2]) == 30 );
   }
   TEST_ASSERT(Memory::total() == memory_);
}


void DArrayTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
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
   
      DArray<Data> u;
      u.allocate(3);
   
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
      for (int i=0; i < capacity; i++ ) {
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
   TEST_ASSERT(Memory::total() == memory_);

}

void DArrayTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         real(v[i]) = (i+1)*10 ;
         imag(v[i]) = (i+1)*10 + 0.1;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(imag(v[0])==10.1);
      TEST_ASSERT(real(v[1])==20.0);
      TEST_ASSERT(imag(v[2])==30.1);
      TEST_ASSERT(v.capacity() == capacity);
   
      DArray<Data> u;
   
      // Note: We do not allocate DArray<Data> u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(imag(u[0]) == 10.1);
      TEST_ASSERT(real(u[1]) == 20.0);
      TEST_ASSERT(imag(u[2]) == 30.1);
      TEST_ASSERT(u.capacity() == 3);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void DArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
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
   
      DArray<Data> u;
      u.allocate(3);
   
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
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void DArrayTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      DArray<Data> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
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
   
      DArray<Data> u;
   
      // u.allocate(3); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
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
   }
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(DArrayTest)
TEST_ADD(DArrayTest, testConstructor)
TEST_ADD(DArrayTest, testAllocate)
TEST_ADD(DArrayTest, testSubscript)
TEST_ADD(DArrayTest, testIterator)
TEST_ADD(DArrayTest, testCopyConstructor)
TEST_ADD(DArrayTest, testAssignment)
TEST_ADD(DArrayTest, testBaseClassReference)
TEST_ADD(DArrayTest, testSerialize1Memory)
TEST_ADD(DArrayTest, testSerialize2Memory)
TEST_ADD(DArrayTest, testSerialize1File)
TEST_ADD(DArrayTest, testSerialize2File)
TEST_END(DArrayTest)

#endif
