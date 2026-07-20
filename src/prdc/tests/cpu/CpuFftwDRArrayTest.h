#ifndef PRDC_CPU_FFTW_DRARRAY_TEST_H
#define PRDC_CPU_FFTW_DRARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cpu/FftwDRArray.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>


using namespace Util;
using namespace Pscf;
using namespace Prdc;

class CpuFftwDRArrayTest : public UnitTest
{
private:

   const static int capacity = 3;

   typedef double Data;

   long int memory_;

public:

   void setUp()
   {  memory_ = Memory::total(); }

   void tearDown() {}
   void testDefaultConstructor();
   void testAllocateConstructor();
   void testAllocate();
   void testSubscript();
   void testSubscriptCmplx();
   void testAssociate();
   void testIterator();
   void testCopyConstructor();
   void testCopyConstructorCmplx();
   void testAssignment();
   void testAssignmentCmplx();
   void testBaseClassReference();

   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();

};


void CpuFftwDRArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<Data> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
      TEST_ASSERT(!v.isOwner());
      TEST_ASSERT(!v.isAssociated());
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testAllocateConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      FftwDRArray<Data> v(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());
      TEST_ASSERT(!v.isAssociated());
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // Deallocate array
      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(Memory::total() == memory_);

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == 0);
   {
      FftwDRArray<Data> v;

      // Allocate array
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());
      TEST_ASSERT(!v.isAssociated());
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // Deallocate array
      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT((int)Memory::total() == memory_);

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      FftwDRArray<Data> v(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }

      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v[2] == 30.0);
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testAssociate()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   FftwDRArray<Data> u;
   {
      // Data owner
      FftwDRArray<Data> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);

      // Data user
      u.associate(v, 1, capacity - 1);
      TEST_ASSERT(u.capacity() == capacity - 1);
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u.isAssociated());
      TEST_ASSERT(!u.isOwner());

      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }

      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 30.0));
      TEST_ASSERT(eq(u[0], 20.0));
      TEST_ASSERT(eq(u[1], 30.0));
      u[1] = 25.0;
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 25.0));
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // v.deallocate(); // Intentional error

      u.dissociate();
      TEST_ASSERT(u.capacity() == 0);
      TEST_ASSERT(!u.isAllocated());
      TEST_ASSERT(!u.isAssociated());
      TEST_ASSERT(!u.isOwner());

      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(!v.isAssociated());
      TEST_ASSERT(!v.isOwner());

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testSubscriptCmplx()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      FftwDRArray< std::complex<Data> > v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i].real((i+1)*10.0);
         v[i].imag((i+1)*10.0 + 0.1);
      }

      TEST_ASSERT(eq(v[0].real(), 10.0));
      TEST_ASSERT(eq(v[1].imag(), 20.1));
      TEST_ASSERT(eq(v[2].real(), 30.0));
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(capacity*sizeof(std::complex<Data>)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      // Data owner
      FftwDRArray<Data> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());
      TEST_ASSERT(!v.isAssociated());
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      FftwDRArray<Data> u(v);
      TEST_ASSERT(u.capacity() == capacity);
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u.isOwner());
      TEST_ASSERT(!u.isAssociated());

      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 30.0));
      TEST_ASSERT(eq(u[0], 10.0));
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(eq(u[2], 30.0));
      u[1] = 25.0;
      TEST_ASSERT(eq(u[1], 25.0));
      TEST_ASSERT(eq(v[1], 20.0));
      tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + 2*capacity*sizeof(Data)));

      u.deallocate();
      tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));
   }
   TEST_ASSERT(Memory::total() == (long int)memory_);
}

void CpuFftwDRArrayTest::testCopyConstructorCmplx()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray< std::complex<Data> > v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );

      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated() );
      for (int i=0; i < capacity; i++ ) {
         v[i].real((i+1)*10.0);
         v[i].imag((i+1)*10.0 + 0.1);
      }

      FftwDRArray< std::complex<Data> > u(v);
      TEST_ASSERT(u.capacity() == capacity);
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(u.isOwner());
      TEST_ASSERT(!u.isAssociated());
      TEST_ASSERT(real(v[0]) == 10.0 );
      TEST_ASSERT(imag(v[1]) == 20.1 );
      TEST_ASSERT(real(v[2]) == 30.0 );
      TEST_ASSERT(real(u[0]) == 10.0 );
      TEST_ASSERT(imag(u[1]) == 20.1 );
      TEST_ASSERT(real(u[2]) == 30 );
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(2*capacity*sizeof(std::complex<Data>)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testAssignment()
{
   printMethod(TEST_FUNC);

   {
      FftwDRArray<Data> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.isAllocated() );
      TEST_ASSERT(v.isOwner() );
      TEST_ASSERT(!v.isAssociated() );

      FftwDRArray<Data> u;
      u.allocate(3);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(u.isOwner() );
      TEST_ASSERT(!u.isAssociated() );

      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10;
      }

      u = v;

      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(u.isOwner() );
      TEST_ASSERT(!u.isAssociated() );
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[2] == 30.0);
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testAssignmentCmplx()
{
   printMethod(TEST_FUNC);

   {
      FftwDRArray< std::complex<Data> > v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3);
      TEST_ASSERT(v.isAllocated());

      FftwDRArray< std::complex<Data> > u;
      u.allocate(3);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );

      for (int i=0; i < capacity; i++ ) {
         v[i].real((i+1)*10.0);
         v[i].imag((i+1)*10.0 + 0.1);
      }

      u  = v;

      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(real(v[0]) == 10.0);
      TEST_ASSERT(imag(v[1]) == 20.1);
      TEST_ASSERT(real(v[2]) == 30.0);
      TEST_ASSERT(real(u[0]) == 10.0);
      TEST_ASSERT(imag(u[1]) == 20.1);
      TEST_ASSERT(real(u[2]) == 30.0);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT((int)Memory::total() == 0);
   {
      FftwDRArray<Data> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }

      ArrayIterator<Data> it;
      v.begin(it);
      TEST_ASSERT(eq(*it, 10.0));
      TEST_ASSERT(!it.isEnd());
      TEST_ASSERT(it.notEnd());
      ++it;
      TEST_ASSERT(eq(*it, 20.0));
      TEST_ASSERT(!it.isEnd());
      TEST_ASSERT(it.notEnd());
      ++it;
      TEST_ASSERT(eq(*it, 30.0));
      ++it;
      TEST_ASSERT(it.isEnd());
      TEST_ASSERT(!it.notEnd());
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(capacity * sizeof(Data)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<Data> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }

      Array<Data>& u = v;
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CpuFftwDRArrayTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
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
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      FftwDRArray<double> u;
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
   
      TEST_ASSERT(u[1] == 20.0);
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
         u[i] = 0.0;
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
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }

}

void CpuFftwDRArrayTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<double> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      FftwDRArray<double> u;
   
      // Note: We do not allocate FftwDRArray<double> u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void CpuFftwDRArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;

      BinaryFileOArchive oArchive;
      openOutputFile("out/FftwDRArray.arx", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      FftwDRArray<double> u;
      u.allocate(3);
   
      BinaryFileIArchive iArchive;
      openInputFile("out/FftwDRArray.arx", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/FftwDRArray.arx", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void CpuFftwDRArrayTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      FftwDRArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;
  
      BinaryFileOArchive oArchive;
      openOutputFile("out/FftwDRArray.arx", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      FftwDRArray<double> u;
   
      // u.allocate(3); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
      BinaryFileIArchive iArchive;
      openInputFile("out/FftwDRArray.arx", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/FftwDRArray.arx", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}

TEST_BEGIN(CpuFftwDRArrayTest)
TEST_ADD(CpuFftwDRArrayTest, testDefaultConstructor)
TEST_ADD(CpuFftwDRArrayTest, testAllocateConstructor)
TEST_ADD(CpuFftwDRArrayTest, testAllocate)
TEST_ADD(CpuFftwDRArrayTest, testSubscript)
TEST_ADD(CpuFftwDRArrayTest, testSubscriptCmplx)
TEST_ADD(CpuFftwDRArrayTest, testAssociate)
TEST_ADD(CpuFftwDRArrayTest, testCopyConstructor)
TEST_ADD(CpuFftwDRArrayTest, testCopyConstructorCmplx)
TEST_ADD(CpuFftwDRArrayTest, testAssignment)
TEST_ADD(CpuFftwDRArrayTest, testAssignmentCmplx)
TEST_ADD(CpuFftwDRArrayTest, testIterator)
TEST_ADD(CpuFftwDRArrayTest, testBaseClassReference)

TEST_ADD(CpuFftwDRArrayTest, testSerialize1Memory)
TEST_ADD(CpuFftwDRArrayTest, testSerialize2Memory)
TEST_ADD(CpuFftwDRArrayTest, testSerialize1File)
TEST_ADD(CpuFftwDRArrayTest, testSerialize2File)

TEST_END(CpuFftwDRArrayTest)

#endif
