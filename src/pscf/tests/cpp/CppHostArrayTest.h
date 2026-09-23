#ifndef PSCF_CPP_HOST_ARRAY_TEST_H
#define PSCF_CPP_HOST_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/HostArray.h>
#include <pscf/backend/cpp/DeviceArray.h>

#include <util/containers/DArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>


using namespace Util;
using namespace Pscf;

class CppHostArrayTest : public UnitTest
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
   void testConversionConstructor();
   void testConversionConstructorCmplx();
   void testAssign();
   void testBaseClassReference();
   void testDArrayOuter();

};


void CppHostArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      HostArray<Data,CPT> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
}

void CppHostArrayTest::testConversionConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      DeviceArray<Data,CPT> w(capacity);
      for (int i=0; i < capacity; i++ ) {
         w[i] = (i+1)*10.0 ;
      }

      HostArray<Data,CPT> v(w);
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v[2] == 30.0);
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostArrayTest::testConversionConstructorCmplx()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      // Initialize device array
      DeviceArray< std::complex<Data>, CPT> w;
      w.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         w[i].real((i+1)*10.0);
         w[i].imag((i+1)*10.0 + 0.1);
      }

      // Copy construct
      HostArray<std::complex<Data>, CPT> v(w);

      // Test elements
      TEST_ASSERT(eq(v[0].real(), 10.0));
      TEST_ASSERT(eq(v[1].imag(), 20.1));
      TEST_ASSERT(eq(v[2].real(), 30.0));
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(capacity*sizeof(std::complex<Data>)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostArrayTest::testAssign()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      // Data owner (device array)
      DeviceArray<Data,CPT> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());

      // Set data in device array
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }

      // Data user (host array)
      HostArray<Data,CPT> u;
      TEST_ASSERT(u.capacity() == 0);
      TEST_ASSERT(!u.isAllocated());

      // Copy device -> host (without prior association)
      u = v;
      TEST_ASSERT(u.capacity() == capacity);
      TEST_ASSERT(u.isAllocated());

      // Test equality after assignment
      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 30.0));
      TEST_ASSERT(eq(u[0], 10.0));
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(eq(u[2], 30.0));

      // Modify on host
      u[1] = 25.0;
      TEST_ASSERT(eq(u[1], 25.0));
      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 25.0));
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // v.deallocate(); // Intentional error

      u.dissociate();
      TEST_ASSERT(u.capacity() == 0);
      TEST_ASSERT(!u.isAllocated());
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());

      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(!v.isOwner());

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   {
      DeviceArray<Data,CPT> w;
      w.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         w[i] = (i+1)*10.0;
      }
      HostArray<Data,CPT> v(w);
      

      Array<Data>& u = v;
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostArrayTest::testDArrayOuter()
{
   printMethod(TEST_FUNC);

   int m = 2;

   // Allocate array of device arrays
   DArray< DeviceArray<Data,CPT> > v;
   v.allocate(m);
   for (int i = 0; i < m; ++i) {
      v[i].allocate(capacity);
   }

   // Allocate and associate array of host arrays
   DArray< HostArray<Data,CPT> > u;
   u.allocate(m);
   for (int i = 0; i < m; ++i) {
      u[i].associate(v[i]);
   }

   // Initialize data on host
   for (int i = 0; i < m; ++i) {
      for (int j=0; j < capacity; j++ ) {
         u[i][j] = (j+1)*10.0 + i;
      }
   }

   // Assign (does nothing for T=CPT)
   for (int i = 0; i < m; ++i) {
      v[i] = u[i];
   }

   // Test equality
   for (int i = 0; i < m; ++i) {
      for (int j=0; j < capacity; j++ ) {
         TEST_ASSERT(eq(v[i][j], (j+1)*10.0 + i));
         TEST_ASSERT(eq(u[i][j], v[i][j]));
      }
   }

   // Dissociate host arrays
   for (int i = 0; i < m; ++i) {
      u[i].dissociate();
   }

   // De-allocate device arrays
   for (int i = 0; i < m; ++i) {
      v[i].deallocate();
   }

}

TEST_BEGIN(CppHostArrayTest)
TEST_ADD(CppHostArrayTest, testDefaultConstructor)
TEST_ADD(CppHostArrayTest, testConversionConstructor)
TEST_ADD(CppHostArrayTest, testConversionConstructorCmplx)
TEST_ADD(CppHostArrayTest, testAssign)
TEST_ADD(CppHostArrayTest, testBaseClassReference)
TEST_ADD(CppHostArrayTest, testDArrayOuter)

TEST_END(CppHostArrayTest)

#endif
