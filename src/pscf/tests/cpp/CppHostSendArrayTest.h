#ifndef PSCF_CPP_HOST_ARRAY_TEST_H
#define PSCF_CPP_HOST_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/HostSendArray.h>
#include <pscf/backend/cpp/DeviceArray.h>

#include <util/containers/DArray.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>


using namespace Util;
using namespace Pscf;

class CppHostSendArrayTest : public UnitTest
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
   void testAssociate();
   void testAssociateCmplx();
   //void testAssign();
   void testBaseClassReference();
   void testDArrayOuter();

};


void CppHostSendArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      HostSendArray<Data,CPT> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
}

void CppHostSendArrayTest::testAssociate()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      DeviceArray<Data,CPT> w(capacity);
      for (int i=0; i < capacity; i++ ) {
         w[i] = (i+1)*10.0 ;
      }

      HostSendArray<Data,CPT> v;
      v.associate(w);
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v[2] == 30.0);
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostSendArrayTest::testAssociateCmplx()
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
      HostSendArray<std::complex<Data>, CPT> v;
      v.associate(w);

      // Test elements
      TEST_ASSERT(eq(v[0].real(), 10.0));
      TEST_ASSERT(eq(v[1].imag(), 20.1));
      TEST_ASSERT(eq(v[2].real(), 30.0));
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(capacity*sizeof(std::complex<Data>)));
   }
   TEST_ASSERT(Memory::total() == memory_);
}

#if 0
void CppHostSendArrayTest::testAssign()
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
      HostSendArray<Data,CPT> u;
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
#endif

void CppHostSendArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   {
      // Data owner (device array)
      DeviceArray<Data,CPT> w(3);

      HostSendArray<Data,CPT> v;
      v.associate(w);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }

      Array<Data>& u = v;
      TEST_ASSERT(u[0] == 10.0);
      TEST_ASSERT(u[2] == 30.0);
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppHostSendArrayTest::testDArrayOuter()
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
   DArray< HostSendArray<Data,CPT> > u;
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

   // Assign (dissociates for T=CPT)
   for (int i = 0; i < m; ++i) {
      v[i] = u[i];
      u[i].dissociate();
   }

   // Test equality
   for (int i = 0; i < m; ++i) {
      for (int j=0; j < capacity; j++ ) {
         TEST_ASSERT(eq(v[i][j], (j+1)*10.0 + i));
      }
   }

   // De-allocate device arrays
   for (int i = 0; i < m; ++i) {
      v[i].deallocate();
   }

}

TEST_BEGIN(CppHostSendArrayTest)
TEST_ADD(CppHostSendArrayTest, testDefaultConstructor)
TEST_ADD(CppHostSendArrayTest, testAssociate)
TEST_ADD(CppHostSendArrayTest, testAssociateCmplx)
//TEST_ADD(CppHostSendArrayTest, testAssign)
TEST_ADD(CppHostSendArrayTest, testBaseClassReference)
TEST_ADD(CppHostSendArrayTest, testDArrayOuter)

TEST_END(CppHostSendArrayTest)

#endif
