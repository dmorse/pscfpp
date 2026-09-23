#ifndef PSCF_CPP_CONST_HOST_ARRAY_TEST_H
#define PSCF_CPP_CONST_HOST_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/ConstHostArray.h>
#include <pscf/backend/cpp/DeviceArray.h>

#include <util/containers/DArray.h>

using namespace Util;
using namespace Pscf;

class CppConstHostArrayTest : public UnitTest
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
   void testAssign();
   void testDArrayOuter();

};


void CppConstHostArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      ConstHostArray<Data,CPT> v;
      TEST_ASSERT(v.size() == 0 );
      TEST_ASSERT(!v.isAllocated());
   }
}

void CppConstHostArrayTest::testConversionConstructor()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      // Data owner
      DeviceArray<Data,CPT> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);

      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }

      // Data user
      DeviceArray<Data,CPT> const & w = v;
      ConstHostArray<Data,CPT> u(w);
      TEST_ASSERT(u.size() == capacity);
      TEST_ASSERT(u.isAllocated());

      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 30.0));

      TEST_ASSERT(eq(u[0], 10.0));
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(eq(u[2], 30.0));

      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // v.deallocate(); // Intentional error

      u.dissociate();
      TEST_ASSERT(u.size() == 0);
      TEST_ASSERT(u.cArray() == nullptr);
      TEST_ASSERT(!u.isAllocated());

      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(!v.isOwner());

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppConstHostArrayTest::testAssign()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   {
      // Data owner (device array)
      DeviceArray<Data,CPT> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }
      DeviceArray<Data,CPT> const & w = v;

      // Data user (host array u)
      ConstHostArray<Data,CPT> u;
      u = w;
      TEST_ASSERT(u.size() == capacity);
      TEST_ASSERT(u.isAllocated());

      TEST_ASSERT(eq(v[0], 10.0));
      TEST_ASSERT(eq(v[1], 20.0));
      TEST_ASSERT(eq(v[2], 30.0));

      TEST_ASSERT(eq(u[0], 10.0));
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(eq(u[2], 30.0));

      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      // v.deallocate(); // Intentional error

      u.dissociate();
      TEST_ASSERT(u.size() == 0);
      TEST_ASSERT(u.cArray() == nullptr);
      TEST_ASSERT(!u.isAllocated());

      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(!v.isOwner());

   }
   TEST_ASSERT(Memory::total() == memory_);
}

void CppConstHostArrayTest::testDArrayOuter()
{
   printMethod(TEST_FUNC);

   int m = 2;

   // Initialize array of device arrays
   DArray< DeviceArray<Data,CPT> > v;
   v.allocate(m);
   for (int i = 0; i < m; ++i) {
      v[i].allocate(capacity);
      for (int j=0; j < capacity; j++ ) {
         v[i][j] = (j+1)*10.0 + i;
      }
   }

   // Associate array of const host arrays
   DArray< ConstHostArray<Data,CPT> > u;
   u.allocate(m);
   for (int i = 0; i < m; ++i) {
      u[i] = v[i];
   }

   // Test equality
   for (int i = 0; i < m; ++i) {
      for (int j=0; j < capacity; j++ ) {
         TEST_ASSERT(eq(u[i][j],(j+1)*10.0 + i));
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

TEST_BEGIN(CppConstHostArrayTest)
TEST_ADD(CppConstHostArrayTest, testDefaultConstructor)
TEST_ADD(CppConstHostArrayTest, testConversionConstructor)
TEST_ADD(CppConstHostArrayTest, testAssign)
TEST_ADD(CppConstHostArrayTest, testDArrayOuter)
TEST_END(CppConstHostArrayTest)

#endif
