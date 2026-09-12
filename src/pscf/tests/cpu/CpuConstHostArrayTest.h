#ifndef PSCF_CPU_CONST_HOST_ARRAY_TEST_H
#define PSCF_CPU_CONST_HOST_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/ConstHostArray.h>
#include <pscf/backend/cpp/DeviceArray.h>

using namespace Util;
using namespace Pscf;

class CpuConstHostArrayTest : public UnitTest
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
   void testAssignFromDevice();

};


void CpuConstHostArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      ConstHostArray<Data,CPT> v;
      TEST_ASSERT(v.size() == 0 );
      TEST_ASSERT(!v.isAssociated());
   }
}

void CpuConstHostArrayTest::testAssignFromDevice()
{
   printMethod(TEST_FUNC);
   TEST_ASSERT(Memory::total() == memory_);
   ConstHostArray<Data,CPT> u;
   {
      // Data owner
      DeviceArray<Data,CPT> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);

      // Data user
      u = v;
      TEST_ASSERT(u.size() == capacity);
      TEST_ASSERT(u.isAssociated());

      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }

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
      TEST_ASSERT(!u.isAssociated());

      v.deallocate();
      TEST_ASSERT(v.capacity() == 0);
      TEST_ASSERT(!v.isAllocated());
      TEST_ASSERT(!v.isAssociated());
      TEST_ASSERT(!v.isOwner());

   }
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(CpuConstHostArrayTest)
TEST_ADD(CpuConstHostArrayTest, testDefaultConstructor)
TEST_ADD(CpuConstHostArrayTest, testAssignFromDevice)
TEST_END(CpuConstHostArrayTest)

#endif
