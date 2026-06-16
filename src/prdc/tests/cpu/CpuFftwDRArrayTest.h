#ifndef PRDC_CPU_FFTW_DRARRAY_TEST_H
#define PRDC_CPU_FFTW_DRARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/FftwDRArray.h>

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
};


void CpuFftwDRArrayTest::testDefaultConstructor()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDRArray<Data> v;
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
      Cpu::FftwDRArray<Data> v(capacity);
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
      Cpu::FftwDRArray<Data> v;

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
      Cpu::FftwDRArray<Data> v(capacity);
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
   Cpu::FftwDRArray<Data> u;
   {
      // Data owner
      Cpu::FftwDRArray<Data> v(capacity);
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
      Cpu::FftwDRArray< std::complex<Data> > v;
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
      Cpu::FftwDRArray<Data> v(capacity);
      TEST_ASSERT(v.capacity() == capacity);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.isOwner());
      TEST_ASSERT(!v.isAssociated());
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }
      long int tot = Memory::total();
      TEST_ASSERT(tot == (long int)(memory_ + capacity*sizeof(Data)));

      Cpu::FftwDRArray<Data> u(v);
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
      Cpu::FftwDRArray< std::complex<Data> > v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );

      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated() );
      for (int i=0; i < capacity; i++ ) {
         v[i].real((i+1)*10.0);
         v[i].imag((i+1)*10.0 + 0.1);
      }

      Cpu::FftwDRArray< std::complex<Data> > u(v);
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
      Cpu::FftwDRArray<Data> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3 );
      TEST_ASSERT(v.isAllocated() );
      TEST_ASSERT(v.isOwner() );
      TEST_ASSERT(!v.isAssociated() );

      Cpu::FftwDRArray<Data> u;
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
      Cpu::FftwDRArray< std::complex<Data> > v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3);
      TEST_ASSERT(v.isAllocated());

      Cpu::FftwDRArray< std::complex<Data> > u;
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
      Cpu::FftwDRArray<Data> v;
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
      Cpu::FftwDRArray<Data> v;
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
TEST_END(CpuFftwDRArrayTest)

#endif
