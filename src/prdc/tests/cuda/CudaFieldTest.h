#ifndef PRDC_CUDA_FIELD_TEST_H
#define PRDC_CUDA_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/Field.h>
#include <prdc/cuda/Field.tpp>
#include <prdc/cpu/Field.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaFieldTest : public UnitTest
{

public:

   void setUp() {};

   void tearDown(){};

   void testFieldConstructor();
   void testFieldAllocate();
   void testFieldDoubleRoundTrip();
   void testRFieldRoundTrip();

};

void CudaFieldTest::testFieldConstructor()
{
   printMethod(TEST_FUNC);
   {
      Cuda::Field<double> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
}

void CudaFieldTest::testFieldAllocate()
{
   printMethod(TEST_FUNC);
   {
      Cuda::Field<double> v;
      int capacity = 3;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
   }
}

void CudaFieldTest::testFieldDoubleRoundTrip()
{
   printMethod(TEST_FUNC);
   {
      int capacity = 3;

      // Initialize host (cpu) data in field vh1
      Cuda::HostField<double> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i] = (i+1)*10.0 ;
      }

      // Copy host field vh1 to device field vd
      Cuda::Field<double> vd;
      vd.allocate(capacity);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      Cuda::HostField<double> vh2;
      vh2 = vd;

      TEST_ASSERT(vh2.capacity() == vh1.capacity());
      TEST_ASSERT(vh2[0] == 10.0);
      TEST_ASSERT(vh2[2] == 30.0);
      for (int i=0; i < capacity; i++ ) {
         TEST_ASSERT(eq(vh2[i], (i+1)*10.0)) ;
         TEST_ASSERT(eq(vh2[i], vh1[i]));
      }

   }
}

#if 1
void CudaFieldTest::testRFieldRoundTrip()
{
   printMethod(TEST_FUNC); {

      IntVec<3> d;
      d[0] = 2;
      d[1] = 1;
      d[2] = 3;
      int capacity = 6;

      // Initialize host (cpu) data in field vh1
      Cuda::HostField<cudaReal> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i] = (i+1)*10.0 ;
      }

      // Copy host field vh1 to device field vd
      Cuda::RField<3> vd;
      vd.allocate(d);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      Cuda::HostField<cudaReal> vh2;
      vh2 = vd;

      TEST_ASSERT(vh2.capacity() == vh1.capacity());
      TEST_ASSERT(vh2[0] == 10.0);
      TEST_ASSERT(vh2[2] == 30.0);
      for (int i = 0; i < capacity; i++ ) {
         TEST_ASSERT(eq(vh2[i], (i+1)*10.0)) ;
      }

   }
}
#endif

TEST_BEGIN(CudaFieldTest)
TEST_ADD(CudaFieldTest, testFieldConstructor)
TEST_ADD(CudaFieldTest, testFieldAllocate)
TEST_ADD(CudaFieldTest, testFieldDoubleRoundTrip)
TEST_ADD(CudaFieldTest, testRFieldRoundTrip)
TEST_END(CudaFieldTest)

#endif
