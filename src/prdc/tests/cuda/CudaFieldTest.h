#ifndef PRDC_CUDA_FIELD_TEST_H
#define PRDC_CUDA_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/Field.h>
#include <prdc/cuda/Field.tpp>
#include <prdc/cuda/RField.tpp>
#include <prdc/cuda/CField.tpp>
#include <prdc/cuda/RFieldDft.tpp>
#include <prdc/cuda/HostField.tpp>

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
   void testCFieldRoundTrip();
   void testRFieldDftRoundTrip();

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

void CudaFieldTest::testCFieldRoundTrip()
{
   printMethod(TEST_FUNC); {

      IntVec<3> d;
      d[0] = 2;
      d[1] = 1;
      d[2] = 3;
      int capacity = 6;

      // Initialize host (cpu) data in field vh1
      Cuda::HostField<cudaComplex> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      Cuda::CField<3> vd;
      vd.allocate(d);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      Cuda::HostField<cudaComplex> vh2;
      vh2 = vd;

      TEST_ASSERT(vh2.capacity() == vh1.capacity());
      TEST_ASSERT(vh2[0].x == 10.0);
      TEST_ASSERT(vh2[0].y ==  3.0);
      TEST_ASSERT(vh2[2].x == 30.0);
      TEST_ASSERT(vh2[2].y == 23.0);
      for (int i = 0; i < capacity; i++ ) {
         TEST_ASSERT(eq(vh2[i].x, vh1[i].x)) ;
         TEST_ASSERT(eq(vh2[i].y, vh1[i].y)) ;
      }

   }
}

void CudaFieldTest::testRFieldDftRoundTrip()
{
   printMethod(TEST_FUNC); {

      IntVec<3> d;
      d[0] = 2;
      d[1] = 1;
      d[2] = 3;

      // Allocate memory on device
      Cuda::RFieldDft<3> vd;
      vd.allocate(d);
      int capacity = vd.capacity();

      // Initialize host (cpu) complex data in field vh1
      Cuda::HostField<cudaComplex> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      Cuda::HostField<cudaComplex> vh2;
      vh2 = vd;

      TEST_ASSERT(vh2.capacity() == vh1.capacity());
      TEST_ASSERT(vh2[0].x == 10.0);
      TEST_ASSERT(vh2[0].y ==  3.0);
      TEST_ASSERT(vh2[2].x == 30.0);
      TEST_ASSERT(vh2[2].y == 23.0);
      for (int i = 0; i < capacity; i++ ) {
         TEST_ASSERT(eq(vh2[i].x, vh1[i].x)) ;
         TEST_ASSERT(eq(vh2[i].y, vh1[i].y)) ;
      }

   }
}

TEST_BEGIN(CudaFieldTest)
TEST_ADD(CudaFieldTest, testFieldConstructor)
TEST_ADD(CudaFieldTest, testFieldAllocate)
TEST_ADD(CudaFieldTest, testFieldDoubleRoundTrip)
TEST_ADD(CudaFieldTest, testRFieldRoundTrip)
TEST_ADD(CudaFieldTest, testCFieldRoundTrip)
TEST_ADD(CudaFieldTest, testRFieldDftRoundTrip)
TEST_END(CudaFieldTest)

#endif
