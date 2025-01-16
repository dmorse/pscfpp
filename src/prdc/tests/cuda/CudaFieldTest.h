#ifndef PRDC_CUDA_FIELD_TEST_H
#define PRDC_CUDA_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/RField.h>
#include <prdc/cuda/CField.h>
#include <prdc/cuda/RFieldDft.h>

#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/GpuResources.h>
#include <pscf/math/IntVec.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaFieldTest : public UnitTest
{

public:

   void setUp() {};

   void tearDown() {};

   void testConstructors();
   void testAllocate();
   void testRFieldRoundTrip();
   void testCFieldRoundTrip();
   void testRFieldDftRoundTrip();

};

void CudaFieldTest::testConstructors()
{
   printMethod(TEST_FUNC);
   {
      Prdc::Cuda::RField<1> r;
      Prdc::Cuda::CField<1> c;
      Prdc::Cuda::RFieldDft<1> d;
      TEST_ASSERT(r.capacity() == 0 );
      TEST_ASSERT(!r.isAllocated() );
      TEST_ASSERT(c.capacity() == 0 );
      TEST_ASSERT(!c.isAllocated() );
      TEST_ASSERT(d.capacity() == 0 );
      TEST_ASSERT(!d.isAllocated() );
   }
}

void CudaFieldTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      Prdc::Cuda::RField<2> r;
      Prdc::Cuda::CField<2> c;
      Prdc::Cuda::RFieldDft<2> d;
      IntVec<2> meshDims, dftDims;
      meshDims[0] = 5;
      meshDims[1] = 7;
      int capacity = meshDims[0] * meshDims[1];
      dftDims = meshDims;
      dftDims[1] = dftDims[1]/2 + 1;
      int dftCapacity = dftDims[0] * dftDims[1];
      r.allocate(meshDims);
      c.allocate(meshDims);
      d.allocate(meshDims);
      
      TEST_ASSERT(r.capacity() == capacity);
      TEST_ASSERT(r.meshDimensions() == meshDims);
      TEST_ASSERT(r.isAllocated());
      TEST_ASSERT(c.capacity() == capacity);
      TEST_ASSERT(c.meshDimensions() == meshDims);
      TEST_ASSERT(c.isAllocated());
      TEST_ASSERT(d.capacity() == dftCapacity);
      TEST_ASSERT(d.meshDimensions() == meshDims);
      TEST_ASSERT(d.dftDimensions() == dftDims);
      TEST_ASSERT(d.isAllocated());
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
      HostDArray<cudaReal> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i] = (i+1)*10.0 ;
      }

      // Copy host field vh1 to device field vd
      Prdc::Cuda::RField<3> vd;
      vd.allocate(d);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostDArray<cudaReal> vh2;
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
      HostDArray<cudaComplex> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      Prdc::Cuda::CField<3> vd;
      vd.allocate(d);
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostDArray<cudaComplex> vh2;
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
      Prdc::Cuda::RFieldDft<3> vd;
      vd.allocate(d);
      int capacity = vd.capacity();

      // Initialize host (cpu) complex data in field vh1
      HostDArray<cudaComplex> vh1;
      vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostDArray<cudaComplex> vh2;
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
TEST_ADD(CudaFieldTest, testConstructors)
TEST_ADD(CudaFieldTest, testAllocate)
TEST_ADD(CudaFieldTest, testRFieldRoundTrip)
TEST_ADD(CudaFieldTest, testCFieldRoundTrip)
TEST_ADD(CudaFieldTest, testRFieldDftRoundTrip)
TEST_END(CudaFieldTest)

#endif
