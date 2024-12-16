#ifndef PSCF_CUDA_ARRAY_TEST_H
#define PSCF_CUDA_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>
#include <pscf/cuda/GpuResources.h>
#include <util/math/Constants.h>

using namespace Util;
using namespace Pscf;

class CudaArrayTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructors()
   {
      printMethod(TEST_FUNC);
      HostDArray<cudaReal> h;
      DeviceArray<cudaReal> d;

      TEST_ASSERT(h.capacity() == 0 );
      TEST_ASSERT(!h.isAllocated() );
      TEST_ASSERT(d.capacity() == 0 );
      TEST_ASSERT(!d.isAllocated() );
   }

   void testAllocate()
   {
      printMethod(TEST_FUNC);

      HostDArray<cudaReal> h;
      DeviceArray<cudaReal> d;

      int capacity = 32;
      h.allocate(capacity);
      d.allocate(capacity);

      TEST_ASSERT(h.capacity() == capacity);
      TEST_ASSERT(h.isAllocated());
      TEST_ASSERT(d.capacity() == capacity);
      TEST_ASSERT(d.isAllocated());
      TEST_ASSERT(d.isOwner());

      h.deallocate();
      d.deallocate(); 
      TEST_ASSERT(h.capacity() == 0);
      TEST_ASSERT(!h.isAllocated());
      TEST_ASSERT(d.capacity() == 0);
      TEST_ASSERT(!d.isAllocated());
   }

   void testAssociate()
   {
      printMethod(TEST_FUNC);

      DeviceArray<cudaReal> d1;
      DeviceArray<cudaReal> d2;

      int capacity = 128;
      d1.allocate(capacity);
      TEST_ASSERT(d1.capacity() == capacity);
      TEST_ASSERT(d1.isAllocated());
      TEST_ASSERT(d1.isOwner());
      
      d2.associate(d1, capacity/4, capacity/2);
      TEST_ASSERT(d2.capacity() == capacity/2);
      TEST_ASSERT(d2.isAllocated());
      TEST_ASSERT(!d2.isOwner());
      TEST_ASSERT(d1.cArray()+(capacity/4) == d2.cArray());
   }

   void testAssignmentOperators()
   {
      printMethod(TEST_FUNC);
      
      int nx = 10;

      // Device arrays
      DeviceArray<cudaReal> d1(nx);
      DeviceArray<cudaReal> d2(nx);

      // Host arrays
      HostDArray<cudaReal> in(nx);
      HostDArray<cudaReal> out1(nx);
      HostDArray<cudaReal> out2(nx);
      HostDArray<cudaReal> out3(nx);

      // Generate data
      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         in[i] = cos(twoPi*double(i)/double(nx));
      }

      // Copy to device, then copy back to host
      d1 = in;
      out1 = d1;

      // Copy on device, then copy back to host
      d2 = d1;
      out2 = d2;

      // Copy on host
      out3 = in;

      // Check that out1, out2, and out3 all match in
      for (int i = 0; i < nx; ++i ) {
         TEST_ASSERT(eq(in[i], out1[i]));
         TEST_ASSERT(eq(in[i], out2[i]));
         TEST_ASSERT(eq(in[i], out3[i]));
      }
   }

};

TEST_BEGIN(CudaArrayTest)
TEST_ADD(CudaArrayTest, testConstructors)
TEST_ADD(CudaArrayTest, testAllocate)
TEST_ADD(CudaArrayTest, testAssociate)
TEST_ADD(CudaArrayTest, testAssignmentOperators)
TEST_END(CudaArrayTest)

#endif