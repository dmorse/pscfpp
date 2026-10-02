#ifndef PSCF_CUDA_ARRAY_TEST_H
#define PSCF_CUDA_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cuda/DeviceArray.h>
#include <pscf/backend/cuda/HostArray.h>
#include <pscf/backend/cuda/ConstHostArray.h>
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
      HostArray<double,CUT> h;
      ConstHostArray<double,CUT> c;
      DeviceArray<double,CUT> d;

      TEST_ASSERT(h.capacity() == 0 );
      TEST_ASSERT(!h.isAllocated() );
      TEST_ASSERT(c.capacity() == 0 );
      TEST_ASSERT(!c.isAllocated() );
      TEST_ASSERT(d.capacity() == 0 );
      TEST_ASSERT(!d.isAllocated() );
   }

   void testAllocate()
   {
      printMethod(TEST_FUNC);

      HostArray<double,CUT> h;
      HostArray<double,CUT> c;
      DeviceArray<double,CUT> d;

      int capacity = 32;
      h.allocate(capacity);
      c.allocate(capacity);
      d.allocate(capacity);

      TEST_ASSERT(h.capacity() == capacity);
      TEST_ASSERT(h.isAllocated());
      TEST_ASSERT(c.capacity() == capacity);
      TEST_ASSERT(c.isAllocated());
      TEST_ASSERT(d.capacity() == capacity);
      TEST_ASSERT(d.isAllocated());
      TEST_ASSERT(d.isOwner());

      h.deallocate();
      c.deallocate();
      d.deallocate(); 
      TEST_ASSERT(h.capacity() == 0);
      TEST_ASSERT(!h.isAllocated());
      TEST_ASSERT(c.capacity() == 0);
      TEST_ASSERT(!c.isAllocated());
      TEST_ASSERT(d.capacity() == 0);
      TEST_ASSERT(!d.isAllocated());
   }

   void testAssociate()
   {
      printMethod(TEST_FUNC);

      DeviceArray<double,CUT> d1;
      DeviceArray<double,CUT> d2;

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
      DeviceArray<double,CUT> d1(nx);
      DeviceArray<double,CUT> d2(nx);

      // Input host arrays
      HostArray<double,CUT> in;
      in.associate(d1);

      // Generate data
      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         in[i] = cos(twoPi*double(i)/double(nx));
      }

      // Output host arrays
      HostArray<double,CUT> host1(nx);
      HostArray<double,CUT> host2(nx);

      // Copy to device, then copy back to host
      d1 = in;
      host1 = d1;

      // Copy from device -> device, then copy back to host
      d2 = d1;
      host2 = d2;


      // Check that host1 and host2 match
      for (int i = 0; i < nx; ++i ) {
         TEST_ASSERT(eq(in[i], host1[i]));
         TEST_ASSERT(eq(in[i], host2[i]));
      }

      // Copy a slice of d1, check that it is correct
      HostArray<double,CUT> host4(nx/2);
      host4.copySlice(d1, 3);
      for (int i = 0; i < nx/2; ++i ) {
         TEST_ASSERT(eq(in[i+3], host4[i]));
      }

      // Copy from d1 to a ConstHostArray
      ConstHostArray<double,CUT> host5;
      //host5.associate(d1);
      host5 = d1;
      for (int i = 0; i < nx; ++i ) {
         TEST_ASSERT(eq(in[i], host5[i]));
      }

   }

   void testConstAssignmentOperators()
   {
      printMethod(TEST_FUNC);
      
      int nx = 10;

      // Device arrays
      DeviceArray<double,CUT> d1(nx);
      DeviceArray<double,CUT> d2(nx);

      // Input array
      HostArray<double,CUT> in;
      in.associate(d1);

      // Generate data
      double twoPi = 2.0*Constants::Pi;
      for (int i=0; i < nx; ++i) {
         in[i] = cos(twoPi*double(i)/double(nx));
      }

      // Copy host -> device, then copy back to host
      d1 = in;

      // Host arrays
      ConstHostArray<double,CUT> out1(nx);
      ConstHostArray<double,CUT> out2;
      ConstHostArray<double,CUT> out3(nx);
      ConstHostArray<double,CUT> out4(nx/2);

      // Copy back to host
      out1 = d1;

      // Copy device -> device, then copy to host
      d2 = d1;
      out2 = d2;

      // Copy directly to host
      out3 = in;

      // Check that out1, out2, and out3 all match in
      for (int i = 0; i < nx; ++i ) {
         TEST_ASSERT(eq(in[i], out1[i]));
         TEST_ASSERT(eq(in[i], out2[i]));
         TEST_ASSERT(eq(in[i], out3[i]));
      }

      // Copy a slice of d1, check that it is correct
      out4.copySlice(d1, 3);
      for (int i = 0; i < nx/2; ++i ) {
         TEST_ASSERT(eq(in[i+3], out4[i]));
      }
   }
};

TEST_BEGIN(CudaArrayTest)
TEST_ADD(CudaArrayTest, testConstructors)
TEST_ADD(CudaArrayTest, testAllocate)
TEST_ADD(CudaArrayTest, testAssociate)
TEST_ADD(CudaArrayTest, testAssignmentOperators)
TEST_ADD(CudaArrayTest, testConstAssignmentOperators)
TEST_END(CudaArrayTest)

#endif
