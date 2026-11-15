#ifndef PSCF_CUDA_FIELD_TEST_H
#define PSCF_CUDA_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cuda/RField.h>
#include <pscf/field/cuda/CField.h>
#include <pscf/field/cuda/RFieldDft.h>
#include <pscf/field/cuda/resources.h>

#include <pscf/backend/cuda/HostRecvArray.h>
#include <pscf/backend/cuda/HostSendArray.h>
#include <pscf/math/IntVec.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/CField.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <cstring>

using namespace Util;
using namespace Pscf;

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
   void testSerialize();

};

void CudaFieldTest::testConstructors()
{
   printMethod(TEST_FUNC);
   {
      RField<1,CUT> r;
      Pscf::CField<1,CUT> c;
      Pscf::RFieldDft<1,CUT> d;
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
      RField<2,CUT> r;
      Pscf::CField<2,CUT> c;
      Pscf::RFieldDft<2,CUT> d;
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

      // Allocate backend array vd
      RField<3,CUT> vd;
      vd.allocate(d);

      // Initialize host (cpu) data in field vh1
      HostSendArray<cudaReal,CUT> vh1;
      vh1.associate(vd);
      //vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i] = (i+1)*10.0 ;
      }

      // Copy host field vh1 to device field vd
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostRecvArray<cudaReal,CUT> vh2;
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

      // Allocate backend array vd
      Pscf::CField<3,CUT> vd;
      vd.allocate(d);

      // Initialize host (cpu) data in field vh1
      HostSendArray<cudaComplex,CUT> vh1;
      //vh1.allocate(capacity);
      vh1.associate(vd);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostRecvArray<cudaComplex,CUT> vh2;
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
      Pscf::RFieldDft<3,CUT> vd;
      vd.allocate(d);
      int capacity = vd.capacity();

      // Initialize host (cpu) complex data in field vh1
      HostSendArray<cudaComplex,CUT> vh1;
      vh1.associate(vd);
      //vh1.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         vh1[i].x = (i+1)*10.0 ;
         vh1[i].y = i*10.0 + 3.0;
      }

      // Copy host field vh1 to device field vd
      vd = vh1;
      TEST_ASSERT(vd.capacity() == vh1.capacity());

      // Copy device field vd to host host field vh2
      HostRecvArray<cudaComplex,CUT> vh2;
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

void CudaFieldTest::testSerialize()
{
   printMethod(TEST_FUNC);

   IntVec<2> d;
   d[0] = 3;
   d[1] = 5;

   // GPU fields with nonuniform data, and C++ fields with the same data
   RField<2,CUT> r(d);
   Pscf::RFieldDft<2,CUT> k(d);
   Pscf::CField<2,CUT> c(d);
   RField<2,CPT> rc;
   Pscf::RFieldDft<2,CPT> kc;
   Pscf::CField<2,CPT> cc;
   rc.allocate(d);
   kc.allocate(d);
   cc.allocate(d);
   HostSendArray<cudaReal,CUT> r_h(r.capacity());
   HostSendArray<cudaComplex,CUT> k_h(k.capacity()), c_h(c.capacity());
   for (int i = 0; i < r.capacity(); ++i) {
      r_h[i] = 0.5*i - 2.0;
      rc[i] = r_h[i];
   }
   for (int i = 0; i < k.capacity(); ++i) {
      k_h[i].x = 1.0 + i;
      k_h[i].y = -0.25*i;
      kc[i][0] = k_h[i].x;
      kc[i][1] = k_h[i].y;
   }
   for (int i = 0; i < c.capacity(); ++i) {
      c_h[i].x = 2.0 - i;
      c_h[i].y = 0.75*i;
      cc[i][0] = c_h[i].x;
      cc[i][1] = c_h[i].y;
   }
   r = r_h;
   k = k_h;
   c = c_h;

   // Save the GPU fields and the C++ fields, each followed by an int,
   // to two archives, which must be identical
   int i1 = 13;
   MemoryOArchive gpuArchive, cppArchive;
   gpuArchive.allocate(4096);
   cppArchive.allocate(4096);
   gpuArchive << r;
   gpuArchive << k;
   gpuArchive << c;
   gpuArchive << i1;
   cppArchive << rc;
   cppArchive << kc;
   cppArchive << cc;
   cppArchive << i1;
   int size = gpuArchive.cursor() - gpuArchive.begin();
   TEST_ASSERT(cppArchive.cursor() - cppArchive.begin() == size);
   TEST_ASSERT(std::memcmp(gpuArchive.begin(), cppArchive.begin(), size)
               == 0);

   // Load the archive into new GPU fields
   RField<2,CUT> r2;
   Pscf::RFieldDft<2,CUT> k2;
   Pscf::CField<2,CUT> c2;
   int i2 = 0;
   MemoryIArchive iArchive;
   iArchive = gpuArchive;
   iArchive >> r2;
   iArchive >> k2;
   iArchive >> c2;
   iArchive >> i2;
   iArchive.release();
   TEST_ASSERT(i2 == 13);

   TEST_ASSERT(r2.capacity() == r.capacity());
   TEST_ASSERT(r2.meshDimensions() == d);
   TEST_ASSERT(k2.capacity() == k.capacity());
   TEST_ASSERT(k2.meshDimensions() == d);
   TEST_ASSERT(k2.dftDimensions() == k.dftDimensions());
   TEST_ASSERT(c2.capacity() == c.capacity());
   TEST_ASSERT(c2.meshDimensions() == d);

   HostRecvArray<cudaReal,CUT> r2_h(r.capacity());
   HostRecvArray<cudaComplex,CUT> k2_h(k.capacity()), c2_h(c.capacity());
   r2_h = r2;
   k2_h = k2;
   c2_h = c2;
   for (int i = 0; i < r.capacity(); ++i) {
      TEST_ASSERT(r2_h[i] == r_h[i]);
   }
   for (int i = 0; i < k.capacity(); ++i) {
      TEST_ASSERT(k2_h[i].x == k_h[i].x);
      TEST_ASSERT(k2_h[i].y == k_h[i].y);
   }
   for (int i = 0; i < c.capacity(); ++i) {
      TEST_ASSERT(c2_h[i].x == c_h[i].x);
      TEST_ASSERT(c2_h[i].y == c_h[i].y);
   }
}

TEST_BEGIN(CudaFieldTest)
TEST_ADD(CudaFieldTest, testConstructors)
TEST_ADD(CudaFieldTest, testAllocate)
TEST_ADD(CudaFieldTest, testRFieldRoundTrip)
TEST_ADD(CudaFieldTest, testCFieldRoundTrip)
TEST_ADD(CudaFieldTest, testRFieldDftRoundTrip)
TEST_ADD(CudaFieldTest, testSerialize)
TEST_END(CudaFieldTest)

#endif
