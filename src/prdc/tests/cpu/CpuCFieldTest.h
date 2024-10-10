#ifndef PRDC_CPU_C_FIELD_TEST_H
#define PRDC_CPU_C_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/CField.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;
using namespace Pscf::Prdc;

class CpuCFieldTest : public UnitTest 
{

private:

   typedef fftw_complex Data;

public:

   void setUp() {}
   void tearDown() {}

   void testConstructor();
   void testAllocateMeshDimension();
   void testSubscript1();
   void testSubscript2();
   void testCopyConstructor();
   void testAssignment();
   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();

};


void CpuCFieldTest::testConstructor()
{

   printMethod(TEST_FUNC);
   {
      Cpu::CField<3> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void CpuCFieldTest::testAllocateMeshDimension()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;
      Cpu::CField<3> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 24);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(v.meshDimensions() == d);
   }
}
 
void CpuCFieldTest::testSubscript1()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = i*10.0 + 3.0;
      }
   
      TEST_ASSERT(eq(v[0][0], 10.0));
      TEST_ASSERT(eq(v[0][1], 3.0));
      TEST_ASSERT(eq(v[2][0], 30.0));
      TEST_ASSERT(eq(v[2][1], 23.0));
   }
}
 
void CpuCFieldTest::testSubscript2()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(capacity == 6);

      Data a, b;
      for (int i=0; i < capacity; i++ ) {
         a[0] = (i+1)*10.0;
         a[1] = i*10.0 + 3.0;
         v[i][0] = a[0];
         v[i][1] = a[1];
      }
 
      b[0] = v[0][0]; 
      b[1] = v[0][1]; 
      TEST_ASSERT(eq(b[0], 10.0));
      TEST_ASSERT(eq(b[1], 3.0));
      TEST_ASSERT(eq(v[2][0], 30.0));
      TEST_ASSERT(eq(v[2][1], 23.0));
   }
}
 
void CpuCFieldTest::testCopyConstructor()
{
   printMethod(TEST_FUNC);

   {
      // Mesh dimensions
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      // Initialize CField
      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(capacity == 6);
      TEST_ASSERT(v.isAllocated() );

      // Preparing test data
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
  
      // Test copy construction of u 
      Cpu::CField<3> u(v);
      TEST_ASSERT(u.capacity() == 6);
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(u.meshDimensions() == v.meshDimensions());
      TEST_ASSERT(eq(v[0][0], 10.0));
      TEST_ASSERT(eq(v[0][1], 3.0));
      TEST_ASSERT(eq(v[2][0], 30.0));
      TEST_ASSERT(eq(v[2][1], 23.0));
      for (int i=0; i < capacity; i++ ) {
         TEST_ASSERT(eq(u[i][0], v[i][0]));
         TEST_ASSERT(eq(u[i][1], v[i][1]));
      }

   }
} 

void CpuCFieldTest::testAssignment()
{
   printMethod(TEST_FUNC);

   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      Cpu::CField<3> u;
      u.allocate(d);
      TEST_ASSERT(u.capacity() == 6);
      TEST_ASSERT(u.isAllocated() );

      // Initialize data   
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
   
      // Assignment operator
      u  = v;
   
      // Test data
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u.capacity() == capacity);
      TEST_ASSERT(eq(v[0][0], 10.0));
      TEST_ASSERT(eq(v[0][1], 3.0));
      TEST_ASSERT(eq(v[2][0], 30.0));
      TEST_ASSERT(eq(v[2][1], 23.0));
      TEST_ASSERT(eq(u[0][0], 10.0));
      TEST_ASSERT(eq(u[0][1], 3.0));
      TEST_ASSERT(eq(u[2][0], 30.0));
      TEST_ASSERT(eq(u[2][1], 23.0));
      for (int i=0; i < capacity; i++ ) {
         TEST_ASSERT(eq(u[i][0], v[i][0]));
         TEST_ASSERT(eq(u[i][1], v[i][1]));
      }
   }
} 

void CpuCFieldTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
      int size = memorySize(v);
     
      int i1 = 13;
      int i2;
   
      MemoryOArchive oArchive;
      oArchive.allocate(size + 12);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
      oArchive << i1;
   
      // Show that v is unchanged by packing
      TEST_ASSERT(eq(v[1][0], 20.0));
      TEST_ASSERT(eq(v[1][1], 13.0));
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::CField<3> u;
      u.allocate(d);
   
      MemoryIArchive iArchive;
      iArchive = oArchive;
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      // Load into u and i2
      iArchive >> u;
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
   
      iArchive >> i2;
      TEST_ASSERT(iArchive.cursor() == iArchive.end());
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
   
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(eq(u[2][0], 30.0));
      TEST_ASSERT(eq(u[2][1], 23.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Release
      iArchive.release();
      TEST_ASSERT(!iArchive.isAllocated());
      TEST_ASSERT(iArchive.begin() == 0);
      TEST_ASSERT(iArchive.cursor() == 0);
      TEST_ASSERT(iArchive.end() == 0);
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size + sizeof(int));
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i][0] = 0.0;
         u[i][1] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      iArchive = oArchive;
      iArchive >> u;
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
   
      iArchive >> i2;
      TEST_ASSERT(iArchive.cursor() == iArchive.end());
      TEST_ASSERT(iArchive.begin() == oArchive.begin());
      TEST_ASSERT(iArchive.end() == oArchive.cursor());
   
      TEST_ASSERT(eq(u[0][0], 10.0));
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(eq(u[2][1], 23.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);

   }

}

void CpuCFieldTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      // Prepare test data
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(eq(v[1][0], 20.0));
      TEST_ASSERT(eq(v[1][1], 13.0));
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::CField<3> u;
   
      // Note: We do not allocate CField u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(u.capacity() == capacity);
   }
}

void CpuCFieldTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      // Initialize test data
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
     
      int i1 = 13;
      int i2;

      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(eq(v[1][0], 20.0));
      TEST_ASSERT(eq(v[1][1], 13.0));
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::CField<3> u;
      u.allocate(d);
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i][0] = 0.0;
         u[i][1] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   }
}

void CpuCFieldTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 1;
      d[1] = 2;
      d[2] = 3;

      Cpu::CField<3> v;
      v.allocate(d);
      int capacity = v.capacity();
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(capacity == 6);

      // Initialize test data
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0;
         v[i][1] = i*10.0 + 3.0;
      }
     
      int i1 = 13;
      int i2;
  
      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(eq(v[1][0], 20.0));
      TEST_ASSERT(eq(v[1][1], 13.0));
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::CField<3> u;
   
      // u.allocate(3); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i][0] = 0.0;
         u[i][1] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1][0], 20.0));
      TEST_ASSERT(eq(u[1][1], 13.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == capacity);
   }
}

TEST_BEGIN(CpuCFieldTest)
TEST_ADD(CpuCFieldTest, testConstructor)
TEST_ADD(CpuCFieldTest, testAllocateMeshDimension)
TEST_ADD(CpuCFieldTest, testSubscript1)
TEST_ADD(CpuCFieldTest, testSubscript2)
TEST_ADD(CpuCFieldTest, testCopyConstructor)
TEST_ADD(CpuCFieldTest, testAssignment)
TEST_ADD(CpuCFieldTest, testSerialize1Memory)
TEST_ADD(CpuCFieldTest, testSerialize2Memory)
TEST_ADD(CpuCFieldTest, testSerialize1File)
TEST_ADD(CpuCFieldTest, testSerialize2File)
TEST_END(CpuCFieldTest)

#endif
