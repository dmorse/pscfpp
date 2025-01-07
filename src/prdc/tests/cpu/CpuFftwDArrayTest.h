#ifndef PRDC_CPU_FFTW_D_ARRAY_TEST_H
#define PRDC_CPU_FFTW_D_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/FftwDArray.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;
using namespace Pscf; 
using namespace Pscf::Prdc; 

class CpuFftwDArrayTest : public UnitTest 
{

private:

   const static int capacity = 3;

public:

   void setUp(){}

   void tearDown(){}

   void testConstructor();
   void testAllocate();
   void testSubscript();
   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();

};


void CpuFftwDArrayTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void CpuFftwDArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
   }
} 

void CpuFftwDArrayTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0 ;
      }
   
      TEST_ASSERT(v[0] == 10.0);
      TEST_ASSERT(v[2] == 30.0);
   }
} 

void CpuFftwDArrayTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
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
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      Cpu::FftwDArray<double> u;
      u.allocate(3);
   
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
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   
      // Release
      iArchive.release();
      TEST_ASSERT(!iArchive.isAllocated());
      TEST_ASSERT(iArchive.begin() == 0);
      TEST_ASSERT(iArchive.cursor() == 0);
      TEST_ASSERT(iArchive.end() == 0);
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size + sizeof(int));
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
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
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }

}

void CpuFftwDArrayTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      Cpu::FftwDArray<double> u;
   
      // Note: We do not allocate Cpu::FftwDArray<double> u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void CpuFftwDArrayTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;

      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      Cpu::FftwDArray<double> u;
      u.allocate(3);
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void CpuFftwDArrayTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      Cpu::FftwDArray<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
      }
     
      int i1 = 13;
      int i2;
  
      BinaryFileOArchive oArchive;
      openOutputFile("out/binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      Cpu::FftwDArray<double> u;
   
      // u.allocate(3); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
      BinaryFileIArchive iArchive;
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
      iArchive.file().close();
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   
      // Clear values of u and i2
      for (int i=0; i < capacity; i++ ) {
         u[i] = 0.0;
      }
      i2 = 0;
   
      // Reload into u and i2
      openInputFile("out/binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}

TEST_BEGIN(CpuFftwDArrayTest)
TEST_ADD(CpuFftwDArrayTest, testConstructor)
TEST_ADD(CpuFftwDArrayTest, testAllocate)
TEST_ADD(CpuFftwDArrayTest, testSubscript)
TEST_ADD(CpuFftwDArrayTest, testSerialize1Memory)
TEST_ADD(CpuFftwDArrayTest, testSerialize2Memory)
TEST_ADD(CpuFftwDArrayTest, testSerialize1File)
TEST_ADD(CpuFftwDArrayTest, testSerialize2File)
TEST_END(CpuFftwDArrayTest)

#endif
