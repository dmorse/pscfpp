#ifndef PSSP_K_MESH_FIELD_TEST_H
#define PSSP_K_MESH_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/RFieldDFT.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;
using namespace Pscf::Pssp;

class RFieldDFTTest : public UnitTest 
{
private:

   const static int capacity = 3;

   typedef double Data;

public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();
   void testAllocate();
   void testAllocate1();
   void testAllocate3();
   void testSubscript();
   void testAssignment();
   //void testSerialize1Memory();
   //void testSerialize2Memory();
   //void testSerialize1File();
   //void testSerialize2File();

};


void RFieldDFTTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void RFieldDFTTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
   }
} 

void RFieldDFTTest::testAllocate1()
{
   printMethod(TEST_FUNC);
   {
      IntVec<1> d;
      d[0] = 3;
      RFieldDFT<1> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 2);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(d == v.meshDimensions());
   }
}
 
void RFieldDFTTest::testAllocate3()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;
      RFieldDFT<3> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 18);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(d == v.meshDimensions());
   }
}
 
void RFieldDFTTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
      }
   
      TEST_ASSERT(v[0][0] == 10.0);
      TEST_ASSERT(v[0][1] == 10.1);
      TEST_ASSERT(v[1][0] == 20.0);
      TEST_ASSERT(v[1][1] == 20.1);
      TEST_ASSERT(v[2][0] == 30.0);
      TEST_ASSERT(v[2][1] == 30.1);
   }
} 

void RFieldDFTTest::testAssignment()
{
   printMethod(TEST_FUNC);

   {
      RFieldDFT<3> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3);
      TEST_ASSERT(v.isAllocated() );
   
      RFieldDFT<3> u;
      u.allocate(3);
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
   
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
      }
   
      u  = v;
   
      TEST_ASSERT(u.capacity() == 3 );
      TEST_ASSERT(u.isAllocated() );
      TEST_ASSERT(v[0][0] == 10.0);
      TEST_ASSERT(v[0][1] == 10.1);
      TEST_ASSERT(v[1][0] == 20.0);
      TEST_ASSERT(v[1][1] == 20.1);
      TEST_ASSERT(u[0][0] == 10.0);
      TEST_ASSERT(u[0][1] == 10.1);
      TEST_ASSERT(u[1][0] == 20.0);
      TEST_ASSERT(u[1][1] == 20.1);
      TEST_ASSERT(u[2][0] == 30.0);
      TEST_ASSERT(u[2][1] == 30.1);
   }
} 

#if 0
void RFieldDFTTest::testSerialize1Memory()
{ 
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
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
      TEST_ASSERT(v[1][0]==20.0);
      TEST_ASSERT(v[1][1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      RFieldDFT<3> u;
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
   
      TEST_ASSERT(u[0][0] == 10.0);
      TEST_ASSERT(u[0][1] == 10.1);
      TEST_ASSERT(u[1][0] == 20.0);
      TEST_ASSERT(u[1][1] == 20.1);
      TEST_ASSERT(u[2][0] == 30.0);
      TEST_ASSERT(u[2][1] == 30.1);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }

}

void RFieldDFTTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(capacity);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
      }
      int size = memorySize(v);
     
      MemoryOArchive oArchive;
      oArchive.allocate(size);
   
      oArchive << v;
      TEST_ASSERT(oArchive.cursor() == oArchive.begin() + size);
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == capacity);
   
      RFieldDFT<3> u;
   
      // Note: We do not allocate RFieldDFT<3> u in this test.
      // This is the main difference from testSerialize1Memory()
   
      MemoryIArchive iArchive;
   
      iArchive = oArchive;
   
      TEST_ASSERT(iArchive.begin()  == oArchive.begin());
      TEST_ASSERT(iArchive.cursor() == iArchive.begin());
   
      iArchive >> u;
   
      TEST_ASSERT(iArchive.cursor() == iArchive.begin() + size);
      TEST_ASSERT(u[0][0] == 10.0);
      TEST_ASSERT(u[0][1] == 10.1);
      TEST_ASSERT(u[1][0] == 20.0);
      TEST_ASSERT(u[1][1] == 20.1);
      TEST_ASSERT(u[2][0] == 30.0);
      TEST_ASSERT(u[2][1] == 30.1);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void RFieldDFTTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
      }
     
      int i1 = 13;
      int i2;

      BinaryFileOArchive oArchive;
      openOutputFile("binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1]==20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      RFieldDFT<3> u;
      u.allocate(3);
   
      BinaryFileIArchive iArchive;
      openInputFile("binary", iArchive.file());
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
      openInputFile("binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(u[1] == 20.0);
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}

void RFieldDFTTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      RFieldDFT<3> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
      }
     
      int i1 = 13;
      int i2;
  
      BinaryFileOArchive oArchive;
      openOutputFile("binary", oArchive.file());
      oArchive << v;
      oArchive << i1;
      oArchive.file().close();
   
      // Show that v is unchanged by packing
      TEST_ASSERT(v[1] == 20.0);
      TEST_ASSERT(v.capacity() == 3);
   
      RFieldDFT<3> u;
   
      // u.allocate(3); -> 
      // Note: We do not allocate first. This is the difference 
      // from the previous test
   
      BinaryFileIArchive iArchive;
      openInputFile("binary", iArchive.file());
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
      openInputFile("binary", iArchive.file());
      iArchive >> u;
      iArchive >> i2;
   
      TEST_ASSERT(eq(u[1], 20.0));
      TEST_ASSERT(i2 == 13);
      TEST_ASSERT(u.capacity() == 3);
   }
}
#endif

TEST_BEGIN(RFieldDFTTest)
TEST_ADD(RFieldDFTTest, testConstructor)
TEST_ADD(RFieldDFTTest, testAllocate)
TEST_ADD(RFieldDFTTest, testAllocate1)
TEST_ADD(RFieldDFTTest, testAllocate3)
TEST_ADD(RFieldDFTTest, testSubscript)
TEST_ADD(RFieldDFTTest, testAssignment)
//TEST_ADD(RFieldDFTTest, testSerialize1Memory)
//TEST_ADD(RFieldDFTTest, testSerialize2Memory)
//TEST_ADD(RFieldDFTTest, testSerialize1File)
//TEST_ADD(RFieldDFTTest, testSerialize2File)
TEST_END(RFieldDFTTest)

#endif
