#ifndef PSSP_R_FIELD_DFT_TEST_H
#define PSSP_R_FIELD_DFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/RFieldDft.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>

using namespace Util;
using namespace Pscf::Pssp;

class RFieldDftTest : public UnitTest 
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
   void testCopyConst();
   //void testSerialize1Memory();
   //void testSerialize2Memory();
   //void testSerialize1File();
   //void testSerialize2File();

};


void RFieldDftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void RFieldDftTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity );
      TEST_ASSERT(v.isAllocated());
   }
} 

void RFieldDftTest::testAllocate1()
{
   printMethod(TEST_FUNC);
   {
      IntVec<1> d;
      d[0] = 3;
      RFieldDft<1> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 2);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(d == v.meshDimensions());
   }
}
 
void RFieldDftTest::testAllocate3()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 2;
      d[1] = 3;
      d[2] = 4;
      RFieldDft<3> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 18);
      TEST_ASSERT(v.isAllocated());
      TEST_ASSERT(d == v.meshDimensions());
   }
}
 
void RFieldDftTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
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

void RFieldDftTest::testCopyConst()
{
   printMethod(TEST_FUNC);
   {
      IntVec<3> d;
      d[0] = 3;
      d[1] = 3;
      d[2] = 2;
      RFieldDft<3> v;
      v.allocate(d);
      TEST_ASSERT(v.capacity() == 18);
      TEST_ASSERT(v.isAllocated());
      
      for(int i = 0; i < v.capacity(); i++)
      {
         v[i][0] = (i + 1) * 10.0;
         v[i][1] = (i + 1) * 10.0 + 0.1;
      }

      RFieldDft<3> u = v;
      TEST_ASSERT(u.isAllocated());
      TEST_ASSERT(u.capacity() == v.capacity());

      for(int i = 0; i < v.capacity(); i++)
      {
         TEST_ASSERT(u[i][0] == v[i][0]);
         TEST_ASSERT(u[i][1] == v[i][1]);
      }
   }
}

void RFieldDftTest::testAssignment()
{
   printMethod(TEST_FUNC);

   {
      RFieldDft<3> v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == 3);
      TEST_ASSERT(v.isAllocated() );
   
      RFieldDft<3> u;
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
void RFieldDftTest::testSerialize1Memory()
{ 
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
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
   
      RFieldDft<3> u;
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

void RFieldDftTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
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
   
      RFieldDft<3> u;
   
      // Note: We do not allocate RFieldDft<3> u in this test.
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

void RFieldDftTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
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
   
      RFieldDft<3> u;
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

void RFieldDftTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      RFieldDft<3> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i][0] = (i+1)*10.0 ;
         v[i][1] = (i+1)*10.0 + 0.1;
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
   
      RFieldDft<3> u;
   
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
#endif

TEST_BEGIN(RFieldDftTest)
TEST_ADD(RFieldDftTest, testConstructor)
TEST_ADD(RFieldDftTest, testAllocate)
TEST_ADD(RFieldDftTest, testAllocate1)
TEST_ADD(RFieldDftTest, testAllocate3)
TEST_ADD(RFieldDftTest, testSubscript)
TEST_ADD(RFieldDftTest, testAssignment)
TEST_ADD(RFieldDftTest, testCopyConst)
//TEST_ADD(RFieldDftTest, testSerialize1Memory)
//TEST_ADD(RFieldDftTest, testSerialize2Memory)
//TEST_ADD(RFieldDftTest, testSerialize1File)
//TEST_ADD(RFieldDftTest, testSerialize2File)
TEST_END(RFieldDftTest)

#endif
