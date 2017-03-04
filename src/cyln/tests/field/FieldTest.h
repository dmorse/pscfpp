#ifndef CYLN_FIELD_TEST_H
#define CYLN_FIELD_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cyln/field/Field.h>

#if 0
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>
#endif

using namespace Util;
using namespace Pscf::Cyln;

class FieldTest : public UnitTest 
{

public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();
   void testAllocate();
   void testSubscript();
   void testSlice();
 
   #if 0
   void testSerialize1Memory();
   void testSerialize2Memory();
   void testSerialize1File();
   void testSerialize2File();
   #endif

   void fillValues(Field<double>& v);
};


void FieldTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      TEST_ASSERT(v.capacity() == 0 );
      TEST_ASSERT(!v.isAllocated() );
   }
} 

void FieldTest::testAllocate()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      int nr = 10;
      int nz = 10;
      int capacity = nr*nz;

      v.allocate(nr, nz);
      TEST_ASSERT(v.capacity() == capacity);
      TEST_ASSERT(v.nr() == nr);
      TEST_ASSERT(v.nz() == nz);
      TEST_ASSERT(v.isAllocated());
   }
} 

void FieldTest::fillValues(Field<double>& v) 
{
   // Fill field values
   int i, j, k;
   k = 0;
   for (i=0; i < v.nr(); i++) {
      for (j=0; j < v.nz(); j++) {
         v[k] = (i+1)*20.0 + j*1.0;
         ++k;
      }
   }
}

void FieldTest::testSubscript()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      int nr = 10;
      int nz = 10;
      int capacity = nr*nz;
      v.allocate(nr, nz);
      TEST_ASSERT(v.capacity() == capacity);

      fillValues(v);
   
      TEST_ASSERT(v[0]  ==  20.0);
      TEST_ASSERT(v[10] ==  40.0);
      TEST_ASSERT(v[11] ==  41.0);
      TEST_ASSERT(v[12] ==  42.0);
      TEST_ASSERT(v[19] ==  49.0);
      TEST_ASSERT(v[20] ==  60.0);
   }
} 

void FieldTest::testSlice()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      int nr = 10;
      int nz = 10;
      int capacity = nr*nz;
      v.allocate(nr, nz);
      TEST_ASSERT(v.capacity() == capacity);

      fillValues(v);
      Array<double> const & slice = v.slice(1);
      TEST_ASSERT(slice[0] ==  40.0);
      TEST_ASSERT(slice[1] ==  41.0);
      TEST_ASSERT(slice[9] ==  49.0);

   }
}
   
#if 0
void FieldTest::testSerialize1Memory()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
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
   
      Field<double> u;
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

void FieldTest::testSerialize2Memory()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
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
   
      Field<double> u;
   
      // Note: We do not allocate Field<double> u in this test.
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

void FieldTest::testSerialize1File()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
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
   
      Field<double> u;
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

void FieldTest::testSerialize2File()
{
   printMethod(TEST_FUNC);
   {
      Field<double> v;
      v.allocate(3);
      for (int i=0; i < capacity; i++ ) {
         v[i] = (i+1)*10.0;
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
   
      Field<double> u;
   
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

TEST_BEGIN(FieldTest)
TEST_ADD(FieldTest, testConstructor)
TEST_ADD(FieldTest, testAllocate)
TEST_ADD(FieldTest, testSubscript)
TEST_ADD(FieldTest, testSlice)

#if 0
TEST_ADD(FieldTest, testSerialize1Memory)
TEST_ADD(FieldTest, testSerialize2Memory)
TEST_ADD(FieldTest, testSerialize1File)
TEST_ADD(FieldTest, testSerialize2File)
#endif

TEST_END(FieldTest)

#endif
