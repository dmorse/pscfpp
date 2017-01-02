#ifdef  UTIL_MPI
#ifndef MEMORY_ARCHIVE_TEST_H
#define MEMORY_ARCHIVE_TEST_H

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include "SerializeTestClass.h"

#include <complex>

using namespace Util;

class MpiMemoryArchiveTest : public UnitTest 
{

private:

   const static size_t capacity = 128;

public:

   void setUp() {}
   void tearDown() {}
   void testSendRecv();
   void testSerializeObject();

};


void MpiMemoryArchiveTest::testSendRecv()
{
   printMethod(TEST_FUNC);

   // Declare variables
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;
   std::string s1, s2;
   Vector a1, a2;
   size_t offset;   

   // Initialize variables
   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   s1 = "My string has spaces";
   a1[0] =  1.3;
   a1[1] = -2.3;
   a1[2] =  3.3;
     
   if (mpiRank() == 0) {

      MemoryOArchive  v;
      v.allocate(capacity);
      TEST_ASSERT(v.capacity() == capacity);
      TEST_ASSERT(v.isAllocated() );
      TEST_ASSERT(v.cursor() == v.begin());

      // Write variables to OArchive v
      v << i1;
      offset = sizeof(int);
      TEST_ASSERT(v.cursor() == v.begin() + offset);
   
      v & d1;
      offset += sizeof(double);
      TEST_ASSERT(v.cursor() == v.begin() + offset);
   
      v << c1;
      offset += sizeof(std::complex<double>);
      TEST_ASSERT(v.cursor() == v.begin() + offset);
   
      v << s1;
      offset += memorySize(s1);
      TEST_ASSERT(v.cursor() == v.begin() + offset);
   
      v << a1;
      offset += memorySize(a1);
      TEST_ASSERT(v.cursor() == v.begin() + offset);


      v.send(communicator(), 1);

   }

   if (mpiRank() == 1) {

      // Read data from IArchive u
      MemoryIArchive u;
      u.allocate(capacity);

      u.recv(communicator(), 0);

      TEST_ASSERT(u.cursor() == u.begin());

      u >> i2;
      offset = sizeof(int);
      TEST_ASSERT(u.cursor() == u.begin() + offset);
      TEST_ASSERT(i1 == i2);
   
      u & d2;
      offset += sizeof(double);
      TEST_ASSERT(u.cursor() == u.begin() + offset);
      TEST_ASSERT(d1 == d2);
   
      u & c2;
      offset += sizeof(std::complex<double>);
      TEST_ASSERT(u.cursor() == u.begin() + offset);
      TEST_ASSERT(c1 == c2);
   
      u >> s2;
      offset += memorySize(s2);
      TEST_ASSERT(u.cursor() == u.begin() + offset);
      TEST_ASSERT(s1 == s2);
      //std::cout << std::endl;
      //std::cout << s2 << std::endl;
   
      u >> a2;
      offset += memorySize(a2);
      TEST_ASSERT(u.cursor() == u.begin() + offset);
      TEST_ASSERT(a1 == a2);
  
   }

}
 

#if 0
void MpiMemoryArchiveTest::testSerializeObject()
{
   printMethod(TEST_FUNC);

   MemoryOArchive v;
   MemoryIArchive u;
   SerializeTestClass a, b;

   a.i =13;
   a.d =26.0;

   int size = memorySize(a);
   v.allocate(size);

   v << a;
   TEST_ASSERT(v.cursor() == v.begin() + size);

   u = v;
   u >> b;

   TEST_ASSERT(b.i == 13);
   TEST_ASSERT(b.d == 26.0);

}
#endif

TEST_BEGIN(MpiMemoryArchiveTest)
TEST_ADD(MpiMemoryArchiveTest, testSendRecv)
//TEST_ADD(MpiMemoryArchiveTest, testSerializeObject)
TEST_END(MpiMemoryArchiveTest)

#endif
#endif
