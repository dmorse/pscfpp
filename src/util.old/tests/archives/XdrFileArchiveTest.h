#ifndef XDR_FILE_ARCHIVE_TEST_H
#define XDR_FILE_ARCHIVE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/archives/XdrFileOArchive.h>
#include <util/archives/XdrFileIArchive.h>
#include "SerializeTestClass.h"

#include <complex>
#include <fstream>

using namespace Util;

class XdrFileArchiveTest : public UnitTest 
{

private:

public:

   XdrFileArchiveTest()
   {}

   void setUp() 
   {}

   void tearDown() {}
   void testOArchiveConstructor();
   void testWriteRead();

};


void XdrFileArchiveTest::testOArchiveConstructor()
{
   printMethod(TEST_FUNC);
   XdrFileOArchive v("xdr");
   fclose(v.file());
}

void XdrFileArchiveTest::testWriteRead()
{
   printMethod(TEST_FUNC);
   //XdrFileOArchive v("xdr");
   XdrFileOArchive v;
   FILE* fo = openFile("xdr", "wb+");
   v.init(fo);

   // Declare variables
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;
   std::string s1, s2;
   Vector a1, a2;
   IntVector b1, b2;

   // Initialize variables
   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   s1 = "My string has spaces";
   a1[0] =  0.5;
   a1[1] = -1.5;
   a1[2] =  3.5;
   b1[0] =  -5;
   b1[1] =  0;
   b1[2] =  7;
  
   // Write variables to OArchive v
   v << i1;
   v & d1;
   v << c1;
   v << s1;
   v << a1;
   v << b1;
   fclose(v.file());

   // Create IArchive u
   //XdrFileIArchive u("xdr");
   XdrFileIArchive u;
   FILE* fi = openFile("xdr", "rb");
   u.init(fi);

   u >> i2;
   TEST_ASSERT(i1 == i2);

   u & d2;
   TEST_ASSERT(d1 == d2);

   u & c2;
   TEST_ASSERT(c1 == c2);

   u & s2;
   TEST_ASSERT(s1 == s2);

   u & a2;
   TEST_ASSERT(a1 == a2);

   u & b2;
   TEST_ASSERT(b1 == b2);

   fclose(u.file());
}

TEST_BEGIN(XdrFileArchiveTest)
TEST_ADD(XdrFileArchiveTest, testOArchiveConstructor)
TEST_ADD(XdrFileArchiveTest, testWriteRead)
TEST_END(XdrFileArchiveTest)

#endif
