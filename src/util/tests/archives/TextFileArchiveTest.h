#ifndef TEXT_FILE_ARCHIVE_TEST_H
#define TEXT_FILE_ARCHIVE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/archives/TextFileOArchive.h>
#include <util/archives/TextFileIArchive.h>
#include "SerializeTestClass.h"

#include <complex>
#include <fstream>

using namespace Util;

class TextFileArchiveTest : public UnitTest 
{

public:

   TextFileArchiveTest() {}

   void setUp() {}
   void tearDown() {}

   void testOArchiveConstructor1();
   void testOArchiveConstructor2();
   void testPack();

};

void TextFileArchiveTest::testOArchiveConstructor1()
{
   printMethod(TEST_FUNC);
   TextFileOArchive  v;
   openOutputFile("text", v.file());
} 

void TextFileArchiveTest::testOArchiveConstructor2()
{
   printMethod(TEST_FUNC);
   TextFileOArchive  v("dummy");
   v.file().close();
} 

void TextFileArchiveTest::testPack()
{
   printMethod(TEST_FUNC);
   TextFileOArchive  v;
   openOutputFile("text", v.file());

   // Declare variables
   int i1, i2;
   int u1, u2;
   double d1, d2;
   std::complex<double> c1, c2;
   std::string s1, s2;
   Vector a1, a2;
   SerializeTestClass o1, o2;
   double b1[4];
   double b2[4];
   double m1[3][3]; // intentionally oversized
   double m2[3][3]; // intentionally oversized

   // Initialize variables
   i1 = 3;
   u1 = 7;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   s1 = "My string has spaces";
   a1[0] =  1.3;
   a1[1] = -2.3;
   a1[2] =  3.3;
   o1.i  =  13;
   o1.d  =  26.0;
   b1[0] = 9.0;
   b1[1] = 8.0;
   b1[2] = 7.0;
   b1[3] = 6.0;
   m1[0][0] = 13.0;
   m1[0][1] = 14.0;
   m1[1][0] = 15.0;
   m1[1][1] = 16.0;
  
   // Write variables to OArchive v
   v << i1;
   v & u1;
   v & d1;
   v << c1;
   v << s1;
   v << a1;
   v << o1;
   v.pack(b1, 4);
   v.pack(m1[0], 2, 2, 3);
   v.file().close();

   TextFileIArchive u;
   openInputFile("text", u.file());

   u >> i2;
   TEST_ASSERT(i1 == i2);

   u & u2;
   TEST_ASSERT(u1 == u2);

   u & d2;
   TEST_ASSERT(d1 == d2);

   u & c2;
   TEST_ASSERT(c1 == c2);

   u >> s2;
   TEST_ASSERT(s1 == s2);

   u >> a2;
   TEST_ASSERT(a1 == a2);

   u >> o2;
   TEST_ASSERT(o2.i == 13);
   TEST_ASSERT(o2.d == 26.0);

   u.unpack(b2, 4);
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(b1[j] == b2[j]);
   }

   u.unpack(m2[0], 2, 2, 3);
   int i, j;
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2; ++j) {
         TEST_ASSERT(eq(m1[i][j], m2[i][j]));
      }
   }

}

TEST_BEGIN(TextFileArchiveTest)
TEST_ADD(TextFileArchiveTest, testOArchiveConstructor1)
TEST_ADD(TextFileArchiveTest, testOArchiveConstructor2)
TEST_ADD(TextFileArchiveTest, testPack)
TEST_END(TextFileArchiveTest)

#endif
