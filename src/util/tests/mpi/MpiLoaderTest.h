#ifndef MPI_LOADER_TEST_H
#define MPI_LOADER_TEST_H

#include <util/mpi/MpiFileIo.h>
#include <util/mpi/MpiLoader.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/space/Vector.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/DMatrix.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>
#include <complex>

using namespace Util;

class MpiLoaderTest : public UnitTest
{

public:

   MpiLoaderTest();
   // void setUp() {}
   // void tearDown() {}
   void testPack();

};

MpiLoaderTest::MpiLoaderTest()
 : UnitTest()
{}

// void setUp() {}
// void tearDown() {}

void MpiLoaderTest::testPack()
{
   printMethod(TEST_FUNC);

   // Declare variables
   int i1, i2;
   double d1, d2;
   std::complex<double> c1, c2;
   std::string s1, s2;
   Vector a1, a2;
   double b1[4];
   double b2[4];
   double m1[3][3]; // Physically 3 x 3, logically 2 x 2
   double m2[3][3]; // Physically 3 x 3, logically 2 x 2
   DArray<double> e1;
   DArray<double> e2;
   e1.allocate(4);
   e2.allocate(4);
   FArray<double, 4> f1;
   FArray<double, 4> f2;
   DMatrix<double> g1;
   DMatrix<double> g2;
   g1.allocate(2, 2);
   g2.allocate(2, 2);

   // Initialize variables
   i1 = 3;
   d1 = 45.0;
   c1 = std::complex<double>(3.0, 4.0);
   s1 = "My string has spaces";
   a1[0] =  1.3;
   a1[1] = -2.3;
   a1[2] =  3.3;

   b1[0] = 9.0;
   b1[1] = 8.0;
   b1[2] = 7.0;
   b1[3] = 6.0;

   m1[0][0] = 13.0;
   m1[0][1] = 14.0;
   m1[1][0] = 15.0;
   m1[1][1] = 16.0;

   e1[0] = 9.0;
   e1[1] = 3.0;
   e1[2] = 7.0;
   e1[3] = 6.0;

   f1[0] = 9.0;
   f1[1] = 3.0;
   f1[2] = 7.0;
   f1[3] = 6.0;
  
   g1(0, 0) = 12.0;
   g1(0, 1) = 14.0;
   g1(1, 0) = 19.0;
   g1(1, 1) = 16.0;

   // Write variable values to OArchive file named "binary"
   if (isIoProcessor()) {
      BinaryFileOArchive  v;
      openOutputFile("binary", v.file());
  
      v << i1;
      v & d1;
      v << s1;
      v << a1;
      v.pack(b1, 4);
      v.pack(m1[0], 2, 2, 3);
      v << e1;
      v << f1;
      v << g1;
      //v << c1;

      v.file().close();
   }

   // Create IArchive u, open file for reading
   BinaryFileIArchive u;
   if (isIoProcessor()) {
      openInputFile("binary", u.file());
   }

   // Construct MpiFileIo and MpiLoader
   MpiFileIo  fileIo_;
   fileIo_.setIoCommunicator(communicator());
   MpiLoader<BinaryFileIArchive> loader(fileIo_, u);

   // Load and test data
   
   loader.load(i2);  // int
   TEST_ASSERT(i1 == i2);

   loader.load(d2);  // double
   TEST_ASSERT(d1 == d2);

   loader.load(s2);   // string
   TEST_ASSERT(s1 == s2);

   loader.load(a2);   // Vector
   TEST_ASSERT(a1 == a2);

   loader.load(b2, 4);    // double C array
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(b1[j] == b2[j]);
   }

   loader.load(m2[0], 2, 2, 3); // double 2D C array 
   int i, j;
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2; ++j) {
         TEST_ASSERT(eq(m1[i][j], m2[i][j]));
      }
   }

   loader.load(e2, 4);  // DArray<double>
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(e1[j] == e2[j]);
   }

   loader.load(f2);  // FArray<double>
   for (int j = 0; j < 4; ++j) {
      TEST_ASSERT(f1[j] == f2[j]);
   }

   loader.load(g2, 2, 2); // DMatrix<double>
   for (i = 0; i < 2; ++i) {
      for (j = 0; j < 2; ++j) {
         TEST_ASSERT(eq(g1(i, j), g2(i, j)));
      }
   }

   //loader.load(u, c2);   // complex
   //TEST_ASSERT(c1 == c2);

}

TEST_BEGIN(MpiLoaderTest)
TEST_ADD(MpiLoaderTest, testPack)
TEST_END(MpiLoaderTest)

#endif
