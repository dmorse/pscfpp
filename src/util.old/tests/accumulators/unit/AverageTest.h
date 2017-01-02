#ifndef AVERAGE_TEST_H
#define AVERAGE_TEST_H

#include <test/UnitTest.h>
#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <util/accumulators/Average.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>
#include <util/archives/MemoryCounter.h>
#include <util/archives/BinaryFileIArchive.h>
#include <util/archives/BinaryFileOArchive.h>

using namespace Util;

class AverageTest : public UnitTest
{

   Average accumulator_;

public:

   AverageTest();
   ~AverageTest();

   void setUp();
   void readData();
   void testReadParam();
   void testSample();
   void testSerialize();
   void testSerializeFile();
   void testSaveLoad();

};

AverageTest::AverageTest(){}
AverageTest::~AverageTest(){}

void AverageTest::setUp() {
   std::ifstream paramFile;
   openInputFile("in/Average", paramFile); 
   accumulator_.readParam(paramFile);
   paramFile.close();
}
   
void AverageTest::readData() 
{
   int i, n;
   double x;
   std::ifstream dataFile;
   openInputFile("in/data", dataFile);
   dataFile >> n;
   for (i = 0; i < n; ++i) {
      dataFile >> x;
      accumulator_.sample(x);
   }
   dataFile.close();
}

void AverageTest::testReadParam() 
{
   printMethod(TEST_FUNC);

   printEndl();      
   accumulator_.writeParam(std::cout);
}

void AverageTest::testSample() 
{
   printMethod(TEST_FUNC);

   readData();

   printEndl();      
   accumulator_.output(std::cout);
}

void AverageTest::testSerialize() 
{
   printMethod(TEST_FUNC);
   printEndl();

   readData();

   MemoryOArchive u;
   int size = memorySize(accumulator_);
   u.allocate(size);

   std::cout << size << std::endl;

   u << accumulator_;
   TEST_ASSERT(u.cursor() == u.begin() + size);

   MemoryIArchive v;
   v = u;

   Average clone;
   v >> clone;
   TEST_ASSERT(v.cursor() == u.begin() + size);

   clone.output(std::cout);
}

void AverageTest::testSerializeFile() 
{
   printMethod(TEST_FUNC);
   printEndl();

   readData();

   BinaryFileOArchive u;
   openOutputFile("binary", u.file());
   u << accumulator_;
   u.file().close();

   Average clone;
   BinaryFileIArchive v;
   openInputFile("binary", v.file());
   v >> clone;
   v.file().close();
   
   clone.output(std::cout);
}

void AverageTest::testSaveLoad() 
{
   printMethod(TEST_FUNC);
   printEndl();

   readData();

   BinaryFileOArchive u;
   openOutputFile("binary", u.file());
   accumulator_.save(u);
   u.file().close();

   Average clone;
   BinaryFileIArchive v;
   openInputFile("binary", v.file());
   clone.load(v);
   v.file().close();
 
   clone.writeParam(std::cout);
   clone.output(std::cout);
}

TEST_BEGIN(AverageTest)
TEST_ADD(AverageTest, testReadParam)
TEST_ADD(AverageTest, testSample)
TEST_ADD(AverageTest, testSerialize)
TEST_ADD(AverageTest, testSerializeFile)
TEST_ADD(AverageTest, testSaveLoad)
TEST_END(AverageTest)

#endif
