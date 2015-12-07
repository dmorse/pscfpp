#ifndef SYSTEM_STUB_TEST_H
#define SYSTEM_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "SystemStub.h"

#include <fstream>

using namespace Pfts;

class SystemStubTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      SystemStub p;
   } 

   void testReadParam() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/System", in);

      SystemStub sys;
      sys.readParam(in);
      sys.writeParam(std::cout);
   }

   void testCompute() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/System", in);

      SystemStub sys;
      sys.readParam(in);
      sys.compute();
   }

};

TEST_BEGIN(SystemStubTest)
TEST_ADD(SystemStubTest, testConstructor)
TEST_ADD(SystemStubTest, testReadParam)
TEST_ADD(SystemStubTest, testCompute)
TEST_END(SystemStubTest)

#endif
