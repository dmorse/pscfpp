#ifndef PSCF_MIXTURE_STUB_TEST_H
#define PSCF_MIXTURE_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MixtureStub.h"

#include <fstream>

using namespace Pscf;

class MixtureStubTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      MixtureStub p;
   } 

   void testReadParam() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/Mixture", in);

      MixtureStub sys;
      sys.readParam(in);
      sys.writeParam(std::cout);
   }

};

TEST_BEGIN(MixtureStubTest)
TEST_ADD(MixtureStubTest, testConstructor)
TEST_ADD(MixtureStubTest, testReadParam)
TEST_END(MixtureStubTest)

#endif
