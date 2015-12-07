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

      SystemStub p;
      p.readParam(in);

      p.writeParam(std::cout);

      #if 0
      for (int i = 0; i < p.nVertex(); ++i) {
         std::cout << p.vertex(i).size() << "\n";
      }

      for (int i = 0; i < p.nSolver(); ++i) {
         std::cout << p.solverId(i)[0] << "  " 
                   << p.solverId(i)[1] << "\n";
      }
      #endif
      
   }

};

TEST_BEGIN(SystemStubTest)
TEST_ADD(SystemStubTest, testConstructor)
TEST_ADD(SystemStubTest, testReadParam)
TEST_END(SystemStubTest)

#endif
