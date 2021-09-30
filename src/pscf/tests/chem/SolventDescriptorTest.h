#ifndef SOLVENT_DESCRIPTOR_TEST_H
#define SOLVENT_DESCRIPTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/SolventDescriptor.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class SolventDescriptorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      SolventDescriptor v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      SolventDescriptor v;
      std::ifstream in;
      openInputFile("in/SolventDescriptor", in);

      in >> v;
      TEST_ASSERT(v.monomerId() == 5);
      TEST_ASSERT(eq(v.size(), 3.0));
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(SolventDescriptorTest)
TEST_ADD(SolventDescriptorTest, testConstructor)
TEST_ADD(SolventDescriptorTest, testReadWrite)
TEST_END(SolventDescriptorTest)

#endif
