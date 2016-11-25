#ifndef BLOCK_DESCRIPTOR_TEST_H
#define BLOCK_DESCRIPTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/BlockDescriptor.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class BlockDescriptorTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      BlockDescriptor v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      BlockDescriptor v;
      std::ifstream in;
      openInputFile("in/BlockDescriptor", in);

      in >> v;
      TEST_ASSERT(v.id() == 5);
      TEST_ASSERT(v.monomerId() == 0);
      TEST_ASSERT(v.vertexId(0) == 3);
      TEST_ASSERT(v.vertexId(1) == 4);
      TEST_ASSERT(eq(v.length(), 2.0));
      // std::cout << v << std::endl ;
   }

};

TEST_BEGIN(BlockDescriptorTest)
TEST_ADD(BlockDescriptorTest, testConstructor)
TEST_ADD(BlockDescriptorTest, testReadWrite)
TEST_END(BlockDescriptorTest)

#endif
