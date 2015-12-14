#ifndef BLOCK_TEST_H
#define BLOCK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/Block.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class BlockTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Block v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Block v;
      std::ifstream in;
      openInputFile("in/Block", in);

      in >> v;
      TEST_ASSERT(v.id() == 5);
      TEST_ASSERT(v.monomerId() == 0);
      TEST_ASSERT(v.vertexId(0) == 3);
      TEST_ASSERT(v.vertexId(1) == 4);
      TEST_ASSERT(eq(v.length(), 2.0));
      // std::cout << v << std::endl ;
   }

};

TEST_BEGIN(BlockTest)
TEST_ADD(BlockTest, testConstructor)
TEST_ADD(BlockTest, testReadWrite)
TEST_END(BlockTest)

#endif
