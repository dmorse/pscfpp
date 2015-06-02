#ifndef BLOCK_TEST_H
#define BLOCK_TEST

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pfts/chem/Block.h>

#include <fstream>

using namespace Pfts;
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
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(BlockTest)
TEST_ADD(BlockTest, testConstructor)
TEST_ADD(BlockTest, testReadWrite)
TEST_END(BlockTest)

#endif
