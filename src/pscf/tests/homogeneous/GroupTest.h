#ifndef GROUP_TEST_H
#define GROUP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/homogeneous/Group.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class GroupTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Homogeneous::Group v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Homogeneous::Group v;
      std::ifstream in;
      openInputFile("in/Group", in);

      in >> v;
      TEST_ASSERT(v.monomerId() == 0);
      TEST_ASSERT(eq(v.size(), 2.0));
      // std::cout << v << std::endl ;
   }

};

TEST_BEGIN(GroupTest)
TEST_ADD(GroupTest, testConstructor)
TEST_ADD(GroupTest, testReadWrite)
TEST_END(GroupTest)

#endif
