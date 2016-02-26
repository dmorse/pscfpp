#ifndef CHI_INTERACTION_TEST_H
#define CHI_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/ChiInteraction.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class ChiInteractionTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      ChiInteraction v;
      v.setNMonomer(2);
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      ChiInteraction v;
      v.setNMonomer(2);
      std::ifstream in;
      openInputFile("in/ChiInteraction", in);

      v.readParam(in);
      // TEST_ASSERT(v.id() == 5);
      v.writeParam(std::cout);
   }

};

TEST_BEGIN(ChiInteractionTest)
TEST_ADD(ChiInteractionTest, testConstructor)
TEST_ADD(ChiInteractionTest, testReadWrite)
TEST_END(ChiInteractionTest)

#endif
