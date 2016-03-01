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

   void testReadWrite() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      ChiInteraction v;
      v.setNMonomer(2);
      std::ifstream in;
      openInputFile("in/ChiInteraction", in);

      v.readParam(in);
      TEST_ASSERT(eq(v.chi(1,0), v.chi(0,1)));
      // TEST_ASSERT(v.id() == 5);
      v.writeParam(std::cout);
   }

   void testComputeW() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      ChiInteraction v;
      v.setNMonomer(2);
      std::ifstream in;
      openInputFile("in/ChiInteraction", in);
      v.readParam(in);

      DArray<double> c;
      DArray<double> w;
      c.allocate(2);
      w.allocate(2);
      double p = 0.4;
      c[0] = 0.3;
      c[1] = 0.7;
      v. computeW(c, p, w);
      TEST_ASSERT(eq(w[0], 1.1));
      TEST_ASSERT(eq(w[1], 0.7));
   }
};

TEST_BEGIN(ChiInteractionTest)
TEST_ADD(ChiInteractionTest, testConstructor)
TEST_ADD(ChiInteractionTest, testReadWrite)
TEST_ADD(ChiInteractionTest, testComputeW)
TEST_END(ChiInteractionTest)

#endif
