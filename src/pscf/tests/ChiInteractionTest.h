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

      TEST_ASSERT(eq(v.chiInverse(0,0), 0.0));
      TEST_ASSERT(eq(v.chiInverse(1,1), 0.0));
      TEST_ASSERT(eq(v.chiInverse(0,1), 1.0));
      TEST_ASSERT(eq(v.chiInverse(1,0), 1.0));
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

      // Test computeW
      DArray<double> c;
      DArray<double> w;
      c.allocate(2);
      w.allocate(2);
      c[0] = 0.3;
      c[1] = 0.7;
      v. computeW(c, w);
      TEST_ASSERT(eq(w[0], 0.7));
      TEST_ASSERT(eq(w[1], 0.3));

      // Test computeC
      w[0] += 0.4;
      w[1] += 0.4;
      double xi;
      v.computeC(w, c, xi);
      TEST_ASSERT(eq(c[0], 0.3));
      TEST_ASSERT(eq(c[1], 0.7));
      TEST_ASSERT(eq(0.4, xi));
   }

};

TEST_BEGIN(ChiInteractionTest)
TEST_ADD(ChiInteractionTest, testConstructor)
TEST_ADD(ChiInteractionTest, testReadWrite)
TEST_ADD(ChiInteractionTest, testComputeW)
TEST_END(ChiInteractionTest)

#endif
