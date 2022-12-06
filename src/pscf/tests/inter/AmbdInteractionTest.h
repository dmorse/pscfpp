#ifndef AMBD_INTERACTION_TEST_H
#define AMBD_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/inter/Interaction.h>
#include <pscf/iterator/AmbdInteraction.h>
#include <util/param/BracketPolicy.h>


#include <fstream>

using namespace Pscf;
//using namespace Util;

class AmbdInteractionTest : public UnitTest 
{

public:

   void setUp()
   { BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}
  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      AmbdInteraction v;
      v.setNMonomer(2);
   } 

   void testReadWrite2() 
   {
      printMethod(TEST_FUNC);

      Interaction v;
      const int nMonomer = 2;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile("in/Interaction2", in);

      v.readParam(in);

      AmbdInteraction u;
      u.setNMonomer(nMonomer);
      u.update(v);

      TEST_ASSERT(eq(u.chi(0,0), 0.5));
      TEST_ASSERT(eq(u.chi(1,1), 1.5));
      TEST_ASSERT(eq(u.chi(0,1), 2.0));
      TEST_ASSERT(eq(u.chi(1,0), 2.0));
      TEST_ASSERT(eq(u.chi(1,0), u.chi(0,1)));

      int i, j, k;
      double sum;
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            TEST_ASSERT(eq(u.chi(i,j), v.chi(i,j)));
            TEST_ASSERT(eq(u.chiInverse(i,j), v.chiInverse(i,j)));
            TEST_ASSERT(eq(u.p(i,j), v.idemp(i,j)));
            sum = 0.0;
            for (k = 0; k < nMonomer; ++k) {
              sum += u.chi(i,k)*u.chiInverse(k,j);
            }
            if (i == j) {
               TEST_ASSERT(eq(sum, 1.0));
            } else {
               TEST_ASSERT(eq(sum, 0.0));
            }
         }
      }
      TEST_ASSERT(eq(v.sum_inv(), u.sumChiInverse()));

   }

   void testReadWrite3() 
   {
      printMethod(TEST_FUNC);

      Interaction v;
      const int nMonomer = 3;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile("in/Interaction3", in);

      v.readParam(in);

      AmbdInteraction u;
      u.setNMonomer(nMonomer);
      u.update(v);

      TEST_ASSERT(eq(u.chi(1,0), u.chi(0,1)));
      TEST_ASSERT(eq(u.chi(1,2), u.chi(2,1)));
      TEST_ASSERT(eq(u.chi(0,2), u.chi(2,0)));

      int i, j, k;
      double sum;
      for (i = 0; i < nMonomer; ++i) {
         for (j = 0; j < nMonomer; ++j) {
            TEST_ASSERT(eq(u.p(i,j), v.idemp(i,j)));
            sum = 0.0;
            for (k = 0; k < nMonomer; ++k) {
               sum += u.chi(i,k)*u.chiInverse(k,j);
            }
            if (i == j) {
               TEST_ASSERT(eq(sum, 1.0));
            } else {
               TEST_ASSERT(eq(sum, 0.0));
            }
         }
      }
      TEST_ASSERT(eq(v.sum_inv(), u.sumChiInverse()));

   }

};

TEST_BEGIN(AmbdInteractionTest)
TEST_ADD(AmbdInteractionTest, testConstructor)
TEST_ADD(AmbdInteractionTest, testReadWrite2)
TEST_ADD(AmbdInteractionTest, testReadWrite3)
TEST_END(AmbdInteractionTest)

#endif
