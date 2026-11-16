#ifndef AMBD_INTERACTION_TEST_H
#define AMBD_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/interaction/Interaction.h>
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

   }

   /*
   * The AM residual uses the inverse of chi and divides by the sum of
   * its elements. A singular chi with 3 monomers aborted through the
   * GSL error handler, chi = 0 with 2 monomers gave a division by a
   * zero determinant, and an invertible chi whose inverse has a zero
   * sum of elements gave a division by zero. update() must throw.
   */
   void testSingularChi()
   {
      printMethod(TEST_FUNC);
      std::ifstream in;

      Interaction v2;
      v2.setNMonomer(2);
      openInputFile("in/Interaction2", in);
      v2.readParam(in);
      in.close();
      v2.setChi(0, 0, 0.0);
      v2.setChi(0, 1, 0.0);
      v2.setChi(1, 1, 0.0);
      TEST_ASSERT(updateThrows(v2));

      // Singular chi: rows 0 and 2 are proportional
      Interaction v3;
      v3.setNMonomer(3);
      openInputFile("in/Interaction3", in);
      v3.readParam(in);
      in.close();
      v3.setChi(0, 0, 0.0);
      v3.setChi(1, 1, 0.0);
      v3.setChi(2, 2, 0.0);
      v3.setChi(0, 1, 30.0);
      v3.setChi(1, 2, 100.0);
      v3.setChi(0, 2, 0.0);
      TEST_ASSERT(updateThrows(v3));

      // Invertible chi (determinant 8000) whose inverse has a zero sum
      // of elements: chi^{-1} (1, 1, 1) = (0.05, -0.1, 0.05)
      v3.setChi(0, 1, 10.0);
      v3.setChi(1, 2, 10.0);
      v3.setChi(0, 2, 40.0);
      TEST_ASSERT(updateThrows(v3));
   }

private:

   /*
   * Return true if AmbdInteraction::update throws an Exception.
   */
   bool updateThrows(Interaction const & v)
   {
      AmbdInteraction u;
      u.setNMonomer(v.nMonomer());
      bool thrown = false;
      try {
         u.update(v);
      } catch (Util::Exception&) {
         thrown = true;
      }
      return thrown;
   }

};

TEST_BEGIN(AmbdInteractionTest)
TEST_ADD(AmbdInteractionTest, testConstructor)
TEST_ADD(AmbdInteractionTest, testReadWrite2)
TEST_ADD(AmbdInteractionTest, testReadWrite3)
TEST_ADD(AmbdInteractionTest, testSingularChi)
TEST_END(AmbdInteractionTest)

#endif
