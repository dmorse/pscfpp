#ifndef CHI_INTERACTION_TEST_H
#define CHI_INTERACTION_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhInteraction.h>
#include <util/param/BracketPolicy.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class FhInteractionTest : public UnitTest 
{

public:

   void setUp()
   {  BracketPolicy::set(BracketPolicy::Optional); }

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      FhInteraction v;
      v.setNMonomer(2);
   } 

   void testReadWrite1() 
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      const int nMonomer = 2;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile("in/FhInteraction", in);

      v.readParam(in);
      if (verbose() > 0){
         printEndl();
         v.writeParam(std::cout);
      }

      TEST_ASSERT(eq(v.chi(0,0), 0.0));
      TEST_ASSERT(eq(v.chi(1,1), 0.0));
      TEST_ASSERT(eq(v.chi(0,1), 2.0));
      TEST_ASSERT(eq(v.chi(1,0), 2.0));
      TEST_ASSERT(eq(v.chi(1,0), v.chi(0,1)));

   }

   void testReadWrite2() 
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      const int nMonomer = 2;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile("in/FhInteraction2", in);

      v.readParam(in);
      if (verbose() > 0){
         printEndl();
         v.writeParam(std::cout);
      }

   }

   void testReadWrite3() 
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      const int nMonomer = 3;
      v.setNMonomer(nMonomer);
      std::ifstream in;
      openInputFile("in/FhInteraction3", in);

      v.readParam(in);
      if (verbose() > 0){
         printEndl();
         v.writeParam(std::cout);
      }

      TEST_ASSERT(eq(v.chi(1,0), v.chi(0,1)));
      TEST_ASSERT(eq(v.chi(1,2), v.chi(2,1)));
      TEST_ASSERT(eq(v.chi(0,2), v.chi(2,0)));

   }

   void testComputeW() 
   {
      printMethod(TEST_FUNC);
      //printEndl();

      FhInteraction v;
      v.setNMonomer(2);
      std::ifstream in;
      openInputFile("in/FhInteraction", in);
      v.readParam(in);

      // Test computeW
      DArray<double> c;
      DArray<double> w;
      c.allocate(2);
      w.allocate(2);
      c[0] = 0.3;
      c[1] = 0.7;
      v. computeW(c, w);
      //TEST_ASSERT(eq(w[0], 0.7));
      //TEST_ASSERT(eq(w[1], 0.3));
      TEST_ASSERT(eq(w[0], 1.4));
      TEST_ASSERT(eq(w[1], 0.6));

      // Test computeC
      w[0] += 0.4;
      w[1] += 0.4;
      double xi;
      v.computeC(w, c, xi);
      TEST_ASSERT(eq(c[0], 0.3));
      TEST_ASSERT(eq(c[1], 0.7));
      TEST_ASSERT(eq(0.4, xi));
   }


   /*
   * A singular chi (one zero chi parameter in a 3-monomer system, as
   * in issue #282) is valid input, but setChi aborted in the GSL error
   * handler. Functions of chi must work, and functions of its inverse
   * must throw an Exception.
   */
   void testSingularChi()
   {
      printMethod(TEST_FUNC);

      FhInteraction v;
      v.setNMonomer(3);
      DMatrix<double> chi;
      chi.allocate(3, 3);
      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            chi(i, j) = 0.0;
         }
      }
      chi(0, 1) = chi(1, 0) = 4.0;
      chi(1, 2) = chi(2, 1) = 6.0;
      v.setChi(chi);

      DArray<double> c, w;
      c.allocate(3);
      w.allocate(3);
      c[0] = 0.25;
      c[1] = 0.5;
      c[2] = 0.25;
      v.computeW(c, w);
      TEST_ASSERT(eq(w[1], 2.5));
      TEST_ASSERT(eq(v.fHelmholtz(c), 1.25));

      bool thrown = false;
      try {
         v.chiInverse(0, 0);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);

      // An inverse computed for an earlier, invertible chi is not reused
      v.setChi(0, 2, 1.0);
      v.setChi(0, 2, 0.0);
      double xi;
      thrown = false;
      try {
         v.computeC(w, c, xi);
      } catch (Exception&) {
         thrown = true;
      }
      TEST_ASSERT(thrown);
   }

};

TEST_BEGIN(FhInteractionTest)
TEST_ADD(FhInteractionTest, testConstructor)
TEST_ADD(FhInteractionTest, testReadWrite1)
TEST_ADD(FhInteractionTest, testReadWrite2)
TEST_ADD(FhInteractionTest, testReadWrite3)
TEST_ADD(FhInteractionTest, testComputeW)
TEST_ADD(FhInteractionTest, testSingularChi)
TEST_END(FhInteractionTest)

#endif
