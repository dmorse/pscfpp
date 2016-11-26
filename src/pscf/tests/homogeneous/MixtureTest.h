#ifndef PSCF_HOMOGENEOUS_MIXTURE_TEST_H
#define PSCF_HOMOGENEOUS_MIXTURE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/homogeneous/Mixture.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class MixtureTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Homogeneous::Mixture mixture;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Homogeneous::Mixture mixture;
      std::ifstream in;
      openInputFile("in/Mixture", in);

      mixture.readParam(in);
      TEST_ASSERT(mixture.nMolecule() == 1);
      TEST_ASSERT(eq(mixture.molecule(0).size(), 5.0));
      TEST_ASSERT(mixture.molecule(0).nClump() == 2);
      TEST_ASSERT(mixture.molecule(0).clump(0).monomerId() == 0);
      TEST_ASSERT(eq(mixture.molecule(0).clump(0).size(), 2.0));
      TEST_ASSERT(mixture.molecule(0).clump(1).monomerId() == 1);
      TEST_ASSERT(eq(mixture.molecule(0).clump(1).size(), 3.0));
      mixture.writeParam(std::cout) ;
   }

};

TEST_BEGIN(MixtureTest)
TEST_ADD(MixtureTest, testConstructor)
TEST_ADD(MixtureTest, testReadWrite)
TEST_END(MixtureTest)

#endif
