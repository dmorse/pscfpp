#ifndef PSCF_MIXTURE_STUB_TEST_H
#define PSCF_MIXTURE_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MixtureStub.h"
#include <pscf/floryHuggins/FhMixture.h>

#include <fstream>
#include <sstream>

using namespace Pscf;

class MixtureStubTest : public UnitTest 
{

public:

   void setUp()
   {
      //setVerbose(1);
   }

   void tearDown()
   {
      setVerbose(0);
   }

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      MixtureStub p;
   } 

   void testReadParam() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/Mixture", in);

      MixtureStub sys;
      sys.readParam(in);

      if (verbose() > 0) {
         std::cout << std::endl;
         sys.writeParam(std::cout);
      }
   }

   /*
   * The Flory-Huggins reference of a diblock with block lengths 1e-9
   * and 2e-9 (lengths may be given in any unit): FhMixture::initialize
   * dropped clumps of size up to 1e-8, so the molecule had no clump.
   */
   void testFhMixtureShortBlocks()
   {
      printMethod(TEST_FUNC);

      std::istringstream in(
         "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1 "
         "Polymer{ type linear nBlock 2 blocks[ 0 1.0E-9 1 2.0E-9 ] "
         "phi 1.0 } vMonomer 1.0 }");
      MixtureStub mix;
      mix.readParam(in);

      FhMixture fhMixture;
      fhMixture.initialize(mix);
      TEST_ASSERT(fhMixture.molecule(0).nClump() == 2);
      TEST_ASSERT(fhMixture.molecule(0).clump(0).monomerId() == 0);
      TEST_ASSERT(fhMixture.molecule(0).clump(1).monomerId() == 1);
      TEST_ASSERT(fhMixture.molecule(0).clump(0).size() == 1.0E-9);
      TEST_ASSERT(fhMixture.molecule(0).clump(1).size() == 2.0E-9);
   }

};

TEST_BEGIN(MixtureStubTest)
TEST_ADD(MixtureStubTest, testConstructor)
TEST_ADD(MixtureStubTest, testReadParam)
TEST_ADD(MixtureStubTest, testFhMixtureShortBlocks)
TEST_END(MixtureStubTest)

#endif
