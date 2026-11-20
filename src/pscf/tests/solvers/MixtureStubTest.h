#ifndef PSCF_MIXTURE_STUB_TEST_H
#define PSCF_MIXTURE_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MixtureStub.h"
#include <pscf/floryHuggins/FhMixture.h>

#include <string>
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
   * A negative block monomerId or a solvent monomerId >= nMonomer was
   * accepted and used as an array index, which crashed the programs.
   */
   void testReadMonomerIdRange()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(!readThrows(mixtureText("0", "0", "0.05", "0.4")));
      TEST_ASSERT(readThrows(mixtureText("-1", "0", "0.05", "0.4")));
      TEST_ASSERT(readThrows(mixtureText("0", "2", "0.05", "0.4")));
   }

   /*
   * A volume fraction outside [0, 1], a solvent size <= 0, or vMonomer
   * <= 0 was accepted, and gave a chemical potential or free energy
   * that is NaN or infinite.
   */
   void testReadValueRange()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(readThrows(mixtureText("0", "0", "0.05", "-0.25", "1.25")));
      TEST_ASSERT(readThrows(mixtureText("0", "0", "0.0", "0.4")));
      TEST_ASSERT(readThrows(mixtureText("0", "0", "0.05", "0.4", "0.6",
                                         "0.0")));
   }

   /*
   * The volume fractions of a closed mixture must add up to 1 within
   * rounding errors (1e-10). A mixture with phi = 0.6 and 0.5 was
   * accepted, and pscf_rpc then converged to a state with a total
   * concentration of 1.1.
   */
   void testReadPhiSum()
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(readThrows(mixtureText("0", "0", "0.05", "0.5")));
      TEST_ASSERT(readThrows(mixtureText("0", "0", "0.05", "0.400001")));

      // Fractions that add up to 1 within a rounding error are kept
      std::istringstream in(mixtureText("0", "0", "0.05", "0.4",
                                        "0.6000000000000001"));
      MixtureStub mix;
      mix.readParam(in);
      TEST_ASSERT(mix.polymer(0).phi() == 0.6000000000000001);
      TEST_ASSERT(mix.solvent(0).phi() == 0.4);
   }

private:

   /*
   * Text of a mixture of a diblock (block monomer ids b0 and 1) and a
   * solvent, with nMonomer = 2.
   */
   std::string mixtureText(std::string b0, std::string solventId,
                           std::string size, std::string phiSolvent,
                           std::string phiPolymer = "0.6",
                           std::string vMonomer = "0.05")
   {
      return "Mixture{ nMonomer 2 monomers[ 1.0 1.0 ] nPolymer 1 "
             "nSolvent 1 Polymer{ type linear nBlock 2 blocks[ " + b0
             + " 2.0 1 3.0 ] phi " + phiPolymer + " } Solvent{ monomerId "
             + solventId + " size " + size + " phi " + phiSolvent
             + " } vMonomer " + vMonomer + " }";
   }

   /*
   * Read a MixtureStub from a string. Return true if this throws.
   */
   bool readThrows(std::string text)
   {
      std::istringstream in(text);
      MixtureStub mix;
      bool thrown = false;
      try {
         mix.readParam(in);
      } catch (Util::Exception&) {
         thrown = true;
      }
      return thrown;
   }

};

TEST_BEGIN(MixtureStubTest)
TEST_ADD(MixtureStubTest, testConstructor)
TEST_ADD(MixtureStubTest, testFhMixtureShortBlocks)
TEST_ADD(MixtureStubTest, testReadParam)
TEST_ADD(MixtureStubTest, testReadMonomerIdRange)
TEST_ADD(MixtureStubTest, testReadValueRange)
TEST_ADD(MixtureStubTest, testReadPhiSum)
TEST_END(MixtureStubTest)

#endif
