#ifndef PSCF_MIXTURE_STUB_TEST_H
#define PSCF_MIXTURE_STUB_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include "MixtureStub.h"

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
TEST_ADD(MixtureStubTest, testReadParam)
TEST_ADD(MixtureStubTest, testReadMonomerIdRange)
TEST_END(MixtureStubTest)

#endif
