#ifndef POLYMER_TYPE_TEST_H
#define POLYMER_TYPE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/PolymerType.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class PolymerTypeTest : public UnitTest 
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

   void testReadWrite() {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/PolymerType", in);

      PolymerType::Enum b;
      in >> b;
      TEST_ASSERT(b == PolymerType::Branched);
      if (verbose() > 0) {
         printEndl();
         std::cout << b << "  ";
      }
      in >> b;
      TEST_ASSERT(b == PolymerType::Linear);
      if (verbose() > 0) {
         std::cout << b << std::endl ;
      }

      // If uncommented out, this one fails to read "Thingy"
      //in >> b;
   }

};

TEST_BEGIN(PolymerTypeTest)
TEST_ADD(PolymerTypeTest, testReadWrite)
TEST_END(PolymerTypeTest)

#endif
